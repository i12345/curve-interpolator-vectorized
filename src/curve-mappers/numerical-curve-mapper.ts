import { AbstractCurveMapper } from "./abstract-curve-mapper";
import { SplineSegmentOptions } from "../core/interfaces";
import { getGaussianQuadraturePointsAndWeights } from "./gauss";
import { derivativeAtT, derivativeAtT_vectorized, evaluateForT, evaluateForT_vectorized } from "../core/spline-segment";
import { magnitude, magnitude_vectorized } from "../core/math";
import { binarySearch, binarySearch_vectorized, binarySearch_vectorized_specific, clamp } from "../core/utils";
import { IntegerNumberArrayLike, NumberArrayLike, arrayLike } from "../core/array";

export interface CurveLengthCalculationOptions extends SplineSegmentOptions {
  /* Gaussian quadrature weights and abscissae */
  gauss?: [number[], number[]];
  /* from t along arc */
  t0?: number,
  /* to t along arc */
  t1?: number,
}

/**
 * This curve mapper implementation uses a numerical integration method (Gauss Legendre)
 * in order to approximate curve segment lengths. For re-parameterization of the curve
 * function in terms of arc length, a number of precalculated lengths (samples) is used
 * to fit a monotone piecewise cubic function using the approach suggested here:
 * https://stackoverflow.com/questions/35275073/uniform-discretization-of-bezier-curve
 */
export class NumericalCurveMapper<VectorArray extends NumberArrayLike> extends AbstractCurveMapper<VectorArray> {
  _nSamples = 21;
  _gauss: number[][];
  _gauss_vectorized: Float64Array;

  /**
   *
   * @param onInvalidateCache callback function to be invoked when cache is invalidated
   * @param nQuadraturePoints the number of Gauss-Legendre Quadrature points to use for arc length approximation
   * @param nInverseSamples the number of arc length samples to use to fit an inverse function for calculating t from arc length
   */
  constructor(nQuadraturePoints = 24, nInverseSamples = 21, onInvalidateCache?: () => void) {
    super(onInvalidateCache);
    this._gauss = getGaussianQuadraturePointsAndWeights(nQuadraturePoints);
    
    this._gauss_vectorized = new Float64Array(2 * this._gauss.length);
    for (let i = 0; i < this._gauss.length; i++) {
      this._gauss_vectorized[(2 * i) + 0] = this._gauss[i][0];
      this._gauss_vectorized[(2 * i) + 1] = this._gauss[i][1];
    }

    this._nSamples = nInverseSamples;
  }

  /**
   * Clear cache
   */
  override _invalidateCache() {
    super._invalidateCache();
    this._cache['arcLengths'] = null;
    this._cache['samples'] = null;
  }

  get arcLengths(): Float64Array {
    if (!this._cache['arcLengths']) {
      this._cache['arcLengths'] = this.computeArcLengths();
    }
    return this._cache['arcLengths'];
  }

  /**
   * Get samples for inverse function from cache if present, otherwise calculate and put
   * in cache for re-use.
   * @param idx curve segment index
   * @returns Lengths, slopes and coefficients for inverse function
   */
  getSamples(idx: number) : [number[], number[], number[], number[]] {
    if (!this.points) return undefined;
    if (!this._cache['samples']) {
      this._cache['samples'] = new Map<number, [number[], number[], number[], number[]]>();
    }
    if (!this._cache['samples'].has(idx)) {
      const samples = this._nSamples;
      const lengths: number[] = [], slopes: number[] = [];
      const coefficients = this.getCoefficients(idx);
      for (let i = 0; i < samples; ++i) {
        const ti = i / (samples - 1);
        lengths.push(this.computeArcLength(idx, 0.0, ti));
        const dtln = magnitude(evaluateForT(derivativeAtT, ti, coefficients));
        let slope = dtln === 0 ? 0 : 1 / dtln;
        // avoid extreme slopes for near linear curve at the segment endpoints (high tension parameter value)
        if (this.tension > 0.95) {
          slope = clamp(slope, -1, 1);
        }
        slopes.push(slope);
      }

      // Precalculate the cubic interpolant coefficients
      const nCoeff = samples - 1;
      const dis = [];  // degree 3 coefficients
      const cis = [];  // degree 2 coefficients
      let li_prev = lengths[0];
      let tdi_prev = slopes[0];
      const step = 1.0 / nCoeff;

      for (let i = 0; i < nCoeff; ++i) {
        const li = li_prev;
        li_prev = lengths[i+1];
        const lDiff = li_prev - li;
        const tdi = tdi_prev;
        const tdi_next = slopes[i+1];
        tdi_prev = tdi_next;
        const si = step / lDiff;
        const di = (tdi + tdi_next - 2 * si) / (lDiff * lDiff);
        const ci = (3 * si - 2 * tdi - tdi_next) / lDiff;
        dis.push(di);
        cis.push(ci);
      }

      this._cache['samples'].set(idx, [lengths, slopes, cis, dis]);
    }
    return this._cache['samples'].get(idx);
  }

  /**
   * Computes the arc length of a curve segment
   * @param index index of curve segment
   * @param t0 calculate length from t
   * @param t1 calculate length to t
   * @returns arc length between t0 and t1
   */
  computeArcLength(index: number, t0 = 0.0, t1 = 1.0) : number {
    if (t0 === t1) return 0;

    const coefficients = this.getCoefficients(index);
    const z = (t1 - t0) * 0.5;

    let sum = 0;
    for (let i = 0; i < this._gauss.length; i++ ) {
      const [T, C] = this._gauss[i];
      const t = z * T + z + t0;
      const dtln = magnitude(evaluateForT(derivativeAtT, t, coefficients));
      sum += C * dtln;
    }
    return z * sum;
  }

  computeArcLength_vectorized<
      IndexArray extends IntegerNumberArrayLike,
      TArray extends NumberArrayLike,
      LengthArray extends NumberArrayLike
    >(
      index: IndexArray,
      t0: TArray,
      t1: TArray,
      length: LengthArray = <LengthArray><unknown>arrayLike(t0),
      skip?: Uint8Array
    ) {
    const n = index.length;

    const prevSkip = skip;
    const currentSkip = new Uint8Array(n);

    if (prevSkip) {
      for (let i = 0; i < n; i++) {
        if (prevSkip[i] !== 0) {
          currentSkip[i] = prevSkip[i];
          continue;
        }

        if (t0[i] === t1[i]) {
          currentSkip[i] = 1;
          length[i] = 0;
        }
      }
    }
    else {
      for (let i = 0; i < n; i++) {
        if (t0[i] === t1[i]) {
          currentSkip[i] = 1;
          length[i] = 0;
        }
      }
    }

    const gauss_vectorized = this._gauss_vectorized
    const gauss_vectorized_length = gauss_vectorized.length
    const gauss_length = gauss_vectorized_length / 2
    let gauss_vectorized_offset: number

    let t0_i: number
    let z_i: number
    const z = arrayLike(t1);
    let sum: number
    let T: number
    let C: number
    const t = arrayLike(t1);
    const coefficients_indices = new Uint32Array(gauss_length * n)
    const insideSkip = new Uint8Array(gauss_length * n);

    let inside_offset = 0;
    let k: number

    for (let i = 0; i < n; i++) {
      if (currentSkip[i] !== 0) {
        for (let k = 0; k < gauss_length; k++)
          insideSkip[inside_offset++] = 1;

        continue;
      }

      t0_i = t0[i];
      z[i] = z_i = (t0_i + t1[i]) * 0.5;

      for (gauss_vectorized_offset = 0; gauss_vectorized_offset < gauss_vectorized_length; gauss_vectorized_offset += 2) {
        T = gauss_vectorized[gauss_vectorized_offset];
        
        t[inside_offset] = z_i * T + z_i + t0_i;
        coefficients_indices[inside_offset] = i;
        inside_offset++;
      }
    }

    const dtln = magnitude_vectorized(this.dimensionality, evaluateForT_vectorized(derivativeAtT_vectorized, t, coefficients_indices, this.dimensionality, this.getCoefficients_vectorized(), undefined, insideSkip));

    inside_offset = 0;

    for (let i = 0; i < n; i++) {
      if (currentSkip[i] !== 0) {
        inside_offset += gauss_length;
        continue;
      }

      sum = 0;

      for (gauss_vectorized_offset = 1; gauss_vectorized_offset < gauss_vectorized_length; gauss_vectorized_offset += 2) {
        C = gauss_vectorized[gauss_vectorized_offset];
        sum += C * dtln[inside_offset++];
      }

      length[i] = z[i] * sum;
    }

    return length;
  }

  /**
   * Calculate a running sum of arc length for mapping a position on the curve (u)
   * to the position at the corresponding curve segment (t).
   * @returns array with accumulated curve segment arc lengths
   */
  computeArcLengths() : Float64Array {
    if (!this.points) return undefined;
    
    const nPoints = this.closed ? this.points.length : this.points.length - 1;
    
    const lengths = new Float64Array(1 + nPoints);
    lengths[0] = 0;

    let tl = 0;
    for (let i = 0; i < nPoints; i++) {
      const length = this.computeArcLength(i);
      tl += length;
      lengths[i + 1] = tl;
    }

    return lengths;
  }

  /**
   * Calculate t from arc length for a curve segment
   * @param idx segment index
   * @param len length
   * @returns time (t) along curve segment matching the input length
   */
  inverse(idx: number, len: number) : number {
    const nCoeff = this._nSamples - 1;
    const step = 1.0 / nCoeff;
    const [lengths, slopes, cis, dis] = this.getSamples(idx);
    const length = lengths[lengths.length - 1];

    if (len >= length) {
      return 1.0;
    }

    if (len <= 0) {
      return 0.0;
    }

    // Find the cubic segment which has 'len'
    const i = Math.max(0, binarySearch(len, lengths));
    const ti = i * step;
    if (lengths[i] === len) {
      return ti;
    }
    const tdi = slopes[i];
    const di = dis[i];
    const ci = cis[i];
    const ld = len - lengths[i];

    return ((di * ld + ci) * ld + tdi) * ld + ti;
  }

  inverse_vectorized<
      IdxArray extends IntegerNumberArrayLike,
      LenArray extends NumberArrayLike,
      ResultArray extends NumberArrayLike = LenArray
    >(
      idx: IdxArray,
      len: LenArray,
      result: ResultArray = <ResultArray><unknown>arrayLike(len),
      skip?: Uint8Array
    ): ResultArray {
    const n = idx.length;
    const idx_count = <number>Math.max.apply(undefined, idx) + 1;
    const uint32_invalid = 0xFFFFFFFF;

    const nCoeff = this._nSamples - 1;
    const step = 1.0 / nCoeff;
    
    const lengths_map = new Array<NumberArrayLike>(idx_count);
    const slopes_map = new Array<NumberArrayLike>(idx_count);
    const cis_map = new Array<NumberArrayLike>(idx_count);
    const dis_map = new Array<NumberArrayLike>(idx_count);
    const length_map = new Uint32Array(idx_count).fill(uint32_invalid);
    
    //TODO: consider if whole curve [0, idx_max] should be computed
    
    let idx_k: number

    if (skip) {
      for (let k = 0; k < n; k++) {
        if (skip[k] !== 0) continue;

        idx_k = idx[k];

        if (length_map[idx_k] === uint32_invalid) {
          const [lengths, slopes, cis, dis] = this.getSamples(idx_k);
          lengths_map[idx_k] = lengths;
          slopes_map[idx_k] = slopes;
          cis_map[idx_k] = cis;
          dis_map[idx_k] = dis;
          length_map[idx_k] = lengths[lengths.length - 1];
        }
      }
    }
    else {
      for (let k = 0; k < n; k++) {
        idx_k = idx[k];

        if (length_map[idx_k] === uint32_invalid) {
          const [lengths, slopes, cis, dis] = this.getSamples(idx_k);
          lengths_map[idx_k] = lengths;
          slopes_map[idx_k] = slopes;
          cis_map[idx_k] = cis;
          dis_map[idx_k] = dis;
          length_map[idx_k] = lengths[lengths.length - 1];
        }
      }
    }

    const prevSkip = skip;
    const currentSkip = new Uint8Array(n);

    for (let k = 0; k < n; k++) {
      if (prevSkip[k] !== 0) {
        currentSkip[k] = prevSkip[k];
        continue;
      }

      if (len[k] >= length_map[idx[k]]) {
        currentSkip[k] = 1;
        result[k] = 1.0;
      }
      else if (len[k] <= 0) {
        currentSkip[k] = 1;
        result[k] = 0.0;
      }
    }

    const i = binarySearch_vectorized_specific(len, idx, lengths_map, undefined, currentSkip)

    let i_k: number
    let lengths_idx_i: number
    let len_k: number

    let ti: number
    let tdi: number
    let di: number
    let ci: number
    let ld: number

    for (let k = 0; k < n; k++) {
      if (currentSkip[k] !== 0) continue;

      i_k = i[k];
      idx_k = idx[k];
      lengths_idx_i = lengths_map[idx_k][i_k];
      len_k = len[k];

      ti = i_k * step;

      if (lengths_idx_i === len_k)
        result[k] = ti;
      else {
        tdi = slopes_map[idx_k][i_k];
        di = dis_map[idx_k][i_k];
        ci = cis_map[idx_k][i_k];
        ld = len_k - lengths_idx_i;
        result[k] = ((di * ld + ci) * ld + tdi) * ld + ti;
      }
    }

    return result;
  }

  /**
   * Get curve length at u
   * @param u normalized uniform position along the spline curve
   * @returns length in point coordinates
   */
  lengthAt(u: number) : number {
    return u * this.arcLengths[this.arcLengths.length - 1];
  }

  lengthAt_vectorized<UArray extends NumberArrayLike, LengthArray extends NumberArrayLike>(u: UArray, length: LengthArray = <LengthArray><unknown>arrayLike(u), skip?: Uint8Array): LengthArray {
    const n = u.length;
    const arcLength_last = this.arcLengths[this.arcLengths.length - 1];

    if (skip) {
      for (let i = 0; i < n; i++) {
        if (skip[i] !== 0) continue;

        length[i] = u[i] * arcLength_last;
      }
    }
    else {
      for (let i = 0; i < n; i++) {
        length[i] = u[i] * arcLength_last;
      }
    }

    return length;
  }

  /**
   * Maps a uniform time along the curve to non-uniform time (t)
   * @param u normalized uniform position along the spline curve
   * @returns t encoding segment index and local time along curve
   */
  getT(u: number) : number {
    const arcLengths = this.arcLengths;
    const il = arcLengths.length;
    const targetArcLength = u * arcLengths[il - 1];

    const i = binarySearch(targetArcLength, arcLengths);
    const ti = i / (il - 1);
    if (arcLengths[i] === targetArcLength) {
      return ti;
    }

    const len = targetArcLength - arcLengths[i];
    const fraction = this.inverse(i, len);
    return (i + fraction) / (il - 1);
  }

  getT_vectorized<UArray extends NumberArrayLike, TArray extends NumberArrayLike>(u: UArray, t: TArray = <TArray><unknown>arrayLike(u), skip?: Uint8Array): TArray {
    const n = u.length

    const prevSkip = skip;
    const currentSkip = new Uint8Array(n);
    if (prevSkip) currentSkip.set(prevSkip);

    const arcLengths = this.arcLengths;
    const il = arcLengths.length;
    const il_minus_1 = il - 1;
    
    const targetArcLengths = new Float64Array(n);
    const arcLength_last = arcLengths[il_minus_1];
    
    if (prevSkip) {
      for (let idx = 0; idx < n; idx++) {
        if (prevSkip[idx] !== 0) continue;

        targetArcLengths[idx] = u[idx] * arcLength_last;
      }
    }
    else {
      for (let idx = 0; idx < n; idx++) {
        targetArcLengths[idx] = u[idx] * arcLength_last;
      }
    }

    const i = binarySearch_vectorized(targetArcLengths, arcLengths, undefined, prevSkip);

    let i_idx: number

    if (prevSkip) {
      for (let idx = 0; idx < n; idx++) {
        if (prevSkip[idx] !== 0) continue;

        i_idx = i[idx];
        if (arcLengths[i_idx] === targetArcLengths[idx]) {
          t[idx] = i_idx / il_minus_1;
          currentSkip[idx] = 1;
        }
      }
    }
    else {
      for (let idx = 0; idx < n; idx++) {
        i_idx = i[idx];
        if (arcLengths[i_idx] === targetArcLengths[idx]) {
          t[idx] = i_idx / il_minus_1;
          currentSkip[idx] = 1;
        }
      }
    }

    const len = arrayLike(targetArcLengths);

    for (let idx = 0; idx < n; idx++)
      len[idx] = targetArcLengths[idx] - arcLengths[i[idx]];

    const fraction = this.inverse_vectorized(i, len, undefined, currentSkip);

    for (let idx = 0; idx < n; idx++) {
      if (currentSkip[idx] !== 0) continue;
      
      t[idx] = (i[idx] + fraction[idx]) / il_minus_1;
    }

    return t
  }

  /**
   * Maps a non-uniform time along the curve to uniform time (u)
   * @param t non-uniform time along curve
   * @returns uniform time along curve
   */
  getU(t: number) : number {
    if (t === 0) return 0;
    if (t === 1) return 1;

    const arcLengths = this.arcLengths;
    const al = arcLengths.length - 1;
    const totalLength = arcLengths[al];

    // need to de-normalize t to find the matching length
    const tIdx = t * al;

    const subIdx = Math.floor(tIdx);
    const l1 = arcLengths[subIdx];

    if (tIdx === subIdx) return l1 / totalLength;

    const t0 = tIdx - subIdx;
    const fraction = this.computeArcLength(subIdx, 0, t0);

    return (l1 + fraction) / totalLength;
  }

  getU_vectorized<
      TArray extends NumberArrayLike,
      UArray extends NumberArrayLike
    >(
      t: TArray,
      u: UArray = <UArray><unknown>arrayLike(t),
      skip?: Uint8Array
    ): UArray {
    const n = t.length;

    const prevSkip = skip;
    const currentSkip = new Uint8Array(n);

    const arcLengths = this.arcLengths;
    const al = arcLengths.length - 1;
    const totalLength = arcLengths[al];

    let t_i: number

    if (prevSkip) {
      for (let i = 0; i < n; i++) {
        if (prevSkip[i] !== 0) {
          currentSkip[i] = prevSkip[i];
          continue;
        }

        t_i = t[i];

        if (t_i === 0 || t_i === 1) {
          currentSkip[i] = 1;
          u[i] = t_i;
        }
      }
    }
    else {
      for (let i = 0; i < n; i++) {
        t_i = t[i];

        if (t_i === 0 || t_i === 1) {
          currentSkip[i] = 1;
          u[i] = t_i;
        }
      }
    }

    let tIdx: number
    let subIdx_i: number
    const subIdx = new Uint32Array(n);

    let l1_i: number
    const l1 = arrayLike(arcLengths);

    const t0 = arrayLike(t);

    for (let i = 0; i < n; i++) {
      if (currentSkip[i] !== 0) continue;

      t_i = t[i];

      // need to de-normalize t to find the matching length
      tIdx = t_i * al;

      subIdx_i = Math.floor(tIdx);
      l1_i = arcLengths[subIdx_i];

      if (tIdx === subIdx_i) {
        u[i] = l1_i / totalLength;
        currentSkip[i] = 1;
      }
      else {
        l1[i] = l1_i;
        subIdx[i] = subIdx_i;
        t0[i] = tIdx - subIdx_i;
      }
    }

    const fraction = this.computeArcLength_vectorized(subIdx, arrayLike(t0), t0);

    for (let i = 0; i < n; i++) {
      if (currentSkip[i] !== 0) continue;

      u[i] = (l1[i] + fraction[i]) / totalLength;
    }

    return u;
  }
}
