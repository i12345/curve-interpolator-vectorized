import { AbstractCurveMapper } from "./abstract-curve-mapper";
import { Vector } from "../core/interfaces";
import { distance, distance_vectorized } from "../core/math";
import { binarySearch, binarySearch_vectorized } from "../core/utils";
import { valueAtT, valueAtT_vectorized } from "../core/spline-segment";
import { NumberArrayLike, arrayLike } from "../core/array";

/**
 * Approximate spline curve by subdividing it into smaller linear
 * line segments. Used to approximate length and mapping between
 * uniform (u) and non-uniform (t) time along curve.
 */
export class SegmentedCurveMapper<VectorArray extends NumberArrayLike> extends AbstractCurveMapper<VectorArray> {

  _subDivisions: number;

  /**
   *
   * @param subDivisions number of sub divisions to use
   * @param onInvalidateCache callback function to be invoked when cache is invalidated
   */
  constructor(subDivisions = 300, onInvalidateCache: () => void = null) {
    super(onInvalidateCache);
    this._subDivisions = subDivisions;
  }

  get arcLengths(): Float64Array {
    if (!this._cache['arcLengths']) {
      this._cache['arcLengths'] = this.computeArcLengths();
    }
    return this._cache['arcLengths'];
  }

  /**
   * Clear cache
   */
  override _invalidateCache() {
    super._invalidateCache();
    this._cache['arcLengths'] = null;
  }

  /**
   * Break curve into segments and return the curve length at each segment index.
   * Used for mapping between t and u along the curve.
   */
  computeArcLengths() {
    const lengths = new Float64Array(1 + this._subDivisions);
    let current: Vector, last = this.evaluateForT(valueAtT, 0);
    let sum = 0;

    lengths[0] = 0;

    for (let p = 1; p <= this._subDivisions; p++) {
      current = this.evaluateForT(valueAtT, p / this._subDivisions);
      sum += distance(current, last);
      lengths[p] = sum;
      last = current;
    }
    return lengths;
  }

  /**
   * Get curve length at u
   * @param u normalized uniform position along the spline curve
   * @returns length in point coordinates
   */
  lengthAt(u: number) {
    const arcLengths = this.arcLengths;
    return u * arcLengths[arcLengths.length - 1];
  }

  lengthAt_vectorized<UArray extends NumberArrayLike, LengthArray extends NumberArrayLike>(u: UArray, length: LengthArray = <LengthArray><unknown>arrayLike(u), skip?: Uint8Array): LengthArray {
    const n = u.length
    const arcLengths = this.arcLengths;
    const arcLength_last = arcLengths[arcLengths.length - 1];

    if (skip) {
      for (let i = 0; i < n; i++) {
        if (skip[i] !== 0) continue;

        length[i] = u[i] * arcLength_last;
      }
    }
    else {
      for (let i = 0; i < n; i++)
        length[i] = u[i] * arcLength_last;
    }
    
    return length;
  }

  /**
   * Maps a uniform time along the curve to non-uniform time (t)
   * @param u normalized uniform position along the spline curve
   * @returns t encoding segment index and local time along curve
   */
  getT(u: number) {
    const arcLengths = this.arcLengths;
    const il = arcLengths.length;
    const targetArcLength = u * arcLengths[il - 1];

    const i = binarySearch(targetArcLength, arcLengths);
    if (arcLengths[i] === targetArcLength) {
      return i / (il - 1);
    }

    // we could get finer grain at lengths, or use simple interpolation between two points
    const lengthBefore = arcLengths[i];
    const lengthAfter = arcLengths[i + 1];
    const segmentLength = lengthAfter - lengthBefore;

    // determine where we are between the 'before' and 'after' points
    const segmentFraction = (targetArcLength - lengthBefore) / segmentLength;

    // add that fractional amount to t
    return (i + segmentFraction) / (il - 1);
  }

  getT_vectorized<UArray extends NumberArrayLike, TArray extends NumberArrayLike>(u: UArray, t: TArray = <TArray><unknown>arrayLike(u), skip?: Uint8Array): TArray {
    const n = u.length;
    const arcLengths = this.arcLengths;
    const il = arcLengths.length;
    const il_minus_1 = il - 1
    const arcLengths_last = arcLengths[il_minus_1];
    const targetArcLength = t;

    if (skip) {
      for (let i = 0; i < n; i++) {
        if (skip[i] !== 0) continue;

        targetArcLength[i] = u[i] * arcLengths_last;
      }
    }
    else
      for (let i = 0; i < n; i++)
        targetArcLength[i] = u[i] * arcLengths_last;

    const i = binarySearch_vectorized(targetArcLength, arcLengths, undefined, skip);

    let i_k: number
    let arcLength_i: number
    let targetArcLength_k: number
    let lengthBefore: number
    let lengthAfter: number
    let segmentLength: number
    let segmentFraction: number

    if (skip) {
      for (let k = 0; k < n; k++) {
        if (skip[k] !== 0) continue;

        i_k = i[k];
        arcLength_i = arcLengths[i_k];
        targetArcLength_k = targetArcLength[k];
      
        if (arcLength_i === targetArcLength_k)
          t[k] = i_k / il_minus_1;
        else {
          lengthBefore = arcLengths[i_k];
          lengthAfter = arcLengths[i_k + 1];
          segmentLength = lengthAfter - lengthBefore;
          segmentFraction = (targetArcLength_k - lengthBefore) / segmentLength;
          t[k] = (i_k + segmentFraction) / il_minus_1;
        }
      }
    }
    else {
      for (let k = 0; k < n; k++) {
        i_k = i[k];
        arcLength_i = arcLengths[i_k];
        targetArcLength_k = targetArcLength[k];
      
        if (arcLength_i === targetArcLength_k)
          t[k] = i_k / il_minus_1;
        else {
          lengthBefore = arcLengths[i_k];
          lengthAfter = arcLengths[i_k + 1];
          segmentLength = lengthAfter - lengthBefore;
          segmentFraction = (targetArcLength_k - lengthBefore) / segmentLength;
          t[k] = (i_k + segmentFraction) / il_minus_1;
        }
      }
    }

    return t
  }

  /**
   * Maps a non-uniform time along the curve to uniform time (u)
   * @param t non-uniform time along curve
   * @returns uniform time along curve
   */
  getU(t: number) {
    if (t === 0) return 0;
    if (t === 1) return 1;

    const arcLengths = this.arcLengths;
    const al = arcLengths.length - 1;
    const totalLength = arcLengths[al];

    // need to denormalize t to find the matching length
    const tIdx = t * al;

    const subIdx = Math.floor(tIdx);
    const l1 = arcLengths[subIdx];

    if (tIdx === subIdx) return l1 / totalLength;

    // measure the length between t0 at subIdx and t
    const t0 = subIdx / al;
    const p0 = this.evaluateForT(valueAtT, t0);
    const p1 = this.evaluateForT(valueAtT, t);
    const l = l1 + distance(p0, p1);

    //const l2 = arcLengths[subIdx + 1];
    //const l = l1 + (tIdx - subIdx) * (l2 - l1);

    return l / totalLength;
  }

  getU_vectorized<TArray extends NumberArrayLike, UArray extends NumberArrayLike>(t: TArray, u: UArray = <UArray><unknown>arrayLike(t), skip?: Uint8Array): UArray {
    const n = t.length;
    const prevSkip = skip
    const currentSkip = new Uint8Array(n).fill(0)

    const arcLengths = this.arcLengths;
    const al = arcLengths.length - 1;
    const totalLength = arcLengths[al];

    let t_i: number
    let tIdx: number
    let subIdx: number
    let l1_i: number

    const l1 = new Float64Array(n)
    const t0 = arrayLike(t);

    if (prevSkip) {
      for (let i = 0; i < n; i++) {
        if (prevSkip[i] === 1) {
          currentSkip[i] = 1;
          continue;
        }

        t_i = t[i];

        if (t_i === 0 || t_i === 1) {
          u[i] = t_i;
          currentSkip[i] = 1;
          continue;
        }
    
        tIdx = t_i * al;
        subIdx = Math.floor(tIdx);
        l1_i = arcLengths[subIdx];

        if (tIdx === subIdx) {
          u[i] = l1_i / totalLength;
          currentSkip[i] = 1;
          continue;
        }

        l1[i] = l1_i
        t0[i] = subIdx / al;
      }
    }
    else {
      for (let i = 0; i < n; i++) {
        t_i = t[i];

        if (t_i === 0 || t_i === 1) {
          u[i] = t_i;
          currentSkip[i] = 1;
          continue;
        }
    
        tIdx = t_i * al;
        subIdx = Math.floor(tIdx);
        l1_i = arcLengths[subIdx];

        if (tIdx === subIdx) {
          u[i] = l1_i / totalLength;
          currentSkip[i] = 1;
          continue;
        }

        l1[i] = l1_i
        t0[i] = subIdx / al;
      }
    }

    const p0 = this.evaluateForT_vectorized(valueAtT_vectorized, t0, undefined, currentSkip)
    const p1 = this.evaluateForT_vectorized(valueAtT_vectorized, t, undefined, currentSkip)
    const distances = distance_vectorized(this.dimensionality, p0, p1, u, currentSkip);

    for (let i = 0; i < n; i++) {
      if (currentSkip[i] === 1) continue;

      u[i] = (l1[i] + distances[i]) / totalLength;
    }

    return u;
  }
}
