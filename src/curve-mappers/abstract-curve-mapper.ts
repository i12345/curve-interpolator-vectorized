import { NumberArrayLike, arrayLike } from "../core/array";
import { CurveMapper, NumArray4, SegmentFunction, SegmentFunction_vectorized, Vector } from "../core/interfaces";
import { getControlPoints, getSegmentIndexAndT, getSegmentIndexAndT_vectorized } from "../core/spline-curve";
import { calculateCoefficients, evaluateForT, evaluateForT_vectorized } from "../core/spline-segment";


/**
 * The curve mapper's main responsibility is to map between normalized
 * curve position (u) to curve segments and segment position (t). Since
 * it requires access to control points and curve parameters, it also keeps
 * this data along with an internal cache. For this reason, the common
 * functionality has been but into this abstract class definition, so that
 * the mapping specific implementation can be held at a minimum by extending
 * this class.
 */
export abstract class AbstractCurveMapper<VectorArray extends NumberArrayLike> implements CurveMapper<VectorArray> {
  _subDivisions: number;
  _cache: object;
  _points: Vector[];
  _points_vectorized: VectorArray
  _alpha = 0.0;
  _tension = 0.5;
  _closed = false;
  _onInvalidateCache: () => void = null;

  _dimensionality: number
  _vectorized_coefficients: Float64Array = null;

  /**
   * AbstractCurveMapper Constructor
   * @param onInvalidateCache callback function to be invoked when cache needs to be reset
   */
  constructor(onInvalidateCache: () => void = null) {
    this._onInvalidateCache = onInvalidateCache;
    this._cache = {
      arcLengths: null,
      coefficients: null,
    };
  }

  /**
   * Clears cache and invoke callback if provided
   * @returns void
   */
  protected _invalidateCache() : void {
    if (!this.points) return;
    this._cache = {
      arcLengths: null,
      coefficients: null,
    };
    this._vectorized_coefficients = null;
    if (this._onInvalidateCache) this._onInvalidateCache();
  }

  /**
   * Returns the curve length in point coordinates from the global
   * curve position u, where u=1 is the full length of the curve.
   * @param u normalized position on curve (0..1)
   */
  abstract lengthAt(u: number): number;
  abstract lengthAt_vectorized<UArray extends NumberArrayLike, LengthArray extends NumberArrayLike>(u: UArray, length?: LengthArray): LengthArray;
  abstract getT(u: number) : number;
  abstract getT_vectorized<UArray extends NumberArrayLike, TArray extends NumberArrayLike>(u: UArray, t?: TArray): TArray;
  abstract getU(t: number) : number;
  abstract getU_vectorized<TArray extends NumberArrayLike, UArray extends NumberArrayLike>(t: TArray, u?: UArray): UArray;

  /**
   * Curve alpha parameter (0=uniform, 0.5=centripetal, 1=chordal)
   */
  get alpha() { return this._alpha; }
  set alpha(alpha: number) {
    if (Number.isFinite(alpha) && alpha !== this._alpha) {
      this._invalidateCache();
      this._alpha = alpha;
    }
  }

  /**
   * Curve tension (0=Catmull-rom, 1=linear)
   */
  get tension() { return this._tension; }
  set tension(tension: number) {
    if (Number.isFinite(tension) && tension !== this._tension) {
      this._invalidateCache();
      this._tension = tension;
    }
  }

  /**
   * Control points for curve
   */
  get points() { return this._points; }
  set points(points: Vector[]) {
    if (!points || points.length < 2) throw Error('At least 2 control points are required!');
    this._points = points;
    this.recompoutePoints_vectorized();
    this._invalidateCache();
  }

  get points_vectorized() {
    return this._points_vectorized
  }

  private recompoutePoints_vectorized() {
    this._dimensionality = this.points[0].length
    if (!this.points.every(point => point.length === this.dimensionality))
      throw new Error("dimensionality of points must match");
    
    //TODO: insert constructor
    this._points_vectorized = <VectorArray><unknown>new Float64Array(this.points.length * this.dimensionality);
    for (let i = 0; i < this.points.length; i++)
      for (let j = 0; j < this.dimensionality; j++)
        this._points_vectorized[(i * this.dimensionality) + j] = this.points[i][j];
  }

  /**
   * Determines whether the curve should be a closed curve or not
   */
  get closed() { return this._closed; }
  set closed(closed: boolean) {
    closed = !!closed;
    if (this._closed !== closed) {
      this._invalidateCache();
      this._closed = closed;
    }
  }

  get dimensionality() {
    return this._dimensionality
  }

  reset() {
    this._invalidateCache();
  }

  /**
   * Evaluate curve segment function at t
   * @param t time along full curve (encodes segment index and segment t)
   * @param target optional target vector
   * @returns vector
   */
  evaluateForT(func: SegmentFunction, t:number, target?: Vector) : Vector {
    const { index, weight } = getSegmentIndexAndT(t, this.points, this.closed);
    const coefficients = this.getCoefficients(index);
    return evaluateForT(func, weight, coefficients, target);
  }

  evaluateForT_vectorized<TArray extends NumberArrayLike>(
      func: SegmentFunction_vectorized<TArray, VectorArray>,
      t: TArray,
      target: VectorArray = <VectorArray><unknown>arrayLike(t, this.dimensionality),
      skip?: Uint8Array
    ): VectorArray {
    const { index, weight } = getSegmentIndexAndT_vectorized(t, this.dimensionality, this.points_vectorized, this.closed, skip);
    const coefficients = this.getCoefficients_vectorized();
    return evaluateForT_vectorized(func, weight, index, this.dimensionality, coefficients, target, skip);
  }

  /**
   * Get the curve function coefficients at the given segment index. The coefficients
   * are calculated once per segment and put in cache until it is invalidated.
   * @param idx segment index
   * @returns coefficients for the curve function at the given segment index
   */
  getCoefficients(idx: number): NumArray4[] {
    if (!this.points) return undefined;
    if (!this._cache['coefficients']) {
      this._cache['coefficients'] = new Map<number, NumArray4[]>();
    }
    if (!this._cache['coefficients'].has(idx)) {
      const [p0, p1, p2, p3] = getControlPoints(idx, this.points, this.closed);
      const coefficients = calculateCoefficients(p0, p1, p2, p3, { tension: this.tension, alpha: this.alpha });
      this._cache['coefficients'].set(idx, coefficients);
    }
    return this._cache['coefficients'].get(idx);
  }

  /**
   * Get the curve function coefficients. The coefficients
   * are calculated for all segments and put in cache until it is invalidated.
   * @returns coefficients packed by [segment index, dimension index, coefficient index]
   */
  getCoefficients_vectorized() {
    if (this._vectorized_coefficients !== null)
      return this._vectorized_coefficients;
    
    const segments = this.closed ? this.points.length : this.points.length - 1;
    
    const vectorized_coefficients = this._vectorized_coefficients = new Float64Array(segments * this.dimensionality * 4)
    let vectorized_coefficients_i = 0

    for (let idx = 0; idx < segments; idx++) {
      const coefficients = this.getCoefficients(idx)
      for (let i = 0; i < this.dimensionality; i++)
        for (let coefficient = 0; coefficient < 4; coefficient++)
          vectorized_coefficients[vectorized_coefficients_i++] = coefficients[i][coefficient]
    }

    return vectorized_coefficients
  }
}
