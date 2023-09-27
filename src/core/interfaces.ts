import { IntegerNumberArrayLike, NumberArrayLike } from "./array";

/**
 * Array of four number items
 */
export type NumArray4 = [number, number, number, number];

/**
 * Either a number array or an object implementing the VectorType interface
 */
export type Vector = (NumberArrayLike | VectorType);

export type SegmentFunction = (t: number, coefficients: NumArray4) => number;
export type SegmentFunction_vectorized<TArray extends NumberArrayLike, VectorArray extends NumberArrayLike> = (t: TArray, coefficients_indices: IntegerNumberArrayLike, dimensionality: number, coefficients_vectorized: Float64Array, results?: VectorArray, skip?: Uint8Array) => VectorArray
export interface CurveMapper<VectorArray extends NumberArrayLike = Float64Array> {
  alpha: number,
  tension: number,
  points: Vector[],
  points_vectorized: VectorArray,
  closed: boolean,

  evaluateForT: (func:SegmentFunction, t:number, target?:VectorType) => Vector,
  evaluateForT_vectorized: <TArray extends NumberArrayLike>(func:SegmentFunction_vectorized<TArray, VectorArray>, t:TArray, target?: VectorArray, skip?: Uint8Array) => VectorArray,
  lengthAt: (u: number) => number,
  lengthAt_vectorized: <UArray extends NumberArrayLike, LengthArray extends NumberArrayLike>(u: UArray, length?: LengthArray, skip?: Uint8Array) => LengthArray,
  getT: (u: number) => number,
  getT_vectorized: <UArray extends NumberArrayLike, TArray extends NumberArrayLike>(u: UArray, t?: TArray, skip?: Uint8Array) => TArray,
  getU: (t: number) => number,
  getU_vectorized: <TArray extends NumberArrayLike, UArray extends NumberArrayLike>(t: TArray, u?: UArray, skip?: Uint8Array) => UArray,
  getCoefficients: (idx: number) => NumArray4[],
  getCoefficients_vectorized: () => Float64Array,
  reset: () => void,
}

/**
 * Any objects that supports indexing values by number may be used as input or return types.
 * See the Point class for an example.
 */
export interface VectorType {
  0: number,
  1: number,
  2?: number,
  3?: number,
  x?: number,
  y?: number,
  z?: number,
  w?: number,
  length: number,
}

export interface CurveParameters {
  /* curve tension (0 = Catmull-Rom curve, 1 = linear curve) */
  tension?: number,
  /* curve velocity vector modifier (0 = uniform, 0.5 = centripetal, 1 = chordal */
  alpha?: number,
}

/**
 * Options required to perform calculations on a curve segment.
 */
export interface SplineSegmentOptions extends CurveParameters {
  knotSequence?: NumArray4,
  target?: Vector,
}

/**
 * Spline Curve characteristics
 */
export interface SplineCurveOptions extends CurveParameters {
  /* flag to set if the curve should be closed or not */
  closed?: boolean,
}

/**
 * Used by the valuesLookup function to set axis, tension etc.
 */
export interface LookupOptions extends SplineCurveOptions {
  axis?: number,
  margin?: number,
  max?: number,
  processRefAxis?: boolean,
}

/**
 * Used by the positions lookup function
 */
export interface PositionLookupOptions extends SplineCurveOptions {
  axis?: number,
  margin?: number,
  max?: number,
}

/**
 * Bounding box interface
 */
export interface BBox {
  min: Vector,
  max: Vector,
}

/**
 * Options to control calculation of bounding box
 */
export interface BBoxOptions extends SplineCurveOptions{
  from?: number,
  to?: number,
}
