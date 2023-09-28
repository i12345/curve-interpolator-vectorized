import { IntegerNumberArrayLike, NumberArrayLike, arrayLike } from "./array";
import { NumArray4, Vector, CurveParameters, SegmentFunction, SegmentFunction_vectorized } from "./interfaces";
import { EPS, getCubicRoots, sumOfSquares } from "./math";
import { clamp } from "./utils";



/**
 * This function will calculate the knot sequence, based on a given value for alpha, for a set of
 * control points for a curve segment. It is used to calculate the velocity vectors, which
 * determines the curvature of the segment. Alpha=0.5 produces a centripetal curve, while
 * alpha=1 produces a chordal curve.
 * @param p0 First control point
 * @param p1 Second control point
 * @param p2 Third control point
 * @param p3 Fourth control point
 * @param alpha alpha value
 * @returns calculated knot sequence to use for curve velocity vector calculations
 */
export function calcKnotSequence(p0 : Vector, p1: Vector, p2: Vector, p3: Vector, alpha = 0) : NumArray4 {
  if (alpha === 0) return [0, 1, 2, 3];

  const deltaT = (u: Vector, v: Vector) : number => Math.pow(sumOfSquares(u, v), 0.5 * alpha);

  const t1 = deltaT(p1, p0);
  const t2 = deltaT(p2, p1) + t1;
  const t3 = deltaT(p3, p2) + t2;

  return [0, t1, t2, t3];
}

/**
 * Calculate coefficients for a curve segment with specified parameters
 * @param p0 control point 1
 * @param p1 control point 2
 * @param p2 control point 3
 * @param p3 control point 4
 * @param options curve parameters
 * @returns coefficients for curve function
 */
export function calculateCoefficients(p0: Vector, p1:Vector, p2:Vector, p3:Vector, options: CurveParameters) : NumArray4[] {
  const tension = Number.isFinite(options.tension) ? options.tension : 0.5;
  const alpha = Number.isFinite(options.alpha) ? options.alpha : null;
  const knotSequence = alpha > 0 ? calcKnotSequence(p0, p1, p2, p3, alpha) : null;
  const coefficientsList = new Array<NumArray4>(p0.length);

  for (let k = 0; k < p0.length; k++) {
    let u = 0, v = 0;
    const v0 = p0[k], v1 = p1[k], v2 = p2[k], v3 = p3[k];
    if (!knotSequence) {
      u = (1 - tension) * (v2 - v0) * 0.5;
      v = (1 - tension) * (v3 - v1) * 0.5;
    } else {
      const [t0, t1, t2, t3] = knotSequence;
      if (t1 - t2 !== 0) {
        if (t0 - t1 !== 0 && t0 - t2 !== 0) {
          u = (1 - tension) * (t2 - t1) * ((v0 - v1) / (t0 - t1) - (v0 - v2) / (t0 - t2) + (v1 - v2) / (t1 - t2));
        }
        if (t1 - t3 !== 0 && t2 - t3 !== 0) {
          v = (1 - tension) * (t2 - t1) * ((v1 - v2) / (t1 - t2) - (v1 - v3) / (t1 - t3) + (v2 - v3) / (t2 - t3));
        }
      }
    }

    const a = (2 * v1 - 2 * v2 + u + v);
    const b = (-3 * v1 + 3 * v2 - 2 * u - v);
    const c = u;
    const d = v1;
    coefficientsList[k] = [a, b, c, d];
  }
  return coefficientsList;
}

/**
 * Calculates vector component for a point along the curve segment at time t
 * @param t time along curve segment
 * @param coefficients coefficients for curve function
 * @returns curve value
 */
export function valueAtT(t: number, coefficients: NumArray4) : number {
  const t2 = t * t;
  const t3 = t * t2;
  const [a, b, c, d] = coefficients;
  return a * t3 + b * t2 + c * t + d;
}

export function valueAtT_vectorized<
    TArray extends NumberArrayLike,
    VectorArray extends NumberArrayLike,
  >(
    t: TArray,
    coefficient_indices: IntegerNumberArrayLike,
    dimensionality: number,
    vectorized_coefficients: Float64Array,
    results: VectorArray = <VectorArray><unknown>arrayLike(t, dimensionality),
    skip?: Uint8Array
): VectorArray {
  const n = t.length
  
  let t_i: number
  let t2: number
  let t3: number

  let coefficients_offset: number
  let a: number
  let b: number
  let c: number
  let d: number
  let dimension: number
  let results_offset = 0
  const dimensionality_times_4 = 4 * dimensionality

  if (skip) {
    for (let i = 0; i < n; i++) {
      if (skip[i] !== 0) {
        results_offset += dimensionality;
        continue;
      }

      t_i = t[i];
      t2 = t_i * t_i;
      t3 = t_i * t2;
      coefficients_offset = dimensionality_times_4 * coefficient_indices[i];
      for (dimension = 0; dimension < dimensionality; dimension++) {
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset++];
        results[results_offset++] = a * t3 + b * t2 + c * t_i + d;
      }
    }
  }
  else {
    for (let i = 0; i < n; i++) {
      t_i = t[i];
      t2 = t_i * t_i;
      t3 = t_i * t2;
      coefficients_offset = dimensionality_times_4 * coefficient_indices[i];
      for (dimension = 0; dimension < dimensionality; dimension++) {
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset++];
        results[results_offset++] = a * t3 + b * t2 + c * t_i + d;
      }
    }
  }

  return results;
}

export const valueAtT_vectorized_custom = <
    TArray extends NumberArrayLike,
    VectorArray extends NumberArrayLike,
  >(result_dimensions: IntegerNumberArrayLike) => (
      t: TArray,
      coefficient_indices: IntegerNumberArrayLike,
      dimensionality: number,
      vectorized_coefficients: Float64Array,
      results: VectorArray = <VectorArray><unknown>arrayLike(t, dimensionality),
      skip?: Uint8Array
  ): VectorArray => {
  const n = t.length
  
  let t_i: number
  let t2: number
  let t3: number

  let coefficients_offset_base: number
  let coefficients_offset: number
  let a: number
  let b: number
  let c: number
  let d: number
  let dimension_i: number
  let dimension: number
  let results_offset = 0
  const dimensionality_times_4 = 4 * dimensionality
  
  const result_dimensionality = result_dimensions.length;
  
  if (skip) {
    for (let i = 0; i < n; i++) {
      if (skip[i] !== 0) {
        results_offset += result_dimensionality
        continue;
      }

      t_i = t[i];
      t2 = t_i * t_i;
      t3 = t_i * t2;
      coefficients_offset_base = dimensionality_times_4 * coefficient_indices[i];
      for (dimension_i = 0; dimension_i < result_dimensionality; dimension_i++) {
        dimension = result_dimensions[dimension_i];
        coefficients_offset = coefficients_offset_base + (4 * dimension);
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset];
        results[results_offset++] = a * t3 + b * t2 + c * t_i + d;
      }
    }
  }
  else {
    for (let i = 0; i < n; i++) {
      t_i = t[i];
      t2 = t_i * t_i;
      t3 = t_i * t2;
      coefficients_offset_base = dimensionality_times_4 * coefficient_indices[i];
      for (dimension_i = 0; dimension_i < result_dimensionality; dimension_i++) {
        dimension = result_dimensions[dimension_i];
        coefficients_offset = coefficients_offset_base + (4 * dimension);
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset];
        results[results_offset++] = a * t3 + b * t2 + c * t_i + d;
      }
    }
  }

  return results;
}

export const axisValueAtT_vectorized = <
      TArray extends NumberArrayLike,
      VectorArray extends NumberArrayLike,
    >(axis: number) => (
      t: TArray,
      coefficient_indices: IntegerNumberArrayLike,
      dimensionality: number,
      vectorized_coefficients: Float64Array,
      results: VectorArray = <VectorArray><unknown>arrayLike(t),
      skip?: Uint8Array
    ): VectorArray => {
    const n = t.length
    
    let t_i: number
    let t2: number
    let t3: number

    let coefficients_offset: number
    let a: number
    let b: number
    let c: number
    let d: number
    const dimensionality_times_4 = 4 * dimensionality;
    const axis_times_4 = 4 * axis;
    
    if (skip) {
      for (let i = 0; i < n; i++) {
        if (skip[i] !== 0) continue;

        t_i = t[i];
        t2 = t_i * t_i;
        t3 = t_i * t2;
        coefficients_offset = (dimensionality_times_4 * coefficient_indices[i]) + axis_times_4;
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset++];
        results[i] = a * t3 + b * t2 + c * t_i + d;
      }
    }
    else {
      for (let i = 0; i < n; i++) {
        t_i = t[i];
        t2 = t_i * t_i;
        t3 = t_i * t2;
        coefficients_offset = (dimensionality_times_4 * coefficient_indices[i]) + axis_times_4;
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset++];
        results[i] = a * t3 + b * t2 + c * t_i + d;
      }
    }

    return results;
  }


export const axisValueAtT_vectorized_custom = <
      TArray extends NumberArrayLike,
      VectorArray extends NumberArrayLike,
    >(axis: number, stride = 1, offset = 0) => (
      t: TArray,
      coefficient_indices: IntegerNumberArrayLike,
      dimensionality: number,
      vectorized_coefficients: Float64Array,
      results: VectorArray = <VectorArray><unknown>arrayLike(t),
      skip?: Uint8Array
    ): VectorArray => {
    const n = t.length
    
    let t_i: number
    let t2: number
    let t3: number

    let coefficients_offset: number
    let a: number
    let b: number
    let c: number
    let d: number
    const dimensionality_times_4 = 4 * dimensionality;
    const axis_times_4 = 4 * axis;
    
    let result_i = offset;
      
    if (skip) {
      for (let i = 0; i < n; i++, result_i += stride) {
        if (skip[i] !== 0) continue;

        t_i = t[i];
        t2 = t_i * t_i;
        t3 = t_i * t2;
        coefficients_offset = (dimensionality_times_4 * coefficient_indices[i]) + axis_times_4;
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset++];
        results[result_i] = a * t3 + b * t2 + c * t_i + d;
      }
    }
    else {
      for (let i = 0; i < n; i++, result_i += stride) {
        t_i = t[i];
        t2 = t_i * t_i;
        t3 = t_i * t2;
        coefficients_offset = (dimensionality_times_4 * coefficient_indices[i]) + axis_times_4;
        a = vectorized_coefficients[coefficients_offset++];
        b = vectorized_coefficients[coefficients_offset++];
        c = vectorized_coefficients[coefficients_offset++];
        d = vectorized_coefficients[coefficients_offset++];
        results[result_i] = a * t3 + b * t2 + c * t_i + d;
      }
    }

    return results;
  }

/**
 * Calculates vector component for the derivative of the curve segment at time t
 * @param t time along curve segment
 * @param coefficients coefficients for curve function
 * @returns derivative (t')
 */
export function derivativeAtT(t: number, coefficients: NumArray4) : number {
  const t2 = t * t;
  const [a, b, c] = coefficients;
  return 3 * a * t2 + 2 * b * t + c;
}

export function derivativeAtT_vectorized<
    TArray extends NumberArrayLike,
    VectorArray extends NumberArrayLike,
  >(
    t: TArray,
    coefficient_indices: IntegerNumberArrayLike,
    dimensionality: number,
    vectorized_coefficients: Float64Array,
    results: VectorArray = <VectorArray><unknown>arrayLike(t)
): VectorArray {
  const n = t.length
  
  let t_i: number
  let t2: number

  let coefficients_offset: number
  let a: number
  let b: number
  let c: number

  for (let i = 0; i < n; i++) {
    t_i = t[i]
    t2 = t_i * t_i
    coefficients_offset = 4 * coefficient_indices[i]
    a = vectorized_coefficients[coefficients_offset + 0]
    b = vectorized_coefficients[coefficients_offset + 1]
    c = vectorized_coefficients[coefficients_offset + 2]
    results[i] = 3 * a * t2 + 2 * b * t_i + c;
  }
  
  return results
}

/**
 * Calculates vector component for the second derivative of the curve segment at time t
 * @param t time along curve segment
 * @param coefficients coefficients for curve function
 * @returns second derivative (t'')
 */
export function secondDerivativeAtT(t: number, coefficients: NumArray4) : number {
  const [a, b] = coefficients;
  return 6 * a * t + 2 * b;
}

export function secondDerivativeAtT_vectorized<
    TArray extends NumberArrayLike,
    VectorArray extends NumberArrayLike,
  >(
    t: TArray,
    coefficient_indices: IntegerNumberArrayLike,
    dimensionality: number,
    vectorized_coefficients: Float64Array,
    results: VectorArray = <VectorArray><unknown>arrayLike(t)
): VectorArray {
  const n = t.length
  
  let t_i: number

  let coefficients_offset: number
  let a: number
  let b: number

  for (let i = 0; i < n; i++) {
    t_i = t[i]
    coefficients_offset = 4 * coefficient_indices[i]
    a = vectorized_coefficients[coefficients_offset + 0]
    b = vectorized_coefficients[coefficients_offset + 1]
    results[i] = 6 * a * t_i + 2 * b;
  }
  
  return results
}

/**
 * Solves the cubic spline equation and return t
 * @param lookup target lookup value
 * @param coefficients lookup axis coefficients
 */
export function findRootsOfT(lookup: number, coefficients: NumArray4): number[] {
  const [a, b, c, d] = coefficients;
  const x = d - lookup;
  if (a === 0 && b === 0 && c === 0 && x === 0) {
    return [0]; // whole segment matches - how to deal with this?
  }
  const roots = getCubicRoots(a, b, c, x);
  return roots.filter(t => t > -EPS && t <= 1 + EPS).map(t => clamp(t, 0, 1));
}

export function findRootsOfT_vectorized<
    VectorArray extends NumberArrayLike,
    RootsArray extends NumberArrayLike
  >(
    lookup: VectorArray,
    root_index: 0 | 1 | 2 | -1,
    coefficients_segment_indices: IntegerNumberArrayLike,
    coefficients_vectorized: Float64Array,
    results: RootsArray = <RootsArray><unknown>arrayLike(lookup),
    skip?: Uint8Array
  ): RootsArray {
  const n = lookup.length;

  let coefficients_offset: number

  let a: number
  let b: number
  let c: number
  let d: number

  let abs_a: number
  let mul_2_a: number
  let sqrt_D: number
  let sqr_B: number
  let mul_27_a_a: number
  let b_div_3_mul_a: number

  let roots_n: 0 | 1 | 2 | 3
  let root_0: number
  let root_1: number
  let root_2: number

  const EPS_negative = -EPS;
  const EPS_positive = 1 + EPS;

  if (skip) {
    for (let i = 0; i < n; i++) {
      if (skip[i] !== 0) continue;

      coefficients_offset = 4 * coefficients_segment_indices[i]
      a = coefficients_vectorized[coefficients_offset++];
      b = coefficients_vectorized[coefficients_offset++];
      c = coefficients_vectorized[coefficients_offset++];
      d = coefficients_vectorized[coefficients_offset++];

      // x = d - lookup[i];
      // d = x
      d -= lookup[i];

      if (a === 0 && b === 0 && c === 0 && d === 0) {
        //   return [0]; // whole segment matches - how to deal with this?
        results[i] = 0;
      }
      else {
        // getCubicRoots(a, b, c, x = d)
        // getCubicRoots(a: number, b: number, c: number, d: number)

        abs_a = Math.abs(a);

        if (abs_a < EPS) { // Quadratic case, ax^2+bx+c=0
          // getQuadRoots(b, c, d);

          if (abs_a < EPS) { // Linear case, ax+b=0
            if (Math.abs(b) < EPS) roots_n = 0; // Degenerate case
            else {
              roots_n = 1;
              root_0 = -c / b;
            }
          }
          else {
            const D = b * b - 4 * a * c;
            if (Math.abs(D) < EPS) {
              roots_n = 1;
              root_0 = -b / (2 * a);
            }
    
            if (D > 0) {
              roots_n = 2;
              sqrt_D = Math.sqrt(D);
              mul_2_a = 2 * a;
              root_0 = (-b + sqrt_D) / mul_2_a;
              root_1 = (-b - sqrt_D) / mul_2_a;
            }
            else roots_n = 0;
          }
        }
        else {
          // Convert to depressed cubic t^3+pt+q = 0 (subst x = t - b/3a)
          // const p = (3 * a * c - b * b) / (3 * a * a);
          // const q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
          const p = (c - b * b) / a;
          sqr_B = b * b;
          mul_27_a_a = 27 * a * a;
          const q = (2 * sqr_B * b - 9 * a * b * c + mul_27_a_a * d) / (mul_27_a_a * a);

          if (Math.abs(p) < EPS) { // p = 0 -> t^3 = -q -> t = -q^1/3
            roots_n = 1;
            root_0 = Math.cbrt(-q);
          } else if (Math.abs(q) < EPS) { // q = 0 -> t^3 + pt = 0 -> t(t^2+p)=0
            root_0 = 0;
            if (p < 0) {
              roots_n = 3;
              root_1 = Math.sqrt(-p);
              // root_2 = -Math.sqrt(-p);
              root_2 = -root_1;
            }
            else {
              roots_n = 1;
            }
          } else {
            const D = q * q / 4 + p * p * p / 27;
            if (Math.abs(D) < EPS) {       // D = 0 -> two roots
              roots_n = 2;
              // root_2 is tmp variable
              root_2 = 1.5 * q / p;
              // root_0 = -1.5 * q / p;
              // root_1 = 3 * q / p;
              root_0 = -root_2;
              root_1 = 2 * root_2;
            } else if (D > 0) {             // Only one real root
              const u = Math.cbrt(-q / 2 - Math.sqrt(D));
              roots_n = 1;
              root_0 = u - p / (3 * u);
            } else { // D < 0, three roots, but needs to use complex numbers/trigonometric solution
              const u = 2 * Math.sqrt(-p / 3);
              const t = Math.acos(3 * q / p / u) / 3;  // D < 0 implies p < 0 and acos argument in [-1..1]
              const k = 2 * Math.PI / 3;
              roots_n = 3;
              root_0 = u * Math.cos(t);
              root_1 = u * Math.cos(t - k);
              root_2 = u * Math.cos(t - 2 * k);
            }
          }
        }

        // Convert back from depressed cubic
        switch (roots_n) {
          case 0:
            break;
          case 1:
            root_0 -= b / (3 * a);
            break;
          case 2:
            b_div_3_mul_a = b / (3 * a);
            root_0 -= b_div_3_mul_a;
            root_1 -= b_div_3_mul_a;
            break;
          case 3:
            b_div_3_mul_a = b / (3 * a);
            root_0 -= b_div_3_mul_a;
            root_1 -= b_div_3_mul_a;
            root_2 -= b_div_3_mul_a;
            break;
        }

        // return roots.filter(t => t > -EPS && t <= 1 + EPS).map(t => clamp(t, 0, 1));
        switch (root_index) {
          case 0:
            switch (roots_n) {
              case 0:
                results[i] = NaN;
                break;
              case 1:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
              case 2:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else
                  results[i] = NaN;
                break;
              case 3:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else if (root_2 > EPS_negative && root_2 <= EPS_positive)
                  results[i] = root_2;
                else
                  results[i] = NaN;
                break;
            }
            break;
      
          case 1:
            switch (roots_n) {
              case 0:
                results[i] = NaN;
                break;
              case 1:
                results[i] = NaN;
                break;
              case 2:
                if (root_0 > EPS_negative && root_0 <= EPS_positive) {
                  if (root_1 > EPS_negative && root_1 <= EPS_positive)
                    results[i] = root_1;
                  else
                    results[i] = NaN;
                }
                else
                  results[i] = NaN;
                break;
              case 3:
                if (root_0 > EPS_negative && root_0 <= EPS_positive) {
                  if (root_1 > EPS_negative && root_1 <= EPS_positive)
                    results[i] = root_1;
                  else if (root_2 > EPS_negative && root_2 <= EPS_positive)
                    results[i] = root_2;
                  else
                    results[i] = NaN;
                }
                else if (root_1 > EPS_negative && root_1 <= EPS_positive) {
                  if (root_2 > EPS_negative && root_2 <= EPS_positive)
                    results[i] = root_2;
                  else
                    results[i] = NaN;
                }
                else
                  results[i] = NaN;
                break;
            }
            break;
      
          case 2:
            if (root_0 > EPS_negative && root_0 <= EPS_positive &&
              root_1 > EPS_negative && root_1 <= EPS_positive &&
              root_2 > EPS_negative && root_2 <= EPS_positive)
              results[i] = root_2;
            else
              results[i] = NaN;
            break;
      
          case -1:
            switch (roots_n) {
              case 0:
                results[i] = NaN;
                break;
              case 1:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
              case 2:
                if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
              case 3:
                if (root_2 > EPS_negative && root_2 <= EPS_positive)
                  results[i] = root_2;
                else if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
            }
            break;
        }
      }
    }
  }
  else {
    for (let i = 0; i < n; i++) {
      coefficients_offset = 4 * coefficients_segment_indices[i]
      a = coefficients_vectorized[coefficients_offset++];
      b = coefficients_vectorized[coefficients_offset++];
      c = coefficients_vectorized[coefficients_offset++];
      d = coefficients_vectorized[coefficients_offset++];
  
      // x = d - lookup[i];
      // d = x
      d -= lookup[i];
  
      if (a === 0 && b === 0 && c === 0 && d === 0) {
        //   return [0]; // whole segment matches - how to deal with this?
        results[i] = 0;
      }
      else {
        // getCubicRoots(a, b, c, x = d)
        // getCubicRoots(a: number, b: number, c: number, d: number)
  
        abs_a = Math.abs(a);
  
        if (abs_a < EPS) { // Quadratic case, ax^2+bx+c=0
          // getQuadRoots(b, c, d);
  
          if (abs_a < EPS) { // Linear case, ax+b=0
            if (Math.abs(b) < EPS) roots_n = 0; // Degenerate case
            else {
              roots_n = 1;
              root_0 = -c / b;
            }
          }
          else {
            const D = b * b - 4 * a * c;
            if (Math.abs(D) < EPS) {
              roots_n = 1;
              root_0 = -b / (2 * a);
            }
      
            if (D > 0) {
              roots_n = 2;
              sqrt_D = Math.sqrt(D);
              mul_2_a = 2 * a;
              root_0 = (-b + sqrt_D) / mul_2_a;
              root_1 = (-b - sqrt_D) / mul_2_a;
            }
            else roots_n = 0;
          }
        }
        else {
          // Convert to depressed cubic t^3+pt+q = 0 (subst x = t - b/3a)
          // const p = (3 * a * c - b * b) / (3 * a * a);
          // const q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
          const p = (c - b * b) / a;
          sqr_B = b * b;
          mul_27_a_a = 27 * a * a;
          const q = (2 * sqr_B * b - 9 * a * b * c + mul_27_a_a * d) / (mul_27_a_a * a);
  
          if (Math.abs(p) < EPS) { // p = 0 -> t^3 = -q -> t = -q^1/3
            roots_n = 1;
            root_0 = Math.cbrt(-q);
          } else if (Math.abs(q) < EPS) { // q = 0 -> t^3 + pt = 0 -> t(t^2+p)=0
            root_0 = 0;
            if (p < 0) {
              roots_n = 3;
              root_1 = Math.sqrt(-p);
              // root_2 = -Math.sqrt(-p);
              root_2 = -root_1;
            }
            else {
              roots_n = 1;
            }
          } else {
            const D = q * q / 4 + p * p * p / 27;
            if (Math.abs(D) < EPS) {       // D = 0 -> two roots
              roots_n = 2;
              // root_2 is tmp variable
              root_2 = 1.5 * q / p;
              // root_0 = -1.5 * q / p;
              // root_1 = 3 * q / p;
              root_0 = -root_2;
              root_1 = 2 * root_2;
            } else if (D > 0) {             // Only one real root
              const u = Math.cbrt(-q / 2 - Math.sqrt(D));
              roots_n = 1;
              root_0 = u - p / (3 * u);
            } else { // D < 0, three roots, but needs to use complex numbers/trigonometric solution
              const u = 2 * Math.sqrt(-p / 3);
              const t = Math.acos(3 * q / p / u) / 3;  // D < 0 implies p < 0 and acos argument in [-1..1]
              const k = 2 * Math.PI / 3;
              roots_n = 3;
              root_0 = u * Math.cos(t);
              root_1 = u * Math.cos(t - k);
              root_2 = u * Math.cos(t - 2 * k);
            }
          }
        }
  
        // Convert back from depressed cubic
        switch (roots_n) {
          case 0:
            break;
          case 1:
            root_0 -= b / (3 * a);
            break;
          case 2:
            b_div_3_mul_a = b / (3 * a);
            root_0 -= b_div_3_mul_a;
            root_1 -= b_div_3_mul_a;
            break;
          case 3:
            b_div_3_mul_a = b / (3 * a);
            root_0 -= b_div_3_mul_a;
            root_1 -= b_div_3_mul_a;
            root_2 -= b_div_3_mul_a;
            break;
        }
  
        // return roots.filter(t => t > -EPS && t <= 1 + EPS).map(t => clamp(t, 0, 1));
        switch (root_index) {
          case 0:
            switch (roots_n) {
              case 0:
                results[i] = NaN;
                break;
              case 1:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
              case 2:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else
                  results[i] = NaN;
                break;
              case 3:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else if (root_2 > EPS_negative && root_2 <= EPS_positive)
                  results[i] = root_2;
                else
                  results[i] = NaN;
                break;
            }
            break;
        
          case 1:
            switch (roots_n) {
              case 0:
                results[i] = NaN;
                break;
              case 1:
                results[i] = NaN;
                break;
              case 2:
                if (root_0 > EPS_negative && root_0 <= EPS_positive) {
                  if (root_1 > EPS_negative && root_1 <= EPS_positive)
                    results[i] = root_1;
                  else
                    results[i] = NaN;
                }
                else
                  results[i] = NaN;
                break;
              case 3:
                if (root_0 > EPS_negative && root_0 <= EPS_positive) {
                  if (root_1 > EPS_negative && root_1 <= EPS_positive)
                    results[i] = root_1;
                  else if (root_2 > EPS_negative && root_2 <= EPS_positive)
                    results[i] = root_2;
                  else
                    results[i] = NaN;
                }
                else if (root_1 > EPS_negative && root_1 <= EPS_positive) {
                  if (root_2 > EPS_negative && root_2 <= EPS_positive)
                    results[i] = root_2;
                  else
                    results[i] = NaN;
                }
                else
                  results[i] = NaN;
                break;
            }
            break;
        
          case 2:
            if (root_0 > EPS_negative && root_0 <= EPS_positive &&
              root_1 > EPS_negative && root_1 <= EPS_positive &&
              root_2 > EPS_negative && root_2 <= EPS_positive)
              results[i] = root_2;
            else
              results[i] = NaN;
            break;
        
          case -1:
            switch (roots_n) {
              case 0:
                results[i] = NaN;
                break;
              case 1:
                if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
              case 2:
                if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
              case 3:
                if (root_2 > EPS_negative && root_2 <= EPS_positive)
                  results[i] = root_2;
                else if (root_1 > EPS_negative && root_1 <= EPS_positive)
                  results[i] = root_1;
                else if (root_0 > EPS_negative && root_0 <= EPS_positive)
                  results[i] = root_0;
                else
                  results[i] = NaN;
                break;
            }
            break;
        }
      }
    }
  }
  
  return results;
}

/**
 * Convenience function for evaluating segment functions for all components of a vector
 * @param func SegmentFunction to evaluate
 * @param t time along curve segment
 * @param coefficients coefficients for curve function (for each component)
 * @param target target vector
 * @returns vector
 */
export function evaluateForT(func: SegmentFunction, t:number, coefficients: NumArray4[], target: Vector = null) : Vector {
  target = target || new Array(coefficients.length);

  for (let k = 0; k < coefficients.length; k++) {
    target[k] = func(t, coefficients[k]);
  }

  return target;
}

/**
 * Symmetric vectorized function for convenience function for evaluating
 * segment functions for all components of a vector. It calls the given
 * SegmentFunction_vectorized
 * @param func SegmentFunction_vectorized to evaluate
 * @param t times along curve segment
 * @param coefficients_segment_indices 
 * @param dimensionality dimensionality of coefficients and result vectors
 * @param coefficients_vectorized coefficients for curve function (for each component)
 * @param target target vector
 * @returns vector
 */
export function evaluateForT_vectorized<
    TArray extends NumberArrayLike,
    VectorArray extends NumberArrayLike,
  >(
    func: SegmentFunction_vectorized<TArray, VectorArray>,
    t: TArray,
    coefficients_segment_indices: IntegerNumberArrayLike,
    dimensionality: number,
    coefficients_vectorized: Float64Array,
    target: VectorArray = <VectorArray><unknown>arrayLike(t, dimensionality),
    skip?: Uint8Array
  ): VectorArray {
  return func(t, coefficients_segment_indices, dimensionality, coefficients_vectorized, target, skip)
}

