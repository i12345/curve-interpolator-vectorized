import { IntegerNumberArrayLike, NumberArrayLike } from "./array";
import { Vector } from "./interfaces";

/**
 * Fill all components of a vector with a value
 * @param v vector
 * @param val fill value
 */
export function fill(v:Vector, val:number) : Vector {
  for (let i = 0; i < v.length; i++) {
    v[i] = val;
  }
  return v;
}

/**
 * Map all components of a vector using provided mapping function.
 * @param v vector
 * @param func mapping function
 */
export function map(v:Vector, func: (c:number, i:number) => number) : Vector {
  for (let i = 0; i < v.length; i++) {
    v[i] = func(v[i], i);
  }
  return v;
}

/**
 * Reduce a vector to a single value using the provided reduce function.
 * @param v vector
 * @param func reduce function
 * @param r initial value
 */
export function reduce(v:Vector, func: (s: number, c:number, i:number) => number, r = 0) : number {
  for (let i = 0; i < v.length; i++) {
    r = func(r, v[i], i);
  }
  return r;
}

/**
 * Copy values from one vector to another. If target is not provided it will be created.
 * @param source source vector
 * @param target target vector
 * @returns vector
 */
export function copyValues(source: Vector, target?: Vector) : Vector {
  target = target || new Array(source.length);
  for (let i = 0; i < source.length; i++) {
    target[i] = source[i];
  }
  return target;
}

/**
 * Reduce the set of coordinates for a curve by eliminating points that are not
 * contributing to the shape of the curve, i.e. multiple points making out a linear
 * segment.
 * @param inputArr set of coordinates
 * @param maxOffset threshold to use for determining if a point is part of a linear line segment
 * @param maxDistance points will not be removed if the distance equals or is greater than the given maxDistance
 */
export function simplify2d(inputArr:number[][], maxOffset = 0.001, maxDistance = 10) : number[][] {
  if (inputArr.length <= 4) return inputArr;
  const [o0, o1] = inputArr[0];
  const arr = inputArr.map(d => [d[0] - o0, d[1] - o1]);
  let [a0, a1] = arr[0];
  const sim = [inputArr[0]];

  for (let i = 1; i + 1 < arr.length; i++) {
    const [t0, t1] = arr[i];
    const [b0, b1] = arr[i + 1];

    if (b0 - t0 !== 0 || b1 - t1 !== 0) {
      // Proximity check
      const proximity =
        Math.abs(a0 * b1 - a1 * b0 + b0 * t1 - b1 * t0 + a1 * t0 - a0 * t1) /
        Math.sqrt((b0 - a0) ** 2 + (b1 - a1) ** 2);

      const dir = [a0 - t0, a1 - t1];
      const len = Math.sqrt(dir[0] ** 2 + dir[1] ** 2);

      if (proximity > maxOffset || len >= maxDistance) {
        sim.push([t0 + o0, t1 + o1]);
        [a0, a1] = [t0, t1];
      }
    }
  }
  const last = arr[arr.length - 1];
  sim.push([last[0] + o0, last[1] + o1]);

  return sim;
}

/**
 * Clamp an input value to min and max
 * @param value input value
 * @param min min value
 * @param max max value
 */
export function clamp(value:number, min = 0, max = 1) : number {
  if (value < min) return min;
  if (value > max) return max;
  return value;
}

export function clamp_vectorized<Array extends NumberArrayLike>(values: Array, min = 0, max = 1): Array {
  let value: number
  for (let i = 0; i < values.length; i++) {
    value = values[i]
    if (value < min) values[i] = min
    else if(value > max) values[i] = max
  }
  return values
}

/**
 * Finds the index in accumulatedValues of the highest value that is less than or equal to targetValue
 * @param targetValue search term
 * @param accumulatedValues array of accumulated values to search in
 * @returns
 */
export function binarySearch(targetValue: number, accumulatedValues: NumberArrayLike) {
  const min = accumulatedValues[0];
  const max = accumulatedValues[accumulatedValues.length - 1];
  if (targetValue >= max) {
    return accumulatedValues.length - 1;
  }

  if (targetValue <= min) {
    return 0;
  }

  let left = 0;
  let right = accumulatedValues.length - 1;

  while (left <= right) {
    const mid = Math.floor((left + right) / 2);
    const lMid = accumulatedValues[mid];

    if (lMid < targetValue) {
      left = mid + 1;
    } else if (lMid > targetValue) {
      right = mid - 1;
    } else {
      return mid;
    }
  }

  return Math.max(0, right);
}

export function binarySearch_vectorized<
    TargetValuesArray extends NumberArrayLike,
    AccumulatedValuesArray extends NumberArrayLike = TargetValuesArray,
    IndicesArray extends IntegerNumberArrayLike = Uint32Array
  >(
    targetValues: TargetValuesArray,
    accumulatedValues: AccumulatedValuesArray,
    results: IndicesArray = <IndicesArray>new Uint32Array(targetValues.length),
    skip?: Uint8Array
  ): IndicesArray {
  const length_minus_1 = accumulatedValues.length - 1
  const min = accumulatedValues[0];
  const max = accumulatedValues[length_minus_1];
  
  let left: number
  let right: number
  let mid: number
  let lMid: number
  let foundTarget: boolean
  let targetValue: number

  if (skip) {
    for (let i = 0; i < targetValues.length; i++) {
      if(skip[i] !== 0) continue

      targetValue = targetValues[i]
    
      if (targetValue >= max) {
        results[i] = length_minus_1
        continue
      }
      else if (targetValue <= min) {
        results[i] = 0
        continue
      }

      left = 0
      right = length_minus_1
      foundTarget = false

      while (left <= right) {
        mid = Math.floor((left + right) / 2);
        lMid = accumulatedValues[mid]

        if (lMid < targetValue)
          left = mid + 1
        else if (lMid > targetValue)
          right = mid - 1
        else {
          results[i] = mid
          foundTarget = true
          break
        }
      }
    
      if (!foundTarget)
        results[i] = Math.max(0, right);
    }
  }
  else {
    for (let i = 0; i < targetValues.length; i++) {
      targetValue = targetValues[i]
    
      if (targetValue >= max) {
        results[i] = length_minus_1
        continue
      }
      else if (targetValue <= min) {
        results[i] = 0
        continue
      }

      left = 0
      right = length_minus_1
      foundTarget = false

      while (left <= right) {
        mid = Math.floor((left + right) / 2);
        lMid = accumulatedValues[mid]

        if (lMid < targetValue)
          left = mid + 1
        else if (lMid > targetValue)
          right = mid - 1
        else {
          results[i] = mid
          foundTarget = true
          break
        }
      }
    
      if (!foundTarget)
        results[i] = Math.max(0, right);
    }
  }

  return results
}

export function binarySearch_vectorized_specific<
    TargetValuesArray extends NumberArrayLike,
    AccumulatedValuesIndexArray extends IntegerNumberArrayLike,
    AccumulatedValuesArray extends NumberArrayLike = TargetValuesArray,
    IndicesArray extends IntegerNumberArrayLike = Uint32Array
  >(
    targetValues: TargetValuesArray,
    accumulatedValues_idx: AccumulatedValuesIndexArray,
    accumulatedValues_arrays: AccumulatedValuesArray[],
    results: IndicesArray = <IndicesArray>new Uint32Array(targetValues.length),
    skip?: Uint8Array
  ): IndicesArray {
  let left: number
  let right: number
  let mid: number
  let lMid: number
  let foundTarget: boolean
  let targetValue: number
  let accumulatedValues: AccumulatedValuesArray

  let length_minus_1: number;
  let min: number;
  let max: number;

  if (skip) {
    for (let i = 0; i < targetValues.length; i++) {
      if (skip[i] !== 0) continue;

      accumulatedValues = accumulatedValues_arrays[accumulatedValues_idx[i]]
      targetValue = targetValues[i]
    
      length_minus_1 = accumulatedValues.length - 1
      min = accumulatedValues[0];
      max = accumulatedValues[length_minus_1];
  
      if (targetValue >= max) {
        results[i] = length_minus_1
        continue
      }
      else if (targetValue <= min) {
        results[i] = 0
        continue
      }

      left = 0
      right = length_minus_1
      foundTarget = false

      while (left <= right) {
        mid = Math.floor((left + right) / 2);
        lMid = accumulatedValues[mid]

        if (lMid < targetValue)
          left = mid + 1
        else if (lMid > targetValue)
          right = mid - 1
        else {
          results[i] = mid
          foundTarget = true
          break
        }
      }
    
      if (!foundTarget)
        results[i] = Math.max(0, right);
    }
  }
  else {
    for (let i = 0; i < targetValues.length; i++) {
      accumulatedValues = accumulatedValues_arrays[accumulatedValues_idx[i]]
      targetValue = targetValues[i]
    
      const length_minus_1 = accumulatedValues.length - 1
      const min = accumulatedValues[0];
      const max = accumulatedValues[length_minus_1];
  
      if (targetValue >= max) {
        results[i] = length_minus_1
        continue
      }
      else if (targetValue <= min) {
        results[i] = 0
        continue
      }

      left = 0
      right = length_minus_1
      foundTarget = false

      while (left <= right) {
        mid = Math.floor((left + right) / 2);
        lMid = accumulatedValues[mid]

        if (lMid < targetValue)
          left = mid + 1
        else if (lMid > targetValue)
          right = mid - 1
        else {
          results[i] = mid
          foundTarget = true
          break
        }
      }
    
      if (!foundTarget)
        results[i] = Math.max(0, right);
    }
  }

  return results
}
