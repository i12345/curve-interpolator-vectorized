export type WritableArrayLike<T> = {
    length: number
    [index: number]: T
}

export type NumberArrayLike = WritableArrayLike<number>
export type IntegerNumberArrayLike = number[] | Uint8Array | Int8Array | Uint16Array | Int16Array | Uint32Array | Int32Array

export function arrayLike<Array extends NumberArrayLike>(array: Array, factor: number | { factor: number } | { divisor: number } = 1): Array {
    const length = typeof factor === 'number' ?
        factor * array.length :
        'factor' in factor ?
            array.length * factor.factor :
            array.length / factor.divisor
    
    return new (<{ new(length: number): Array }>array.constructor)(length)
}