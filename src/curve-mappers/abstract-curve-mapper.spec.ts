import 'mocha';
import * as sinon from 'sinon';
import { expect } from 'chai';
import { AbstractCurveMapper } from './abstract-curve-mapper';
import { points, points3d } from '../../test/test-data';
import { NumberArrayLike } from '../core/array';


class TestMapper extends AbstractCurveMapper<Float64Array> {
  lengthAt(u: number): number {
    throw new Error('Method not implemented.');
  }
  lengthAt_vectorized<UArray extends NumberArrayLike, LengthArray extends NumberArrayLike>(u: UArray, length?: LengthArray): LengthArray {
    throw new Error('Method not implemented.');
  }
  getT(u: number): number {
    throw new Error('Method not implemented.');
  }
  getT_vectorized<UArray extends NumberArrayLike, TArray extends NumberArrayLike>(u: UArray, t?: TArray): TArray {
    throw new Error('Method not implemented.');
  }
  getU(t: number): number {
    throw new Error('Method not implemented.');
  }
  getU_vectorized<TArray extends NumberArrayLike, UArray extends NumberArrayLike>(t: TArray, u?: UArray): UArray {
    throw new Error('Method not implemented.');
  }
}

describe('abstract-curve-mapper.ts', () => {
  it('should be able to instantiate class', () => {
    const mapper = new TestMapper();
    expect(mapper).to.not.be.null;
  });

  it('should be able set parameters', () => {
    const mapper = new TestMapper();

    // test defaults
    expect(mapper.alpha).to.equal(0);
    expect(mapper.tension).to.equal(0.5);
    expect(mapper.closed).to.be.false;
    expect(mapper.points).to.be.undefined;

    // should not be allowed to pass less than 3 control points
    expect(() => mapper.points = []).to.throw;

    mapper.alpha = 0.5;
    mapper.tension = 0;
    mapper.closed = true;
    mapper.points = [[], [], []];

    expect(mapper.alpha).to.equal(0.5);
    expect(mapper.tension).to.equal(0);
    expect(mapper.closed).to.be.true;
    expect(mapper.points).to.deep.equal([[], [], []]);
  });

  it('should be able to calculate and cache coefficients and invalidate cache if parameters are changed', () => {
    const mapper = new TestMapper();
    mapper.points = points;
    const coefficients = mapper.getCoefficients(1);
    expect(coefficients).to.be.instanceOf(Array);
    expect(mapper._cache['coefficients'].has(1));

    mapper.alpha = 0; // same as default
    expect(mapper._cache['coefficients']).to.not.be.null;

    mapper.alpha = 0.5;
    expect(mapper._cache['coefficients']).to.be.null;
  });

  it('should invoke callback when cache is invalidated', () => {
    const callback = sinon.spy();
    const mapper = new TestMapper(callback);
    mapper.points = [[],[]];
    mapper.alpha = 0.5;
    expect(callback.called).to.be.true;
  });


});
