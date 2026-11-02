'use strict';
const assert = require('node:assert/strict');
const metrics = require('../assets/js/metric-demo.js');
const approx = (got, want, tol = 1e-10) => {
  assert.ok(Math.abs(got - want) <= tol, `expected ${want}, got ${got}`);
};

// Midranks and one-vs-rest U. Includes an exact expression tie.
assert.deepEqual(metrics.averageRanks([3, 2, 2, 5]), [3, 1.5, 1.5, 4]);
const gt = [[1, 4, 2, 3], [2, 5, 3, 1], [1, 1, 7, 8], [2, 2, 6, 7]];
const pred = [[2, 4, 2, 3], [3, 4, 3, 1], [1, 1, 6, 7], [2, 1, 7, 8]];
const group = ['a', 'a', 'b', 'b'];
assert.deepEqual(metrics.mannWhitneyU(gt, group)['a'], [2, 4, 0, 0]);
assert.deepEqual(metrics.mannWhitneyU(gt, group)['b'], [2, 0, 4, 4]);
const result = metrics.evaluate(gt, pred, group, 2);
assert.equal(result.clusters, 2);
assert.equal(result.k, 2);
assert.ok(result.ndcg >= 0 && result.ndcg <= 1);
assert.ok(result.spearman >= -1 && result.spearman <= 1);

// sklearn-style tie averaging: tied prediction means half relevant, half not.
approx(metrics.binaryNdcg([3, 3, 1], [3, 1, 2], 1), 0.5);
approx(metrics.binaryNdcg([4, 3, 2], [4, 3, 2], 2), 1);
approx(metrics.spearman([1, 2, 2, 3], [1, 2, 2, 3]), 1);
const identical = metrics.evaluate(gt, gt, group, 2);
approx(identical.spearman, 1);
approx(identical.ndcg, 1);
approx(identical.overlap, 1);

const toy = metrics.makeToyData(0, 123);
const perfect = metrics.evaluate(toy.gt, toy.pred, toy.clusters, 10);
approx(perfect.spearman, 1);
approx(perfect.ndcg, 1);
approx(perfect.overlap, 1);
console.log('All IDER metric demo tests passed.');
