import { describe, it, expect } from 'vitest';
import { calculateDecomposition } from '../utils/scoreDecomposition';

describe('scoreDecomposition', () => {
  it('correctly calculates weighted components summing to pre_penalty', () => {
    const candidate = {
      struct_sim: 0.5,
      target_sim: 0.8,
      side_effect_sim: 0.3,
      hubness_percentile: 50.0,
      score: 0.464, // (0.2 + 0.32 + 0.06) / (1 + 0.25) = 0.58 / 1.25 = 0.464
    };

    const decomp = calculateDecomposition(candidate);

    expect(decomp.wStruct).toBeCloseTo(0.2, 8);
    expect(decomp.wTarget).toBeCloseTo(0.32, 8);
    expect(decomp.wSideEffect).toBeCloseTo(0.06, 8);

    const sumWeighted = decomp.wStruct + decomp.wTarget + decomp.wSideEffect;
    expect(decomp.prePenalty).toBeCloseTo(sumWeighted, 8);
    expect(decomp.prePenalty).toBeCloseTo(0.58, 8);
  });

  it('correctly applies hub penalty factor: pre_penalty / (1 + (hub/100)^2)', () => {
    const candidate = {
      struct_sim: 0.65,
      target_sim: 0.45,
      side_effect_sim: 0.25,
      hubness_percentile: 57.4768,
      score: 0.3705,
    };

    const decomp = calculateDecomposition(candidate);

    const expectedPrePenalty = 0.4 * 0.65 + 0.4 * 0.45 + 0.2 * 0.25; // 0.26 + 0.18 + 0.05 = 0.49
    expect(decomp.prePenalty).toBeCloseTo(expectedPrePenalty, 8);

    const pctRank = 57.4768 / 100.0;
    const expectedPenaltyFactor = 1.0 + pctRank * pctRank;
    expect(decomp.penaltyFactor).toBeCloseTo(expectedPenaltyFactor, 8);

    const expectedFinal = expectedPrePenalty / expectedPenaltyFactor;
    expect(decomp.finalScore).toBeCloseTo(expectedFinal, 8);
  });

  it('handles zero or missing hubness percentile with penaltyFactor=1.0', () => {
    const candidate = {
      struct_sim: 0.4,
      target_sim: 0.2,
      side_effect_sim: 0.1,
      hubness_percentile: null,
      score: 0.26,
    };

    const decomp = calculateDecomposition(candidate);
    expect(decomp.penaltyFactor).toBe(1.0);
    expect(decomp.finalScore).toBeCloseTo(0.26, 8);
  });
});
