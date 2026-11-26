/**
 * Score decomposition helpers.
 * Score calculation from Phase 4.4 hand-weighted pipeline:
 *   pre_penalty = 0.4 * struct_sim + 0.4 * target_sim + 0.2 * side_effect_sim
 *   penalty_factor = 1.0 + (hubness_percentile / 100.0) ^ 2
 *   final_score = pre_penalty / penalty_factor
 *
 * Values are rounded to 4 decimal places for display only; full float used for bar widths.
 */

export function calculateDecomposition(candidate) {
  const s = candidate.struct_sim ?? 0;
  const t = candidate.target_sim ?? 0;
  const e = candidate.side_effect_sim ?? 0;

  const wStruct = 0.4 * s;
  const wTarget = 0.4 * t;
  const wSideEffect = 0.2 * e;

  const prePenalty = wStruct + wTarget + wSideEffect;

  const hubPct = candidate.hubness_percentile;
  const hasHubPenalty = hubPct !== null && hubPct !== undefined;
  const penaltyFactor = hasHubPenalty ? 1.0 + Math.pow(hubPct / 100.0, 2) : 1.0;

  const finalScore = prePenalty / penaltyFactor;

  return {
    rawStruct: s,
    rawTarget: t,
    rawSideEffect: e,
    wStruct,
    wTarget,
    wSideEffect,
    prePenalty,
    penaltyFactor,
    finalScore,
    hubPct,
    hasHubPenalty,
  };
}
