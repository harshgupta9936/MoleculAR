import { calculateDecomposition } from '../utils/scoreDecomposition';

/**
 * StackedBar component.
 * Displays weighted similarity components:
 *   - Structure (40% weight, #0077BB)
 *   - Target (40% weight, #EE7733)
 *   - Side Effects (20% weight, #009988)
 *
 * Values are rounded to 4 decimal places for display only; full float used for bar widths.
 * A vertical marker line indicates the final penalised score position.
 */
export function StackedBar({ candidate }) {
  const decomp = calculateDecomposition(candidate);

  // Maximum theoretical pre-penalty score is 1.0 (0.4 + 0.4 + 0.2)
  const maxScore = 1.0;
  const structWidthPct = Math.min(100, Math.max(0, (decomp.wStruct / maxScore) * 100));
  const targetWidthPct = Math.min(100, Math.max(0, (decomp.wTarget / maxScore) * 100));
  const seWidthPct = Math.min(100, Math.max(0, (decomp.wSideEffect / maxScore) * 100));

  // The penalised final score position on the 0-1 scale
  const finalScorePct = Math.min(100, Math.max(0, (candidate.score / maxScore) * 100));

  const ariaDescription = `Decomposition for ${candidate.candidate}: Structure weighted ${decomp.wStruct.toFixed(4)}, Target weighted ${decomp.wTarget.toFixed(4)}, Side Effects weighted ${decomp.wSideEffect.toFixed(4)}, Pre-penalty sum ${decomp.prePenalty.toFixed(4)}, Final score ${candidate.score.toFixed(4)}`;

  return (
    <div
      className="stacked-bar-container"
      tabIndex={0}
      role="img"
      aria-label={ariaDescription}
    >
      <div className="stacked-bar-track">
        {/* Structure segment (0.4 * struct_sim) */}
        {structWidthPct > 0 && (
          <div
            className="bar-segment bar-segment-struct"
            style={{ width: `${structWidthPct}%` }}
            title={`Structure: raw ${decomp.rawStruct.toFixed(4)} × 0.4 = ${decomp.wStruct.toFixed(4)}`}
          />
        )}

        {/* Target segment (0.4 * target_sim) */}
        {targetWidthPct > 0 && (
          <div
            className="bar-segment bar-segment-target"
            style={{ width: `${targetWidthPct}%` }}
            title={`Target: raw ${decomp.rawTarget.toFixed(4)} × 0.4 = ${decomp.wTarget.toFixed(4)}`}
          />
        )}

        {/* Side Effect segment (0.2 * side_effect_sim) */}
        {seWidthPct > 0 && (
          <div
            className="bar-segment bar-segment-se"
            style={{ width: `${seWidthPct}%` }}
            title={`Side Effects: raw ${decomp.rawSideEffect.toFixed(4)} × 0.2 = ${decomp.wSideEffect.toFixed(4)}`}
          />
        )}

        {/* Penalised Final Score Marker */}
        <div
          className="bar-marker"
          style={{ left: `${finalScorePct}%` }}
          title={`Penalised Score: ${candidate.score.toFixed(4)} (pre-penalty: ${decomp.prePenalty.toFixed(4)}, penalty factor: ${decomp.penaltyFactor.toFixed(4)})`}
        />
      </div>

      {/* Tooltip on hover / focus */}
      <div className="stacked-bar-tooltip" role="tooltip">
        <div className="tooltip-row">
          <span className="tooltip-dot struct-dot" />
          <span>Structure:</span>
          <span className="mono-num">{decomp.rawStruct.toFixed(4)} × 0.4 = <strong>{decomp.wStruct.toFixed(4)}</strong></span>
        </div>
        <div className="tooltip-row">
          <span className="tooltip-dot target-dot" />
          <span>Target:</span>
          <span className="mono-num">{decomp.rawTarget.toFixed(4)} × 0.4 = <strong>{decomp.wTarget.toFixed(4)}</strong></span>
        </div>
        <div className="tooltip-row">
          <span className="tooltip-dot se-dot" />
          <span>Side Effects:</span>
          <span className="mono-num">{decomp.rawSideEffect.toFixed(4)} × 0.2 = <strong>{decomp.wSideEffect.toFixed(4)}</strong></span>
        </div>
        <div className="tooltip-divider" />
        <div className="tooltip-row summary-row">
          <span>Pre-penalty sum:</span>
          <span className="mono-num"><strong>{decomp.prePenalty.toFixed(4)}</strong></span>
        </div>
        {decomp.hasHubPenalty && (
          <div className="tooltip-row hub-row">
            <span>Hub penalty (P{decomp.hubPct?.toFixed(0)}):</span>
            <span className="mono-num">÷ {decomp.penaltyFactor.toFixed(4)}</span>
          </div>
        )}
        <div className="tooltip-row final-row">
          <span>Final score:</span>
          <span className="mono-num final-val"><strong>{candidate.score.toFixed(4)}</strong></span>
        </div>
      </div>
    </div>
  );
}
