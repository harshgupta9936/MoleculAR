import { useState } from 'react';
import { StackedBar } from './StackedBar';
import { HubnessBadge } from './HubnessBadge';
import { CoverageBadge } from './CoverageBadge';
import { ExpandedDetail } from './ExpandedDetail';

export function ResultRow({
  candidate,
  rank,
  isCompared = false,
  onToggleCompare,
  canCompare = true,
}) {
  const [isExpanded, setIsExpanded] = useState(false);

  return (
    <div className={`result-row-card ${isExpanded ? 'is-expanded' : ''} ${isCompared ? 'is-compared' : ''}`}>
      <div
        className="result-row"
        onClick={() => setIsExpanded((prev) => !prev)}
        role="button"
        tabIndex={0}
        onKeyDown={(e) => {
          if (e.key === 'Enter' || e.key === ' ') {
            e.preventDefault();
            setIsExpanded((prev) => !prev);
          }
        }}
        aria-expanded={isExpanded}
        aria-label={`Rank ${rank}: ${candidate.candidate}, score ${candidate.score.toFixed(4)}`}
      >
        {/* Compare Checkbox / Button */}
        <div className="col-compare" onClick={(e) => e.stopPropagation()}>
          <button
            type="button"
            className={`btn-compare-pin ${isCompared ? 'is-pinned' : ''}`}
            disabled={!isCompared && !canCompare}
            onClick={() => onToggleCompare && onToggleCompare(candidate)}
            title={
              isCompared
                ? 'Remove from compare tray'
                : canCompare
                ? 'Pin to compare tray (up to 3)'
                : 'Compare tray is full (3/3)'
            }
            aria-label={`Pin ${candidate.candidate} to compare tray`}
            aria-pressed={isCompared}
          >
            {isCompared ? '★' : '☆'}
          </button>
        </div>

        <div className="col-rank">
          <span className="rank-badge">#{rank}</span>
        </div>

        <div className="col-name">
          <div className="candidate-name-row">
            <span className="candidate-name">{candidate.candidate}</span>
            {candidate.is_duplicate_ingredient && (
              <span
                className="badge-duplicate"
                title="Same molecular skeleton as the query (salt, ester or stereoisomer, for example an enantiomer). The evaluation counts such pairs as duplicates and reports results with them neutralised. Stereoisomers can be distinct drugs. This badge covers duplicates within the evaluation pool only."
              >
                Salt or stereoisomer form
              </span>
            )}
            <HubnessBadge percentile={candidate.hubness_percentile} />
          </div>
          <CoverageBadge
            hasTargets={candidate.has_targets}
            hasSideEffects={candidate.has_side_effects}
          />
        </div>

        <div className="col-score">
          <span className="score-mono">{candidate.score.toFixed(4)}</span>
        </div>

        <div className="col-bar" onClick={(e) => e.stopPropagation()}>
          <StackedBar candidate={candidate} />
        </div>

        <div className="col-expand">
          <button
            type="button"
            className="expand-toggle-btn"
            aria-label={isExpanded ? 'Collapse details' : 'Expand details'}
            tabIndex={-1}
          >
            <svg
              className={`chevron-icon ${isExpanded ? 'open' : ''}`}
              viewBox="0 0 16 16"
              width="14"
              height="14"
              fill="currentColor"
              aria-hidden="true"
            >
              <path fillRule="evenodd" d="M1.646 4.646a.5.5 0 0 1 .708 0L8 10.293l5.646-5.647a.5.5 0 0 1 .708.708l-6 6a.5.5 0 0 1-.708 0l-6-6a.5.5 0 0 1 0-.708z" />
            </svg>
          </button>
        </div>
      </div>

      {isExpanded && (
        <ExpandedDetail candidate={candidate} />
      )}
    </div>
  );
}
