import { useState } from 'react';
import { API_BASE_URL } from '../utils/api';

export function CompareTray({
  pinnedCandidates,
  onRemove,
  onClear,
}) {
  const [isCollapsed, setIsCollapsed] = useState(false);

  if (!pinnedCandidates || pinnedCandidates.length === 0) {
    return null;
  }

  // Calculate shared intersection of targets across pinned candidates
  const commonTargets = pinnedCandidates.reduce((acc, c, idx) => {
    const set = new Set(c.shared_targets || []);
    if (idx === 0) return set;
    return new Set([...acc].filter((x) => set.has(x)));
  }, new Set());

  // Calculate shared intersection of side effects across pinned candidates
  const commonSideEffects = pinnedCandidates.reduce((acc, c, idx) => {
    const set = new Set(c.shared_side_effects || []);
    if (idx === 0) return set;
    return new Set([...acc].filter((x) => set.has(x)));
  }, new Set());

  return (
    <div className={`compare-tray ${isCollapsed ? 'is-collapsed' : ''}`} role="region" aria-label="Candidate Comparison Tray">
      <div className="tray-header">
        <div className="tray-title-group">
          <span className="tray-badge">{pinnedCandidates.length}/3</span>
          <span className="tray-title">Compare Candidates</span>
        </div>
        <div className="tray-header-actions">
          <button
            type="button"
            className="tray-btn-text"
            onClick={onClear}
            title="Clear all compared candidates"
          >
            Clear All
          </button>
          <button
            type="button"
            className="tray-btn-icon"
            onClick={() => setIsCollapsed((prev) => !prev)}
            aria-label={isCollapsed ? 'Expand compare tray' : 'Collapse compare tray'}
          >
            {isCollapsed ? '▲' : '▼'}
          </button>
        </div>
      </div>

      {!isCollapsed && (
        <div className="tray-body">
          <div className="compare-grid">
            {pinnedCandidates.map((cand) => (
              <div key={cand.candidate} className="compare-card">
                <div className="compare-card-top">
                  <div>
                    <div className="compare-name">{cand.candidate}</div>
                    {cand.is_duplicate_ingredient && (
                      <span
                        className="badge-duplicate"
                        title="Same molecular skeleton as the query (salt, ester or stereoisomer, for example an enantiomer). The evaluation counts such pairs as duplicates and reports results with them neutralised. Stereoisomers can be distinct drugs. This badge covers duplicates within the evaluation pool only."
                        style={{ display: 'inline-block', marginTop: '0.2rem' }}
                      >
                        Salt or stereoisomer form
                      </span>
                    )}
                  </div>
                  <button
                    type="button"
                    className="compare-remove-btn"
                    onClick={() => onRemove(cand.candidate)}
                    aria-label={`Remove ${cand.candidate} from compare`}
                  >
                    ×
                  </button>
                </div>

                {/* 2D Structure */}
                <div className="compare-structure">
                  <img
                    src={`${API_BASE_URL}/api/structure?name=${encodeURIComponent(cand.candidate)}`}
                    alt={`2D Structure of ${cand.candidate}`}
                    className="compare-svg"
                    loading="lazy"
                  />
                </div>

                {/* Metrics */}
                <div className="compare-metrics">
                  <div className="compare-metric-row">
                    <span className="metric-label">Score</span>
                    <span className="metric-val mono-num score-highlight">{cand.score.toFixed(4)}</span>
                  </div>
                  <div className="compare-metric-row">
                    <span className="metric-label">Structure</span>
                    <span className="metric-val mono-num">{cand.struct_sim.toFixed(4)}</span>
                  </div>
                  <div className="compare-metric-row">
                    <span className="metric-label">Target</span>
                    <span className="metric-val mono-num">{cand.target_sim.toFixed(4)}</span>
                  </div>
                  <div className="compare-metric-row">
                    <span className="metric-label">Side Effect</span>
                    <span className="metric-val mono-num">{cand.side_effect_sim.toFixed(4)}</span>
                  </div>
                  {cand.hubness_percentile !== null && cand.hubness_percentile !== undefined && (
                    <div className="compare-metric-row">
                      <span className="metric-label">Hubness</span>
                      <span className="metric-val mono-num">{Number(cand.hubness_percentile).toFixed(1)}%</span>
                    </div>
                  )}
                </div>
              </div>
            ))}
          </div>

          {/* Shared Overlap Analysis if >= 2 items pinned */}
          {pinnedCandidates.length >= 2 && (
            <div className="compare-overlap-section">
              <div className="overlap-title">Overlap Across Compared Candidates:</div>
              <div className="overlap-groups">
                <div className="overlap-group">
                  <span className="overlap-label">Common Targets ({commonTargets.size}):</span>
                  {commonTargets.size > 0 ? (
                    <div className="chip-list">
                      {[...commonTargets].map((t, idx) => (
                        <span key={idx} className="feature-chip chip-target">{t}</span>
                      ))}
                    </div>
                  ) : (
                    <span className="overlap-empty">None in common</span>
                  )}
                </div>
                <div className="overlap-group">
                  <span className="overlap-label">Common Side Effects ({commonSideEffects.size}):</span>
                  {commonSideEffects.size > 0 ? (
                    <div className="chip-list">
                      {[...commonSideEffects].map((se, idx) => (
                        <span key={idx} className="feature-chip chip-se">{se}</span>
                      ))}
                    </div>
                  ) : (
                    <span className="overlap-empty">None in common</span>
                  )}
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
