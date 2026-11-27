import { useState } from 'react';
import { API_BASE_URL } from '../utils/api';

/**
 * ExpandedDetail view for a candidate.
 * Shows:
 *   - 2D structure SVG fetched from /api/structure?name=<drug_name>
 *   - Shared protein targets chips
 *   - Shared side-effect chips
 *   - Component similarity breakdown
 *
 * NOTE: The Known Indication column is strictly hidden. No indication data is rendered.
 */
export function ExpandedDetail({ candidate }) {
  const [imgError, setImgError] = useState(false);

  const structUrl = `${API_BASE_URL}/api/structure?name=${encodeURIComponent(candidate.candidate)}`;

  return (
    <div className="expanded-detail">
      <div className="detail-grid">
        {/* Structure Preview Card */}
        <div className="structure-card">
          <div className="card-header-sm">2D Chemical Structure</div>
          <div className="structure-img-wrapper">
            {!imgError ? (
              <img
                src={structUrl}
                alt={`2D structure of ${candidate.candidate}`}
                className="structure-svg"
                loading="lazy"
                onError={() => setImgError(true)}
              />
            ) : (
              <div className="structure-placeholder">Structure preview not available</div>
            )}
          </div>
          <div className="structure-caption">2D Structure: {candidate.candidate}</div>
        </div>

        {/* Similarities & Shared Biological Features */}
        <div className="features-card">
          {/* Detailed Similarities */}
          <div className="card-header-sm">Similarity Breakdown</div>
          <div className="sim-bars">
            <div className="sim-item">
              <div className="sim-label">
                <span>Structure Similarity (ECFP4)</span>
                <span className="mono-num">{candidate.struct_sim.toFixed(4)}</span>
              </div>
              <div className="sim-track">
                <div
                  className="sim-fill sim-struct"
                  style={{ width: `${Math.min(100, candidate.struct_sim * 100)}%` }}
                />
              </div>
            </div>

            <div className="sim-item">
              <div className="sim-label">
                <span>Target Similarity (Jaccard)</span>
                <span className="mono-num">{candidate.target_sim.toFixed(4)}</span>
              </div>
              <div className="sim-track">
                <div
                  className="sim-fill sim-target"
                  style={{ width: `${Math.min(100, candidate.target_sim * 100)}%` }}
                />
              </div>
            </div>

            <div className="sim-item">
              <div className="sim-label">
                <span>Side-Effect Similarity (Jaccard)</span>
                <span className="mono-num">{candidate.side_effect_sim.toFixed(4)}</span>
              </div>
              <div className="sim-track">
                <div
                  className="sim-fill sim-se"
                  style={{ width: `${Math.min(100, candidate.side_effect_sim * 100)}%` }}
                />
              </div>
            </div>
          </div>

          {/* Shared Targets */}
          <div className="shared-section">
            <div className="section-title">
              Shared Targets ({candidate.shared_targets?.length || 0})
            </div>
            {candidate.shared_targets && candidate.shared_targets.length > 0 ? (
              <div className="chip-list">
                {candidate.shared_targets.map((tgt, i) => (
                  <span key={i} className="feature-chip chip-target" title={`Shared target: ${tgt}`}>
                    {tgt}
                  </span>
                ))}
              </div>
            ) : (
              <div className="empty-features">No overlapping protein targets recorded.</div>
            )}
          </div>

          {/* Shared Side Effects */}
          <div className="shared-section">
            <div className="section-title">
              Shared Side Effects ({candidate.shared_side_effects?.length || 0})
            </div>
            {candidate.shared_side_effects && candidate.shared_side_effects.length > 0 ? (
              <div className="chip-list">
                {candidate.shared_side_effects.map((se, i) => (
                  <span key={i} className="feature-chip chip-se" title={`Shared side effect: ${se}`}>
                    {se}
                  </span>
                ))}
              </div>
            ) : (
              <div className="empty-features">No overlapping side-effect profiles recorded.</div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}
