import { useState, useEffect, Fragment } from 'react';

export function HowToRead({ isOpen, onClose }) {
  const [metrics, setMetrics] = useState(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!isOpen) return;

    fetch('/metrics.json')
      .then((res) => (res.ok ? res.json() : null))
      .then((data) => {
        setMetrics(data);
        setLoading(false);
      })
      .catch((err) => {
        console.error('Failed to load metrics.json:', err);
        setLoading(false);
      });
  }, [isOpen]);

  if (!isOpen) return null;

  return (
    <div className="modal-backdrop" onClick={onClose} role="presentation">
      <div
        className="how-to-read-modal"
        role="dialog"
        aria-modal="true"
        aria-labelledby="how-to-read-title"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="modal-header">
          <div>
            <h2 id="how-to-read-title" className="modal-title">
              How to Read MoleculAR Rankings & Scores
            </h2>
            <p className="modal-subtitle">
              Methodology, score decomposition, and duplicate-neutral benchmark performance.
            </p>
          </div>
          <button
            type="button"
            className="modal-close-btn"
            onClick={onClose}
            aria-label="Close modal"
          >
            ×
          </button>
        </div>

        <div className="modal-content-scroll">
          {/* Section 1: Tri-Modal Similarity Breakdown */}
          <section className="guide-section">
            <h3 className="section-heading">1. Tri-Modal Similarity Components</h3>
            <div className="guide-cards-grid">
              <div className="guide-card card-struct">
                <div className="card-badge badge-struct">40% Weight</div>
                <h4>Structure Similarity (ECFP4)</h4>
                <p>
                  2048-bit chirality-aware circular extended-connectivity fingerprints (ECFP4) comparing 2D chemical scaffolds via Tanimoto similarity coefficient.
                </p>
              </div>

              <div className="guide-card card-target">
                <div className="card-badge badge-target">40% Weight</div>
                <h4>Target Similarity (Jaccard)</h4>
                <p>
                  Overlap of curated high-confidence human protein targets (UniProt accessions) retrieved from ChEMBL 37, scored using Jaccard intersection over union.
                </p>
              </div>

              <div className="guide-card card-se">
                <div className="card-badge badge-se">20% Weight</div>
                <h4>Side-Effect Similarity (Jaccard)</h4>
                <p>
                  Phenotypic profile overlap of clinical adverse drug reactions and phenotypic side effects extracted from SIDER 4.1, computed via Jaccard index.
                </p>
              </div>
            </div>
          </section>

          {/* Section 2: Hubness Bias Correction */}
          <section className="guide-section">
            <h3 className="section-heading">2. Hub Bias Penalty Formula</h3>
            <p className="guide-p">
              In high-dimensional chemical and biological spaces, certain "hub" compounds naturally exhibit high similarity to many unrelated molecules. To suppress false-positive hub bias, raw pre-penalty scores are adjusted by the compound's database hubness percentile rank:
            </p>
            <div className="formula-box">
              <code>
                pre_penalty = 0.40 × Struct + 0.40 × Target + 0.20 × SideEffect<br />
                final_score = pre_penalty / (1 + (hubness_percentile / 100)²)
              </code>
            </div>
            <p className="guide-p-muted">
              Compounds at the 50th percentile receive a 1.25× divisor; compounds at the 80th percentile receive a 1.64× divisor, penalizing ubiquitous matches while preserving genuine multi-modal similarity.
            </p>
            <p className="guide-p-muted" style={{ marginTop: '0.5rem' }}>
              Salt or stereoisomer form badges identify candidate pairs that share the same molecular skeleton (salt, ester or stereoisomer) within the evaluation pool only.
            </p>
          </section>

          {/* Section 3: Empirical Benchmark Performance */}
          <section className="guide-section">
            <h3 className="section-heading">3. Evaluated Performance Metrics</h3>
            <p className="guide-summary-sentence">
              Strong at finding drugs with the same indication; much weaker across drug classes, where a popularity-only baseline can do as well or better.
            </p>

            {loading ? (
              <div className="metrics-loading">Loading benchmark metrics...</div>
            ) : metrics && metrics.protocols ? (
              <div className="metrics-comparison-table-wrapper">
                <table className="metrics-table">
                  <thead>
                    <tr>
                      <th scope="col">Protocol</th>
                      <th scope="col">Candidate Method</th>
                      <th scope="col">Hit@10 [95% CI]</th>
                      <th scope="col">MRR [95% CI]</th>
                    </tr>
                  </thead>
                  <tbody>
                    {Object.values(metrics.protocols).map((p) => (
                      <Fragment key={p.key}>
                        <tr className="protocol-row">
                          <td colSpan="4">
                            <strong>{p.label} (N={p.n_queries})</strong>
                          </td>
                        </tr>
                        {p.methods.map((m) => (
                          <tr key={m.id} className={m.isDefault ? 'default-method-row' : ''}>
                            <td className="protocol-name">{p.key.split('_')[0].toUpperCase()}</td>
                            <td>
                              <strong>{m.label}</strong> {m.isDefault && <span className="active-tag">(Current System)</span>}
                            </td>
                            <td className="mono-num">
                              <span className="font-bold">{(m.hit10.mean * 100).toFixed(1)}%</span>{' '}
                              <span className="ci-text">
                                [{(m.hit10.ci_low * 100).toFixed(1)}%, {(m.hit10.ci_high * 100).toFixed(1)}%]
                              </span>
                            </td>
                            <td className="mono-num">
                              <span className="font-bold">{m.mrr.mean.toFixed(4)}</span>{' '}
                              <span className="ci-text">
                                [{m.mrr.ci_low.toFixed(4)}, {m.mrr.ci_high.toFixed(4)}]
                              </span>
                            </td>
                          </tr>
                        ))}
                      </Fragment>
                    ))}
                  </tbody>
                </table>
              </div>
            ) : (
              <div className="metrics-error">Benchmark metrics unavailable.</div>
            )}
          </section>

          {/* Section 4: Provenance Notice */}
          <div className="provenance-footer">
            <span className="provenance-tag">
              Evaluation Git Tag: <strong>{metrics?.provenance?.git_tag || 'eval-v2'}</strong>
            </span>
            <span className="provenance-sep">•</span>
            <span className="provenance-tag">
              Database SHA-256: <code>{metrics?.provenance?.db_sha256 || 'd88a2333743af2dbe859b3717764b6a645b9786caa6f059cc9558e29c7ff0da7'}</code>
            </span>
          </div>
        </div>
      </div>
    </div>
  );
}
