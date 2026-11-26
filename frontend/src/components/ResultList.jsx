import { useState, useMemo } from 'react';
import { ResultRow } from './ResultRow';
import { SkeletonTable } from './SkeletonRow';
import { EmptyState } from './EmptyState';
import { ExportButton } from './ExportButton';
import { CoverageBadge } from './CoverageBadge';
import { HubnessBadge } from './HubnessBadge';
import { WeightSliders, DEFAULT_WEIGHTS } from './WeightSliders';
import { CompareTray } from './CompareTray';
import { NetworkView } from './NetworkView';
import { ScatterView } from './ScatterView';

const PAGE_SIZE = 25;

export function ResultList({
  mode,
  query,
  data,
  loading,
  error,
  onRetry,
  onOpenHowToRead,
}) {
  const [currentPage, setCurrentPage] = useState(1);
  const [prevKey, setPrevKey] = useState(`${mode}:${query}`);
  const [activeTab, setActiveTab] = useState('list'); // 'list' | 'network' | 'scatter'
  const [showSliders, setShowSliders] = useState(false);
  const [weights, setWeights] = useState(DEFAULT_WEIGHTS);
  const [pinnedCandidates, setPinnedCandidates] = useState([]);

  // Reset pagination and weights on new search
  const currentKey = `${mode}:${query}`;
  if (currentKey !== prevKey) {
    setPrevKey(currentKey);
    setCurrentPage(1);
    setWeights(DEFAULT_WEIGHTS);
  }

  // Client-side exploratory re-ranking if weights or penalty differ from defaults
  const candidates = useMemo(() => {
    const rawList = data?.results;
    if (!rawList || rawList.length === 0) return [];

    const isDefault =
      Math.abs(weights.wStruct - DEFAULT_WEIGHTS.wStruct) < 1e-4 &&
      Math.abs(weights.wTarget - DEFAULT_WEIGHTS.wTarget) < 1e-4 &&
      Math.abs(weights.wSideEffect - DEFAULT_WEIGHTS.wSideEffect) < 1e-4 &&
      weights.useHubCorrection === DEFAULT_WEIGHTS.useHubCorrection;

    if (isDefault) return rawList;

    return [...rawList]
      .map((c) => {
        const s = c.struct_sim || 0;
        const t = c.target_sim || 0;
        const e = c.side_effect_sim || 0;
        const pre = weights.wStruct * s + weights.wTarget * t + weights.wSideEffect * e;

        const hub = c.hubness_percentile;
        const pen =
          weights.useHubCorrection && hub !== null && hub !== undefined
            ? 1.0 + Math.pow(hub / 100.0, 2)
            : 1.0;

        return {
          ...c,
          score: pre / pen,
        };
      })
      .sort((a, b) => b.score - a.score);
  }, [data?.results, weights]);

  const handleToggleCompare = (candidate) => {
    setPinnedCandidates((prev) => {
      const exists = prev.some((p) => p.candidate === candidate.candidate);
      if (exists) {
        return prev.filter((p) => p.candidate !== candidate.candidate);
      }
      if (prev.length >= 3) return prev;
      return [...prev, candidate];
    });
  };

  const handleRemoveCompare = (candName) => {
    setPinnedCandidates((prev) => prev.filter((p) => p.candidate !== candName));
  };

  const handleClearCompare = () => {
    setPinnedCandidates([]);
  };

  if (loading) {
    return <SkeletonTable count={10} />;
  }

  if (error || !data || !query) {
    return <EmptyState mode={mode} error={error} onRetry={onRetry} />;
  }

  const totalCount = candidates.length;

  if (totalCount === 0) {
    return (
      <div className="empty-state">
        <div className="empty-icon">🔍</div>
        <h3>No candidates found</h3>
        <p>No similar approved drugs were found matching "{query}".</p>
      </div>
    );
  }

  const totalPages = Math.max(1, Math.ceil(totalCount / PAGE_SIZE));
  const pageStart = (currentPage - 1) * PAGE_SIZE;
  const pageEnd = Math.min(pageStart + PAGE_SIZE, totalCount);
  const pagedCandidates = candidates.slice(pageStart, pageEnd);

  return (
    <div className="results-wrapper">
      {/* Query Drug / Target Summary Banner */}
      <div className="query-summary-banner">
        <div className="query-meta-left">
          <span className="query-badge-type">{mode === 'drug' ? 'Query Drug' : 'Target Protein'}</span>
          <h2 className="query-title">{mode === 'drug' ? data.query_drug : data.query_target}</h2>
          {mode === 'drug' && data.query_hubness_percentile !== null && data.query_hubness_percentile !== undefined && (
            <HubnessBadge percentile={data.query_hubness_percentile} />
          )}
          <CoverageBadge
            hasTargets={data.query_has_targets}
            hasSideEffects={data.query_has_side_effects}
            size="md"
          />
        </div>

        <div className="query-meta-right">
          <button
            type="button"
            className="btn-how-to-read"
            onClick={onOpenHowToRead}
            title="Read methodology, scoring formula & benchmark performance"
          >
            ℹ️ How to Read
          </button>
          <ExportButton
            mode={mode}
            query={query}
            queryEnvelope={data}
            candidates={candidates}
          />
        </div>
      </div>

      {/* Navigation & View Selection Toolbar */}
      <div className="view-toolbar">
        <div className="tab-group" role="tablist" aria-label="Candidate Visualization Views">
          <button
            type="button"
            role="tab"
            aria-selected={activeTab === 'list'}
            className={`tab-btn ${activeTab === 'list' ? 'is-active' : ''}`}
            onClick={() => setActiveTab('list')}
          >
            📋 Candidate List
          </button>
          <button
            type="button"
            role="tab"
            aria-selected={activeTab === 'network'}
            className={`tab-btn ${activeTab === 'network' ? 'is-active' : ''}`}
            onClick={() => setActiveTab('network')}
          >
            🕸️ Network View
          </button>
          <button
            type="button"
            role="tab"
            aria-selected={activeTab === 'scatter'}
            className={`tab-btn ${activeTab === 'scatter' ? 'is-active' : ''}`}
            onClick={() => setActiveTab('scatter')}
          >
            📊 Scatter Distribution
          </button>
        </div>

        <button
          type="button"
          className={`btn-toggle-sliders ${showSliders ? 'is-open' : ''}`}
          onClick={() => setShowSliders((prev) => !prev)}
          aria-expanded={showSliders}
        >
          ⚙️ {showSliders ? 'Hide Re-Ranking Sliders' : 'Exploratory Sliders'}
        </button>
      </div>

      {/* Exploratory Weight Sliders Panel */}
      {showSliders && (
        <WeightSliders
          weights={weights}
          onChange={setWeights}
          onReset={() => setWeights(DEFAULT_WEIGHTS)}
          candidateCount={totalCount}
        />
      )}

      {/* Active Tab View */}
      {activeTab === 'network' ? (
        <NetworkView
          queryName={mode === 'drug' ? data.query_drug : data.query_target}
          candidates={candidates}
          onSelectCandidate={handleToggleCompare}
        />
      ) : activeTab === 'scatter' ? (
        <ScatterView
          candidates={candidates}
          onSelectCandidate={handleToggleCompare}
        />
      ) : (
        /* List View */
        <>
          <div className="results-header-bar">
            <div className="results-count">
              Showing <strong>{pageStart + 1}–{pageEnd}</strong> of <strong>{totalCount}</strong> candidates
            </div>

            <div className="bar-legend">
              <span className="legend-item">
                <span className="legend-box legend-struct" /> Structure (40%)
              </span>
              <span className="legend-item">
                <span className="legend-box legend-target" /> Target (40%)
              </span>
              <span className="legend-item">
                <span className="legend-box legend-se" /> Side Effects (20%)
              </span>
              <span className="legend-item">
                <span className="legend-line" /> Final Score
              </span>
            </div>
          </div>

          <div className="table-header-row">
            <div className="col-compare">Pin</div>
            <div className="col-rank">Rank</div>
            <div className="col-name">Candidate Drug</div>
            <div className="col-score">Score</div>
            <div className="col-bar">Decomposition & Penalty</div>
            <div className="col-expand">Detail</div>
          </div>

          <div className="results-rows">
            {pagedCandidates.map((cand, idx) => (
              <ResultRow
                key={cand.candidate || idx}
                candidate={cand}
                rank={pageStart + idx + 1}
                isCompared={pinnedCandidates.some((p) => p.candidate === cand.candidate)}
                onToggleCompare={handleToggleCompare}
                canCompare={pinnedCandidates.length < 3}
              />
            ))}
          </div>

          {totalPages > 1 && (
            <div className="pagination-bar">
              <button
                type="button"
                className="pagination-btn"
                disabled={currentPage === 1}
                onClick={() => setCurrentPage((p) => Math.max(1, p - 1))}
              >
                ← Previous
              </button>
              <span className="pagination-status">
                Page <strong>{currentPage}</strong> of <strong>{totalPages}</strong>
              </span>
              <button
                type="button"
                className="pagination-btn"
                disabled={currentPage === totalPages}
                onClick={() => setCurrentPage((p) => Math.min(totalPages, p + 1))}
              >
                Next →
              </button>
            </div>
          )}
        </>
      )}

      {/* Floating / Docked Compare Tray */}
      <CompareTray
        pinnedCandidates={pinnedCandidates}
        onRemove={handleRemoveCompare}
        onClear={handleClearCompare}
      />
    </div>
  );
}
