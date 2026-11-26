export function EmptyState({ mode, error, onRetry }) {
  if (error) {
    return (
      <div className="empty-state error-state" role="alert">
        <div className="empty-icon error-icon">⚠️</div>
        <h3>Unable to complete search</h3>
        <p className="error-message">{error}</p>
        {onRetry && (
          <button className="btn-retry" onClick={onRetry}>
            Retry Search
          </button>
        )}
      </div>
    );
  }

  const isDrug = mode === 'drug';

  return (
    <div className="empty-state" role="status">
      <div className="empty-icon">🔬</div>
      <h3>{isDrug ? 'Search Approved Drugs' : 'Search Target Proteins'}</h3>
      <p>
        {isDrug
          ? 'Enter an approved drug name to rank similar candidates by 2D chemical structure, shared protein targets, and phenotypic side effects.'
          : 'Enter a UniProt accession ID (e.g. P35354 for PTGS2 / COX-2) to identify approved drugs that interact with the protein, ranked by selectivity.'}
      </p>
      <div className="empty-hint">
        Select a quick example chip above or type in the search bar to begin.
      </div>
    </div>
  );
}
