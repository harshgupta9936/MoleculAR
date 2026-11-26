export function SkeletonRow() {
  return (
    <div className="result-row skeleton-row" aria-hidden="true">
      <div className="col-rank">
        <div className="skeleton-box skeleton-rank" />
      </div>
      <div className="col-name">
        <div className="skeleton-box skeleton-title" />
        <div className="skeleton-box skeleton-sub" />
      </div>
      <div className="col-score">
        <div className="skeleton-box skeleton-score" />
      </div>
      <div className="col-bar">
        <div className="skeleton-box skeleton-bar" />
      </div>
      <div className="col-expand">
        <div className="skeleton-box skeleton-btn" />
      </div>
    </div>
  );
}

export function SkeletonTable({ count = 10 }) {
  return (
    <div className="results-container skeleton-container" aria-label="Loading candidates...">
      {Array.from({ length: count }).map((_, i) => (
        <SkeletonRow key={i} />
      ))}
    </div>
  );
}
