export function CoverageBadge({ hasTargets, hasSideEffects, size = 'sm' }) {
  if (hasTargets && hasSideEffects) {
    return null;
  }

  const badges = [];

  if (!hasTargets) {
    badges.push({
      key: 'no-targets',
      label: 'No target data',
      detail: 'Score driven by structure and side effects only',
      type: 'warning',
    });
  }

  if (!hasSideEffects) {
    badges.push({
      key: 'no-se',
      label: 'No side-effect data',
      detail: 'Score driven by structure and targets only',
      type: 'warning',
    });
  }

  return (
    <div className={`coverage-badges coverage-badges-${size}`} role="status">
      {badges.map((b) => (
        <span
          key={b.key}
          className="coverage-badge"
          title={`${b.label}: ${b.detail}`}
          aria-label={`${b.label}: ${b.detail}`}
        >
          <svg className="badge-icon" viewBox="0 0 16 16" width="12" height="12" fill="currentColor" aria-hidden="true">
            <path d="M8 1a7 7 0 1 0 0 14A7 7 0 0 0 8 1zm0 10.5a.75.75 0 1 1 0-1.5.75.75 0 0 1 0 1.5zM8.75 4.5v4a.75.75 0 0 1-1.5 0v-4a.75.75 0 0 1 1.5 0z" />
          </svg>
          {b.label}
        </span>
      ))}
    </div>
  );
}
