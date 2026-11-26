export function HubnessBadge({ percentile }) {
  if (percentile === null || percentile === undefined) {
    return null;
  }

  const pct = Number(percentile);
  const isHigh = pct >= 80;
  const isModerate = pct >= 50 && pct < 80;

  const tooltipText = `Hubness: ${pct.toFixed(1)}th percentile. This compound matches many candidates across the database; score is penalised by 1 + (${(pct/100).toFixed(2)})² to reduce hub bias.`;

  return (
    <span
      className={`hubness-badge ${isHigh ? 'hub-high' : isModerate ? 'hub-med' : 'hub-low'}`}
      title={tooltipText}
      aria-label={tooltipText}
      tabIndex={0}
    >
      <span className="hub-dot" aria-hidden="true" />
      Hub: {pct.toFixed(0)}%
    </span>
  );
}
