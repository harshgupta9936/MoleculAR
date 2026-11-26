import { useState, useMemo } from 'react';
import * as d3Scale from 'd3-scale';

export function ScatterView({
  candidates = [],
  onSelectCandidate,
}) {
  const [showTable, setShowTable] = useState(false);
  const [hoveredPoint, setHoveredPoint] = useState(null);

  const width = 760;
  const height = 440;
  const margin = { top: 30, right: 30, bottom: 50, left: 60 };

  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  const xScale = useMemo(() => {
    return d3Scale.scaleLinear().domain([0, 1]).range([0, innerWidth]);
  }, [innerWidth]);

  const yScale = useMemo(() => {
    return d3Scale.scaleLinear().domain([0, 1]).range([innerHeight, 0]);
  }, [innerHeight]);

  return (
    <div className="scatter-view-panel" role="region" aria-label="Structure vs Target similarity scatter plot">
      <div className="view-sub-header">
        <div className="view-meta">
          <span className="view-title">Multi-Dimensional Similarity Distribution</span>
          <span className="view-hint">
            X = Structure similarity (ECFP4) • Y = Target similarity (Jaccard) • Point Size = Side Effect similarity
          </span>
        </div>
        <div className="view-actions">
          <button
            type="button"
            className="btn-toggle-view"
            onClick={() => setShowTable((prev) => !prev)}
            aria-pressed={showTable}
          >
            {showTable ? 'Switch to Scatter Plot' : 'Switch to Accessible Table'}
          </button>
        </div>
      </div>

      {showTable ? (
        /* Keyboard Accessible Table Alternative */
        <div className="accessible-table-wrapper" tabIndex={0} role="region" aria-label="Scatter plot data table alternative">
          <table className="accessible-table">
            <caption className="sr-only">Structure and Target Similarity Distribution Data Table</caption>
            <thead>
              <tr>
                <th scope="col">Rank</th>
                <th scope="col">Candidate</th>
                <th scope="col">Score</th>
                <th scope="col">Structure Sim (X)</th>
                <th scope="col">Target Sim (Y)</th>
                <th scope="col">Side Effect Sim (Radius)</th>
                <th scope="col">Hubness Percentile</th>
              </tr>
            </thead>
            <tbody>
              {candidates.map((c, idx) => (
                <tr key={c.candidate}>
                  <td className="mono-num">#{idx + 1}</td>
                  <td><strong>{c.candidate}</strong></td>
                  <td className="mono-num">{c.score.toFixed(4)}</td>
                  <td className="mono-num">{c.struct_sim.toFixed(4)}</td>
                  <td className="mono-num">{c.target_sim.toFixed(4)}</td>
                  <td className="mono-num">{c.side_effect_sim.toFixed(4)}</td>
                  <td className="mono-num">
                    {c.hubness_percentile !== null && c.hubness_percentile !== undefined
                      ? `${Number(c.hubness_percentile).toFixed(1)}%`
                      : '—'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : (
        /* SVG Interactive Scatter Plot */
        <div className="scatter-svg-container">
          <svg
            viewBox={`0 0 ${width} ${height}`}
            className="scatter-svg"
            role="img"
            aria-label="Scatter plot showing structure similarity vs target similarity for all candidates"
          >
            <g transform={`translate(${margin.left}, ${margin.top})`}>
              {/* Grid Lines */}
              {[0, 0.2, 0.4, 0.6, 0.8, 1.0].map((tick) => (
                <g key={tick} className="grid-tick">
                  {/* Horizontal grid line */}
                  <line
                    x1={0}
                    y1={yScale(tick)}
                    x2={innerWidth}
                    y2={yScale(tick)}
                    stroke="var(--border, #e2e8f0)"
                    strokeDasharray="3 3"
                  />
                  <text
                    x={-10}
                    y={yScale(tick)}
                    dy=".32em"
                    textAnchor="end"
                    fontSize="10"
                    fill="var(--text-muted, #94a3b8)"
                    className="mono-num"
                  >
                    {tick.toFixed(1)}
                  </text>

                  {/* Vertical grid line */}
                  <line
                    x1={xScale(tick)}
                    y1={0}
                    x2={xScale(tick)}
                    y2={innerHeight}
                    stroke="var(--border, #e2e8f0)"
                    strokeDasharray="3 3"
                  />
                  <text
                    x={xScale(tick)}
                    y={innerHeight + 18}
                    textAnchor="middle"
                    fontSize="10"
                    fill="var(--text-muted, #94a3b8)"
                    className="mono-num"
                  >
                    {tick.toFixed(1)}
                  </text>
                </g>
              ))}

              {/* Axes lines */}
              <line x1={0} y1={innerHeight} x2={innerWidth} y2={innerHeight} stroke="var(--border, #94a3b8)" strokeWidth="1.5" />
              <line x1={0} y1={0} x2={0} y2={innerHeight} stroke="var(--border, #94a3b8)" strokeWidth="1.5" />

              {/* Axis Labels */}
              <text
                x={innerWidth / 2}
                y={innerHeight + 40}
                textAnchor="middle"
                fontSize="12"
                fontWeight="600"
                fill="var(--text-secondary, #475569)"
              >
                Structure Similarity (ECFP4) →
              </text>
              <text
                transform="rotate(-90)"
                x={-innerHeight / 2}
                y={-42}
                textAnchor="middle"
                fontSize="12"
                fontWeight="600"
                fill="var(--text-secondary, #475569)"
              >
                Target Similarity (Jaccard) →
              </text>

              {/* Data Points */}
              {candidates.map((c, i) => {
                const cx = xScale(c.struct_sim || 0);
                const cy = yScale(c.target_sim || 0);
                const r = Math.max(5, Math.min(18, 5 + (c.side_effect_sim || 0) * 16));
                const isHovered = hoveredPoint?.candidate === c.candidate;

                return (
                  <circle
                    key={c.candidate}
                    cx={cx}
                    cy={cy}
                    r={isHovered ? r + 3 : r}
                    fill={isHovered ? 'var(--accent, #2563EB)' : 'var(--c-struct, #0077BB)'}
                    fillOpacity={0.75}
                    stroke={isHovered ? '#ffffff' : 'var(--surface, #ffffff)'}
                    strokeWidth={isHovered ? 2.5 : 1.5}
                    className="scatter-point"
                    tabIndex={0}
                    role="button"
                    aria-label={`Rank ${i + 1}: ${c.candidate}, Struct: ${c.struct_sim.toFixed(4)}, Target: ${c.target_sim.toFixed(4)}, Side Effect: ${c.side_effect_sim.toFixed(4)}`}
                    onMouseEnter={() => setHoveredPoint({ ...c, rank: i + 1, cx, cy })}
                    onMouseLeave={() => setHoveredPoint(null)}
                    onFocus={() => setHoveredPoint({ ...c, rank: i + 1, cx, cy })}
                    onBlur={() => setHoveredPoint(null)}
                    onClick={() => {
                      if (onSelectCandidate) onSelectCandidate(c);
                    }}
                  />
                );
              })}
            </g>
          </svg>

          {/* Interactive Tooltip */}
          {hoveredPoint && (
            <div
              className="scatter-tooltip"
              style={{
                left: `${margin.left + hoveredPoint.cx}px`,
                top: `${margin.top + hoveredPoint.cy - 12}px`,
              }}
              role="tooltip"
            >
              <div className="tooltip-title">
                #{hoveredPoint.rank} {hoveredPoint.candidate}
              </div>
              <div className="tooltip-stat">Score: <strong className="mono-num">{hoveredPoint.score.toFixed(4)}</strong></div>
              <div className="tooltip-stat">Structure: <span className="mono-num">{hoveredPoint.struct_sim.toFixed(4)}</span></div>
              <div className="tooltip-stat">Target: <span className="mono-num">{hoveredPoint.target_sim.toFixed(4)}</span></div>
              <div className="tooltip-stat">Side Effects: <span className="mono-num">{hoveredPoint.side_effect_sim.toFixed(4)}</span></div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
