import { useState, useMemo, useRef } from 'react';
import * as d3Force from 'd3-force';

export function NetworkView({
  queryName,
  candidates = [],
  onSelectCandidate,
}) {
  const [showTable, setShowTable] = useState(false);
  const [selectedNode, setSelectedNode] = useState(null);
  const svgRef = useRef(null);

  const { nodes, links } = useMemo(() => {
    if (!candidates || candidates.length === 0) {
      return { nodes: [], links: [] };
    }

    const width = 800;
    const height = 500;
    const centerX = width / 2;
    const centerY = height / 2;

    // Center node for query drug
    const queryNode = {
      id: 'QUERY',
      name: queryName,
      isQuery: true,
      x: centerX,
      y: centerY,
      fx: centerX,
      fy: centerY,
      radius: 22,
    };

    // Candidate nodes (top 25 candidates for clean graph readability)
    const candNodes = candidates.slice(0, 25).map((c, i) => {
      // Determine dominant component for node coloring
      const s = 0.4 * (c.struct_sim || 0);
      const t = 0.4 * (c.target_sim || 0);
      const e = 0.2 * (c.side_effect_sim || 0);

      let dominant = 'struct';
      if (t >= s && t >= e) dominant = 'target';
      else if (e >= s && e >= t) dominant = 'se';

      return {
        id: `CAND_${i}`,
        name: c.candidate,
        candidateData: c,
        rank: i + 1,
        score: c.score,
        dominant,
        radius: Math.max(10, Math.min(18, 8 + c.score * 18)),
      };
    });

    const allNodes = [queryNode, ...candNodes];
    const allLinks = candNodes.map((cn) => ({
      source: 'QUERY',
      target: cn.id,
      score: cn.score,
    }));

    const simulation = d3Force
      .forceSimulation(allNodes)
      .force(
        'link',
        d3Force
          .forceLink(allLinks)
          .id((d) => d.id)
          .distance((d) => 180 - Math.min(120, d.score * 150))
      )
      .force('charge', d3Force.forceManyBody().strength(-200))
      .force('collide', d3Force.forceCollide().radius((d) => d.radius + 10))
      .force('center', d3Force.forceCenter(centerX, centerY))
      .stop();

    // Run simulation synchronously to achieve stable layout
    for (let i = 0; i < 150; ++i) simulation.tick();

    return { nodes: allNodes, links: allLinks };
  }, [queryName, candidates]);

  const dominantColorMap = {
    struct: 'var(--c-struct, #0077BB)',
    target: 'var(--c-target, #EE7733)',
    se: 'var(--c-se, #009988)',
  };

  return (
    <div className="network-view-panel" role="region" aria-label="Candidate similarity network visualization">
      {/* View Header & Accessibility Toggle */}
      <div className="view-sub-header">
        <div className="view-meta">
          <span className="view-title">Similarity Network</span>
          <span className="view-hint">
            Query at center; edge weight & distance reflect candidate similarity score
          </span>
        </div>
        <div className="view-actions">
          <button
            type="button"
            className="btn-toggle-view"
            onClick={() => setShowTable((prev) => !prev)}
            aria-pressed={showTable}
          >
            {showTable ? 'Switch to Graph View' : 'Switch to Accessible Table'}
          </button>
        </div>
      </div>

      {showTable ? (
        /* Keyboard Accessible Table Alternative */
        <div className="accessible-table-wrapper" tabIndex={0} role="region" aria-label="Network data table alternative">
          <table className="accessible-table">
            <caption className="sr-only">Candidate Similarity Network Data Alternative</caption>
            <thead>
              <tr>
                <th scope="col">Rank</th>
                <th scope="col">Candidate</th>
                <th scope="col">Score</th>
                <th scope="col">Dominant Dimension</th>
                <th scope="col">Structure (40%)</th>
                <th scope="col">Target (40%)</th>
                <th scope="col">Side Effects (20%)</th>
                <th scope="col">Hubness</th>
              </tr>
            </thead>
            <tbody>
              {candidates.slice(0, 25).map((c, i) => {
                const s = 0.4 * (c.struct_sim || 0);
                const t = 0.4 * (c.target_sim || 0);
                const e = 0.2 * (c.side_effect_sim || 0);
                let dom = 'Structure';
                if (t >= s && t >= e) dom = 'Target';
                else if (e >= s && e >= t) dom = 'Side Effects';

                return (
                  <tr key={c.candidate}>
                    <td className="mono-num">#{i + 1}</td>
                    <td><strong>{c.candidate}</strong></td>
                    <td className="mono-num">{c.score.toFixed(4)}</td>
                    <td>{dom}</td>
                    <td className="mono-num">{c.struct_sim.toFixed(4)}</td>
                    <td className="mono-num">{c.target_sim.toFixed(4)}</td>
                    <td className="mono-num">{c.side_effect_sim.toFixed(4)}</td>
                    <td className="mono-num">
                      {c.hubness_percentile !== null && c.hubness_percentile !== undefined
                        ? `${Number(c.hubness_percentile).toFixed(1)}%`
                        : '—'}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      ) : (
        /* SVG Interactive Network Graph */
        <div className="network-svg-container">
          <svg
            ref={svgRef}
            viewBox="0 0 800 500"
            className="network-svg"
            role="img"
            aria-label={`Interactive network graph with center ${queryName} and ${nodes.length - 1} connected candidates`}
          >
            {/* Links */}
            <g className="network-links" aria-hidden="true">
              {links.map((link, idx) => {
                const source = typeof link.source === 'object' ? link.source : nodes.find((n) => n.id === link.source);
                const target = typeof link.target === 'object' ? link.target : nodes.find((n) => n.id === link.target);
                if (!source || !target) return null;

                const strokeWidth = Math.max(1, Math.min(5, link.score * 8));
                return (
                  <line
                    key={idx}
                    x1={source.x}
                    y1={source.y}
                    x2={target.x}
                    y2={target.y}
                    stroke="var(--border, #cbd5e1)"
                    strokeWidth={strokeWidth}
                    strokeOpacity={0.65}
                  />
                );
              })}
            </g>

            {/* Nodes */}
            <g className="network-nodes">
              {nodes.map((node) => {
                const isSelected = selectedNode?.name === node.name;
                if (node.isQuery) {
                  return (
                    <g
                      key={node.id}
                      transform={`translate(${node.x}, ${node.y})`}
                      className="node-query"
                      tabIndex={0}
                      role="button"
                      aria-label={`Query compound: ${node.name}`}
                    >
                      <circle
                        r={node.radius}
                        fill="var(--accent, #2563EB)"
                        stroke="#ffffff"
                        strokeWidth="3"
                        className="pulse-circle"
                      />
                      <text
                        dy=".3em"
                        textAnchor="middle"
                        fill="#ffffff"
                        fontSize="10"
                        fontWeight="700"
                        pointerEvents="none"
                      >
                        QUERY
                      </text>
                      <text
                        dy="2.4em"
                        textAnchor="middle"
                        fill="var(--text-primary, #0f172a)"
                        fontSize="11"
                        fontWeight="700"
                        pointerEvents="none"
                      >
                        {node.name}
                      </text>
                    </g>
                  );
                }

                const color = dominantColorMap[node.dominant] || 'var(--accent)';
                return (
                  <g
                    key={node.id}
                    transform={`translate(${node.x}, ${node.y})`}
                    className={`node-candidate ${isSelected ? 'is-selected' : ''}`}
                    onClick={() => {
                      setSelectedNode(node);
                      if (onSelectCandidate) onSelectCandidate(node.candidateData);
                    }}
                    onKeyDown={(e) => {
                      if (e.key === 'Enter' || e.key === ' ') {
                        e.preventDefault();
                        setSelectedNode(node);
                        if (onSelectCandidate) onSelectCandidate(node.candidateData);
                      }
                    }}
                    tabIndex={0}
                    role="button"
                    aria-label={`Rank ${node.rank}: ${node.name}, score ${node.score.toFixed(4)}`}
                  >
                    <circle
                      r={node.radius}
                      fill={color}
                      stroke={isSelected ? 'var(--text-primary)' : '#ffffff'}
                      strokeWidth={isSelected ? '3' : '1.5'}
                    />
                    <text
                      dy="1.7em"
                      textAnchor="middle"
                      fill="var(--text-secondary, #475569)"
                      fontSize="9"
                      fontWeight="600"
                      className="node-label"
                      pointerEvents="none"
                    >
                      {node.name.length > 10 ? `${node.name.slice(0, 9)}…` : node.name}
                    </text>
                  </g>
                );
              })}
            </g>
          </svg>

          {/* Node Detail Callout on Click */}
          {selectedNode && !selectedNode.isQuery && (
            <div className="network-node-popover" role="dialog" aria-label="Selected Node Details">
              <div className="popover-header">
                <strong>#{selectedNode.rank} {selectedNode.name}</strong>
                <button
                  type="button"
                  className="popover-close"
                  onClick={() => setSelectedNode(null)}
                  aria-label="Close detail popover"
                >
                  ×
                </button>
              </div>
              <div className="popover-body">
                <div>Score: <strong className="mono-num">{selectedNode.score.toFixed(4)}</strong></div>
                <div>Dominant Dimension: <strong>{selectedNode.dominant.toUpperCase()}</strong></div>
                <div className="popover-sims">
                  <span>Struct: {selectedNode.candidateData?.struct_sim?.toFixed(4)}</span>
                  <span>Target: {selectedNode.candidateData?.target_sim?.toFixed(4)}</span>
                  <span>SE: {selectedNode.candidateData?.side_effect_sim?.toFixed(4)}</span>
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
