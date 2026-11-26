export default function ExploratoryTargetMode() {
  return (
    <div className="mode-container exploratory-mode">
      <div className="badge-container">
        <span className="badge exploratory-badge">⚠ Exploratory Mode</span>
      </div>
      
      <p className="description">
        Search by novel disease target or family. This mode uses untested connectivity mapping 
        and broad target-family inferences. Hit rates are unknown.
      </p>

      <form className="search-form" onSubmit={(e) => e.preventDefault()}>
        <input 
          type="text" 
          disabled
          placeholder="e.g., Viral RNA polymerase (Coming in Phase 9)"
          className="search-input disabled"
        />
        <button disabled className="search-btn disabled">
          Search Targets
        </button>
      </form>

      <div className="pending-message">
        <h3>Backend pending Phase 9 (Target-Family Mapping)</h3>
        <p>This endpoint (/api/target-search) will be enabled once the manual target-family dataset expansion is complete.</p>
      </div>
    </div>
  );
}
