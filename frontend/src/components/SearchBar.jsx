import { useState, useEffect, useRef } from 'react';

const DRUG_EXAMPLES = [
  { label: 'Sildenafil', query: 'Sildenafil' },
  { label: 'Aspirin', query: 'Aspirin' },
  { label: 'Metformin', query: 'Metformin' },
  { label: 'Minoxidil', query: 'Minoxidil' },
  { label: 'Thalidomide', query: 'Thalidomide' },
  { label: 'Duloxetine', query: 'Duloxetine' },
];

const TARGET_EXAMPLES = [
  { label: 'P35354 (PTGS2 / COX-2)', query: 'P35354' },
  { label: 'P28482 (MAPK1 / ERK2)', query: 'P28482' },
  { label: 'P00533 (EGFR)', query: 'P00533' },
  { label: 'P10636 (MAPT / Tau)', query: 'P10636' },
];

export function SearchBar({
  mode,
  query,
  onSearch,
  recentSearches = [],
  onClearRecent,
  loading = false,
}) {
  const [inputValue, setInputValue] = useState(query || '');
  const [prevQuery, setPrevQuery] = useState(query);
  const [showRecent, setShowRecent] = useState(false);
  const [activeIndex, setActiveIndex] = useState(-1);
  const containerRef = useRef(null);

  if (query !== prevQuery) {
    setPrevQuery(query);
    setInputValue(query || '');
  }

  // Filter recent searches for the active mode
  const filteredRecent = recentSearches.filter((r) => r.mode === mode);

  // Close recent searches dropdown on outside click
  useEffect(() => {
    function handleClickOutside(event) {
      if (containerRef.current && !containerRef.current.contains(event.target)) {
        setShowRecent(false);
      }
    }
    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, []);

  const handleSubmit = (e) => {
    e?.preventDefault();
    const trimmed = inputValue.trim();
    if (trimmed) {
      setShowRecent(false);
      onSearch(trimmed);
    }
  };

  const handleSelectRecent = (q) => {
    setInputValue(q);
    setShowRecent(false);
    onSearch(q);
  };

  const handleKeyDown = (e) => {
    if (!showRecent || filteredRecent.length === 0) return;

    if (e.key === 'ArrowDown') {
      e.preventDefault();
      setActiveIndex((prev) => (prev < filteredRecent.length - 1 ? prev + 1 : 0));
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      setActiveIndex((prev) => (prev > 0 ? prev - 1 : filteredRecent.length - 1));
    } else if (e.key === 'Enter' && activeIndex >= 0) {
      e.preventDefault();
      handleSelectRecent(filteredRecent[activeIndex].query);
    } else if (e.key === 'Escape') {
      setShowRecent(false);
    }
  };

  const examples = mode === 'drug' ? DRUG_EXAMPLES : TARGET_EXAMPLES;

  return (
    <div className="search-section" ref={containerRef}>
      <form onSubmit={handleSubmit} className="search-form">
        <div className="search-input-wrapper">
          <svg
            className="search-icon"
            viewBox="0 0 16 16"
            width="16"
            height="16"
            fill="currentColor"
            aria-hidden="true"
          >
            <path d="M11.742 10.344a6.5 6.5 0 1 0-1.397 1.398h-.001c.03.04.062.078.098.115l3.85 3.85a1 1 0 0 0 1.415-1.414l-3.85-3.85a1.007 1.007 0 0 0-.115-.1zM12 6.5a5.5 5.5 0 1 1-11 0 5.5 5.5 0 0 1 11 0z" />
          </svg>

          <input
            type="text"
            className="search-input"
            value={inputValue}
            onChange={(e) => {
              setInputValue(e.target.value);
              setShowRecent(true);
            }}
            onFocus={() => setShowRecent(true)}
            onKeyDown={handleKeyDown}
            placeholder={
              mode === 'drug'
                ? 'Enter approved drug name (e.g., Sildenafil, Aspirin)...'
                : 'Enter UniProt accession ID (e.g., P35354)...'
            }
            aria-label={mode === 'drug' ? 'Drug Search Query' : 'Target UniProt Query'}
            autoComplete="off"
            spellCheck="false"
          />

          {inputValue && (
            <button
              type="button"
              className="search-clear-btn"
              onClick={() => {
                setInputValue('');
                onSearch('');
              }}
              aria-label="Clear search input"
            >
              ×
            </button>
          )}

          <button
            type="submit"
            className="btn-search-submit"
            disabled={loading || !inputValue.trim()}
          >
            {loading ? 'Searching...' : 'Search'}
          </button>
        </div>

        {/* Recent Searches Dropdown */}
        {showRecent && filteredRecent.length > 0 && (
          <div className="recent-dropdown" role="listbox">
            <div className="recent-header">
              <span>Recent Searches</span>
              {onClearRecent && (
                <button
                  type="button"
                  className="btn-clear-recent"
                  onClick={(e) => {
                    e.stopPropagation();
                    onClearRecent();
                  }}
                >
                  Clear all
                </button>
              )}
            </div>
            {filteredRecent.map((item, idx) => (
              <div
                key={idx}
                className={`recent-item ${idx === activeIndex ? 'is-active' : ''}`}
                onClick={() => handleSelectRecent(item.query)}
                role="option"
                aria-selected={idx === activeIndex}
              >
                <span className="recent-query">{item.query}</span>
                <span className="recent-tag">{item.mode === 'drug' ? 'Drug' : 'Target'}</span>
              </div>
            ))}
          </div>
        )}
      </form>

      {/* Quick Example Chips */}
      <div className="example-chips">
        <span className="chips-label">Examples:</span>
        {examples.map((ex, i) => (
          <button
            key={i}
            type="button"
            className="chip-btn"
            onClick={() => {
              setInputValue(ex.query);
              onSearch(ex.query);
            }}
          >
            {ex.label}
          </button>
        ))}
      </div>
    </div>
  );
}
