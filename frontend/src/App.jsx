import { useState } from 'react';
import { useUrlState } from './hooks/useUrlState';
import { useSearch } from './hooks/useSearch';
import { SearchBar } from './components/SearchBar';
import { ResultList } from './components/ResultList';
import { HowToRead } from './components/HowToRead';

function App() {
  const [urlState, setUrlState] = useUrlState();
  const { mode, q, topN } = urlState;
  const [showHowToRead, setShowHowToRead] = useState(false);

  const {
    data,
    loading,
    error,
    recentSearches,
    clearRecent,
    retry,
  } = useSearch(mode, q, topN);

  const handleModeChange = (newMode) => {
    if (newMode !== mode) {
      setUrlState({ mode: newMode, q: '' });
    }
  };

  const handleSearch = (newQuery) => {
    setUrlState({ q: newQuery });
  };

  return (
    <div className="app-container">
      {/* Top Header */}
      <header className="header">
        <div className="header-brand">
          <div className="brand-logo" aria-hidden="true">
            <svg viewBox="0 0 24 24" width="24" height="24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="3" />
              <circle cx="19" cy="5" r="2" />
              <circle cx="5" cy="19" r="2" />
              <path d="M10.5 13.5L6.5 17.5M13.5 10.5L17.5 6.5" />
            </svg>
          </div>
          <div>
            <h1 className="brand-title">MoleculAR</h1>
            <p className="brand-subtitle">
              {mode === 'drug'
                ? 'Rank approved drugs by similarity of structure, protein targets and side effects.'
                : 'List approved drugs that act on a given protein (UniProt ID), ranked by selectivity.'}
            </p>
          </div>
        </div>

        {/* Mode Navigation Tabs */}
        <nav className="mode-toggle" aria-label="Search Mode Navigation">
          <button
            type="button"
            className={`mode-btn ${mode === 'drug' ? 'active' : ''}`}
            onClick={() => handleModeChange('drug')}
            aria-selected={mode === 'drug'}
            role="tab"
          >
            Drug Search
          </button>
          <button
            type="button"
            className={`mode-btn ${mode === 'target' ? 'active' : ''}`}
            onClick={() => handleModeChange('target')}
            aria-selected={mode === 'target'}
            role="tab"
          >
            Target Search
          </button>
        </nav>
      </header>

      {/* Main Search and Results Area */}
      <main className="main-content">
        <SearchBar
          mode={mode}
          query={q}
          onSearch={handleSearch}
          recentSearches={recentSearches}
          onClearRecent={clearRecent}
          loading={loading}
        />

        <ResultList
          mode={mode}
          query={q}
          data={data}
          loading={loading}
          error={error}
          onRetry={retry}
          onOpenHowToRead={() => setShowHowToRead(true)}
        />
      </main>

      {/* How to Read Modal */}
      <HowToRead
        isOpen={showHowToRead}
        onClose={() => setShowHowToRead(false)}
      />

      {/* Persistent Research Tool Footer */}
      <footer className="footer">
        <p>Research tool. Rankings show similarity, not proven treatments. Not medical advice.</p>
        <p>Research demo, non-commercial. Data: ChEMBL 37 (CC BY-SA 3.0), DrugCentral (CC BY-SA 4.0), SIDER 4.1 (CC BY-NC-SA 4.0).</p>
      </footer>
    </div>
  );
}

export default App;
