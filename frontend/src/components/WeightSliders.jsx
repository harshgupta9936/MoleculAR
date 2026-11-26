import { useMemo } from 'react';

const DEFAULT_WEIGHTS = {
  wStruct: 0.4,
  wTarget: 0.4,
  wSideEffect: 0.2,
  useHubCorrection: true,
};

export function WeightSliders({
  weights,
  onChange,
  onReset,
  candidateCount,
}) {
  const isDefault = useMemo(() => {
    return (
      Math.abs(weights.wStruct - DEFAULT_WEIGHTS.wStruct) < 1e-4 &&
      Math.abs(weights.wTarget - DEFAULT_WEIGHTS.wTarget) < 1e-4 &&
      Math.abs(weights.wSideEffect - DEFAULT_WEIGHTS.wSideEffect) < 1e-4 &&
      weights.useHubCorrection === DEFAULT_WEIGHTS.useHubCorrection
    );
  }, [weights]);

  const handleSliderChange = (key, val) => {
    onChange({
      ...weights,
      [key]: parseFloat(val),
    });
  };

  const handleToggleHub = () => {
    onChange({
      ...weights,
      useHubCorrection: !weights.useHubCorrection,
    });
  };

  return (
    <div className="exploratory-panel" role="region" aria-label="Exploratory weight sliders">
      <div className="exploratory-header">
        <div className="exploratory-title-row">
          <span className="exploratory-title">Exploratory Re-Ranking</span>
          {!isDefault && (
            <span className="exploratory-warning-badge" role="status">
              Exploratory setting, not the evaluated default
            </span>
          )}
        </div>
        <p className="exploratory-notice">
          Re-orders only these {candidateCount} returned candidates, not the full database
        </p>
      </div>

      <div className="sliders-grid">
        {/* Structure Slider */}
        <div className="slider-control">
          <div className="slider-label-row">
            <span className="slider-label">
              <span className="slider-indicator struct-indicator" /> Structure Weight
            </span>
            <span className="slider-value mono-num">{weights.wStruct.toFixed(2)}</span>
          </div>
          <input
            type="range"
            min="0"
            max="1"
            step="0.05"
            value={weights.wStruct}
            onChange={(e) => handleSliderChange('wStruct', e.target.value)}
            className="range-input range-struct"
            aria-label="Structure weight"
          />
        </div>

        {/* Target Slider */}
        <div className="slider-control">
          <div className="slider-label-row">
            <span className="slider-label">
              <span className="slider-indicator target-indicator" /> Target Weight
            </span>
            <span className="slider-value mono-num">{weights.wTarget.toFixed(2)}</span>
          </div>
          <input
            type="range"
            min="0"
            max="1"
            step="0.05"
            value={weights.wTarget}
            onChange={(e) => handleSliderChange('wTarget', e.target.value)}
            className="range-input range-target"
            aria-label="Target weight"
          />
        </div>

        {/* Side Effect Slider */}
        <div className="slider-control">
          <div className="slider-label-row">
            <span className="slider-label">
              <span className="slider-indicator se-indicator" /> Side Effect Weight
            </span>
            <span className="slider-value mono-num">{weights.wSideEffect.toFixed(2)}</span>
          </div>
          <input
            type="range"
            min="0"
            max="1"
            step="0.05"
            value={weights.wSideEffect}
            onChange={(e) => handleSliderChange('wSideEffect', e.target.value)}
            className="range-input range-se"
            aria-label="Side effect weight"
          />
        </div>

        {/* Hub Penalty Toggle & Reset */}
        <div className="slider-actions">
          <label className="hub-toggle-label">
            <input
              type="checkbox"
              checked={weights.useHubCorrection}
              onChange={handleToggleHub}
              className="hub-checkbox"
            />
            <span>Apply Hub Bias Penalty</span>
          </label>

          {!isDefault && (
            <button
              type="button"
              onClick={onReset}
              className="btn-reset-weights"
              title="Reset weights to 0.40 / 0.40 / 0.20 with hub penalty enabled"
            >
              Reset to Evaluated Default
            </button>
          )}
        </div>
      </div>
    </div>
  );
}

export { DEFAULT_WEIGHTS };
