package fourier_spectrum

import (
	td "github.com/geoport/GoQuakeLib/TestData"
	"math"
	"testing"

	np "github.com/geoport/numpy4go/vectors"
)

func compareOutput(actual, expected float64) bool {
	actual = np.Round(actual, 5).(float64)
	expected = np.Round(expected, 5).(float64)
	diff := math.Abs(actual - expected)
	return np.Round(diff, 1).(float64) < 1e-10
}

func TestFourierSpectrum(t *testing.T) {
	// Run analysis
	acc := td.TestMotion["Accelerations"].([]float64)
	timeStep := td.TestMotion["TimeStep"].(float64)
	frequencies, fourierAmplitudes, powerAmplitudes := FourierSpectrum(acc, timeStep)

	expectedMinFrequency := 0.0
	expectedMaxFrequency := 9.994
	expectedMinFourierAmplitude := 0.0
	expectedMaxFourierAmplitude := 0.001
	expectedMinPowerAmplitude := 0.0
	expectedMaxPowerAmplitude := 0.0

	if !compareOutput(np.Min(frequencies), expectedMinFrequency) {
		t.Errorf("Expected min(Frequency) = %f, got %f", expectedMinFrequency, np.Min(frequencies))
	}

	if !compareOutput(np.Max(frequencies), expectedMaxFrequency) {
		t.Errorf("Expected max(Frequency) = %f, got %f", expectedMaxFrequency, np.Round(np.Max(frequencies), 3))
	}

	if !compareOutput(np.Min(fourierAmplitudes), expectedMinFourierAmplitude) {
		t.Errorf("Expected min(FourierAmplitude) = %f, got %f", expectedMinFourierAmplitude, np.Min(fourierAmplitudes))
	}

	if !compareOutput(np.Max(fourierAmplitudes), expectedMaxFourierAmplitude) {
		t.Errorf(
			"Expected max(FourierAmplitude) = %f, got %f", expectedMaxFourierAmplitude,
			np.Round(np.Max(fourierAmplitudes), 3),
		)
	}

	if !compareOutput(np.Min(powerAmplitudes), expectedMinPowerAmplitude) {
		t.Errorf("Expected min(powerAmplitude) = %f, got %f", expectedMinPowerAmplitude, np.Min(powerAmplitudes))
	}

	if !compareOutput(np.Max(powerAmplitudes), expectedMaxPowerAmplitude) {
		t.Errorf(
			"Expected max(powerAmplitude) = %f, got %f", expectedMaxPowerAmplitude,
			np.Round(np.Max(powerAmplitudes), 3),
		)
	}

}
