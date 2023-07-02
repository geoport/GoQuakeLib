package Scaling

import (
	td "github.com/geoport/GoQuakeLib/TestData"
	ds "github.com/geoport/GoQuakeLib/design_spectrums"
	rs "github.com/geoport/GoQuakeLib/response_spectra"
	"testing"

	np "github.com/geoport/numpy4go/vectors"
)

var (
	targetSpectrum, periods = ds.GetSpectrumByTBDY(0.02, 8, 0.4, 0.2, false)
	responseSpectrum        = rs.ResponseSpectra(td.TestMotion["Accelerations"].([]float64), 0.05, periods, 0.05)
	spectralAccelerations   = [][]float64{responseSpectrum.SpectralAccelerations}
)

func TestScaleByMSE(t *testing.T) {
	expectedScaleFactor := 0.631
	actualScaleFactor, _ := ScaleByMSE(periods, targetSpectrum, spectralAccelerations, 1, 6)

	if expectedScaleFactor != np.Round(actualScaleFactor[0], 3) {
		t.Errorf("Expected %f, got %f", expectedScaleFactor, actualScaleFactor)
	}
}

func TestScaleByPeriodRange(t *testing.T) {
	expectedScaleFactor := 1.476
	actualScaleFactor, _ := ScaleByPeriodRange(periods, targetSpectrum, spectralAccelerations, 1, 6)

	if expectedScaleFactor != np.Round(actualScaleFactor[0], 3) {
		t.Errorf("Expected %f, got %f", expectedScaleFactor, actualScaleFactor)
	}
}

func TestScaleBySinglePeriod(t *testing.T) {
	expextedScaleFactor := 0.701
	actualScaleFactor, _ := ScaleBySinglePeriod(periods, targetSpectrum, spectralAccelerations, 3)

	if expextedScaleFactor != np.Round(actualScaleFactor[0], 3) {
		t.Errorf("Expected %f, got %f", expextedScaleFactor, actualScaleFactor)
	}
}
