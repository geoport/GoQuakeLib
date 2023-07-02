package design_spectrums

import (
	td "github.com/geoport/GoQuakeLib/TestData"
	np "github.com/geoport/numpy4go/vectors"
	"reflect"
	"testing"
)

func TestInsertPeriods(t *testing.T) {
	periods := []float64{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
	startPeriod := 0.15
	endPeriod := 0.95
	expected := []float64{0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0}
	actual := InsertPeriods(periods, startPeriod, endPeriod)
	if reflect.DeepEqual(actual, expected) == false {
		t.Errorf("InsertPeriods() = %v, want %v", actual, expected)
	}
}

func TestGetSpectrumByTBDY(t *testing.T) {
	expected := np.Round(td.TestTBDYSpectrum, 5).([]float64)
	actual, _ := GetSpectrumByTBDY(0.02, 8, 0.4, 0.2, false)
	diff := np.Round(np.Substract(actual, expected), 3).([]float64)
	if reflect.DeepEqual(diff, make([]float64, len(diff))) == false {
		t.Errorf("GetSpectrumByTBDY() = %v, want %v", actual, expected)
	}
}
