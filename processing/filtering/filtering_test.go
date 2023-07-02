package Filtering

import (
	np "github.com/geoport/numpy4go/vectors"
	"testing"
)

func TestLinearFilter(t *testing.T) {
	signal1 := np.Loadtxt("..\\..\\TestData\\RSN10_IMPVALL.BG_C-ELC000.txt", 0, false)[0]
	_, filtered1 := FilterSignal(signal1, []float64{10}, 2, "lowpass", "butterworth", 0.005)
	max1 := np.Round(np.Max(filtered1), 3)
	if max1 != 0.029 {
		t.Errorf("Expected 0.03 got %d", max1)
	}

	signal2 := np.Loadtxt("..\\..\\TestData\\RSN10_IMPVALL.BG_C-ELC090.txt", 0, false)[0]
	_, filtered2 := FilterSignal(signal2, []float64{10}, 2, "highpass", "butterworth", 0.005)
	max2 := np.Round(np.Max(filtered2), 3)
	if max2 != 0.006 {
		t.Errorf("Expected 0.006 got %d", max2)
	}

	signal3 := np.Loadtxt("..\\..\\TestData\\RSN11_NWCALIF.AB_B-FRN224.txt", 0, false)[0]
	_, filtered3 := FilterSignal(signal3, []float64{10, 50}, 2, "bandstop", "butterworth", 0.005)
	max3 := np.Round(np.Max(filtered3), 3)
	if max3 != 0.105 {
		t.Errorf("Expected 0.105 got %d", max3)
	}

	signal4 := np.Loadtxt("..\\..\\TestData\\RSN15_KERN_TAF111.txt", 0, false)[0]
	_, filtered4 := FilterSignal(signal4, []float64{10, 50}, 2, "bandpass", "butterworth", 0.005)
	max4 := np.Round(np.Max(filtered4), 3)
	if max4 != 0.078 {
		t.Errorf("Expected 0.078 got %d", max4)
	}
}
