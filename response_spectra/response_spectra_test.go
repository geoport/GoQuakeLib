package response_spectra

import (
	td "github.com/geoport/GoQuakeLib/TestData"
	"reflect"
	"testing"

	np "github.com/geoport/numpy4go/vectors"
)

func TestGetTimeSeries(t *testing.T) {
	expectedXA1 := np.Round(td.TestTimeSeries["xa_0:10"], 8)
	expectedXV1 := np.Round(td.TestTimeSeries["xv_0:10"], 8)
	expectedXD1 := np.Round(td.TestTimeSeries["xd_0:10"], 8)
	expectedXA2 := np.Round(td.TestTimeSeries["xa_20:10"], 8)
	expectedXV2 := np.Round(td.TestTimeSeries["xv_20:10"], 8)
	expectedXD2 := np.Round(td.TestTimeSeries["xd_20:10"], 8)
	testAcceleration := td.TestMotion["Accelerations"].([]float64)
	output := GetTimeSeries(
		td.TestConst, td.TestOmega2, len(testAcceleration), len(td.TestPeriods), testAcceleration, 0.005,
	)
	outputXA1 := np.Round(output["xa"][0][:10], 8)
	outputXV1 := np.Round(output["xv"][0][:10], 8)
	outputXD1 := np.Round(output["xd"][0][:10], 8)

	outputXA2 := np.Round(output["xa"][20][:10], 8)
	outputXV2 := np.Round(output["xv"][20][:10], 8)
	outputXD2 := np.Round(output["xd"][20][:10], 8)

	if !reflect.DeepEqual(outputXA1, expectedXA1) {
		t.Errorf("xa(case1) Expected %v, got %v", expectedXA1, outputXA1)
	}
	if !reflect.DeepEqual(outputXV1, expectedXV1) {
		t.Errorf("xv(case1) Expected %v, got %v", expectedXV1, outputXV1)
	}
	if !reflect.DeepEqual(outputXD1, expectedXD1) {
		t.Errorf("xd(case1) Expected %v, got %v", expectedXD1, outputXD1)
	}
	if !reflect.DeepEqual(outputXA2, expectedXA2) {
		t.Errorf("xa(case2) Expected %v, got %v", expectedXA2, outputXA2)
	}
	if !reflect.DeepEqual(outputXV2, expectedXV2) {
		t.Errorf("xv(case2) Expected %v, got %v", expectedXV2, outputXV2)
	}
	if !reflect.DeepEqual(outputXD2, expectedXD2) {
		t.Errorf("xd(case2) Expected %v, got %v", expectedXD2, outputXD2)
	}

}

func TestCalcResponseSpectra(t *testing.T) {
	expectedSpectralAcc := np.Round(td.TestSpectra.SpectralAccelerations, 5)
	expectedSpectralVel := np.Round(td.TestSpectra.SpectralVelocities, 5)
	expectedSpectralDisp := np.Round(td.TestSpectra.SpectralDisplacements, 5)
	expectedPseudoAcc := np.Round(td.TestSpectra.PseudoAccelerations, 5)
	expectedPseudoVel := np.Round(td.TestSpectra.PseudoVelocities, 5)

	testAcceleration := td.TestMotion["Accelerations"].([]float64)
	output := ResponseSpectra(testAcceleration, 0.005, td.TestPeriods, 0.05)
	outputSpectralAcc := np.Round(output.SpectralAccelerations, 5)
	outputSpectralVel := np.Round(output.SpectralVelocities, 5)
	outputSpectralDisp := np.Round(output.SpectralDisplacements, 5)
	outputPseudoAcc := np.Round(output.PseudoAccelerations, 5)
	outputPseudoVel := np.Round(output.PseudoVelocities, 5)

	if !reflect.DeepEqual(expectedSpectralAcc, outputSpectralAcc) {
		t.Errorf("spectral_acceleration Expected %v, got %v", expectedSpectralAcc, outputSpectralAcc)
	}
	if !reflect.DeepEqual(expectedSpectralVel, outputSpectralVel) {
		t.Errorf("spectral_velocity Expected %v, got %v", expectedSpectralVel, outputSpectralVel)
	}
	if !reflect.DeepEqual(expectedSpectralDisp, outputSpectralDisp) {
		t.Errorf("spectral_displacement Expected %v, got %v", expectedSpectralDisp, outputSpectralDisp)
	}
	if !reflect.DeepEqual(expectedPseudoAcc, outputPseudoAcc) {
		t.Errorf("pseudo_acceleration Expected %v, got %v", expectedPseudoAcc, outputPseudoAcc)
	}
	if !reflect.DeepEqual(expectedPseudoVel, outputPseudoVel) {
		t.Errorf("pseudo_velocity Expected %v, got %v", expectedPseudoVel, outputPseudoVel)
	}
}
