package ground_motion_parameters

import (
	td "github.com/geoport/GoQuakeLib/TestData"
	rs "github.com/geoport/GoQuakeLib/response_spectra"
	ts "github.com/geoport/GoQuakeLib/time_series"
	"testing"

	np "github.com/geoport/numpy4go/vectors"
)

var testMotion = ts.MotionData{
	Accelerations: td.TestMotion["Accelerations"].([]float64),
	TimeStep:      td.TestMotion["TimeStep"].(float64),
	AccUnit:       td.TestMotion["AccUnit"].(string),
	VelUnit:       td.TestMotion["VelUnit"].(string),
	DispUnit:      td.TestMotion["DispUnit"].(string),
	Times:         td.TestMotion["Times"].([]float64),
}

var testPeriods = np.Arange(0, 4, 0.02)
var testSpectra = rs.ResponseSpectra(testMotion.Accelerations, testMotion.TimeStep, testPeriods, 0.05)

var gmp = GMPData{}

func TestGMP_CalcPGA(t *testing.T) {
	expectedPGA := 0.16076
	expectedPGATime := 13.35
	gmp.CalcPGA(testMotion)
	outputPGA := np.Round(gmp.Pga, 5)
	outputPGATime := np.Round(gmp.PgaTime, 2)
	if outputPGA != expectedPGA {
		t.Errorf("Expected PGA = %f, got %f", expectedPGA, outputPGA)
	}
	if outputPGATime != expectedPGATime {
		t.Errorf("Expected PGA Time = %f, got %f", expectedPGATime, outputPGATime)
	}
}

func TestGMP_CalcPGV(t *testing.T) {
	testMotion.FromAcceleration()
	expectedPGV := 29.37936
	expectedPGVTime := 9.7
	gmp.CalcPGV(testMotion)
	outputPGV := np.Round(gmp.Pgv, 5)
	outputPGVTime := np.Round(gmp.PgvTime, 2)
	if outputPGV != expectedPGV {
		t.Errorf("Expected PGV = %f, got %f", expectedPGV, outputPGV)
	}
	if outputPGVTime != expectedPGVTime {
		t.Errorf("Expected PGV Time = %f, got %f", expectedPGVTime, outputPGVTime)
	}
}

func TestGMP_CalcPGD(t *testing.T) {
	testMotion.FromAcceleration()
	expectedPGD := 33.87
	expectedPGDTime := 12.
	gmp.CalcPGD(testMotion)
	outputPGD := np.Round(gmp.Pgd, 2)
	outputPGDTime := np.Round(gmp.PgdTime, 2)
	if outputPGD != expectedPGD {
		t.Errorf("Expected PGD = %f, got %f", expectedPGD, outputPGD)
	}
	if outputPGDTime != expectedPGDTime {
		t.Errorf("Expected PGD Time = %f, got %f", expectedPGDTime, outputPGDTime)
	}
}

func TestGMP_CalcHousnerIntensity(t *testing.T) {
	expected := 104.42
	gmp.CalcHousnerIntensity(testSpectra)
	output := np.Round(gmp.HousnerIntensity, 2)
	if output != expected {
		t.Errorf("Expected HI = %f, got %f", expected, output)
	}
}

func TestGMP_CalcSustainedMaxAcceleration(t *testing.T) {
	expected := 0.13839
	gmp.CalcSustainedMaxAcceleration(testMotion)
	output := np.Round(gmp.SustainedMaxAcceleration, 5)
	if output != expected {
		t.Errorf("Expected SMA = %f, got %f", expected, output)
	}
}

func TestGMP_CalcSustainedMaxVelocity(t *testing.T) {
	testMotion.FromAcceleration()
	expected := 28.30913
	gmp.CalcSustainedMaxVelocity(testMotion)
	output := np.Round(gmp.SustainedMaxVelocity, 5)
	if output != expected {
		t.Errorf("Expected SMA = %f, got %f", expected, output)
	}
}

func TestGMP_CalcEffectiveDesignAcceleration(t *testing.T) {
	expected := 0.16586
	gmp.CalcEffectiveDesignAcceleration(testMotion)
	output := np.Round(gmp.EffectiveDesignAcceleration, 5)
	if output != expected {
		t.Errorf("Expected EDA = %f, got %f", expected, output)
	}
}

func TestGMP_CalcAccelerationSpectrumIntensity(t *testing.T) {
	expected := 0.11127
	gmp.CalcAccelerationSpectrumIntensity(testSpectra)
	output := np.Round(gmp.AccelerationSpectrumIntensity, 5)
	if output != expected {
		t.Errorf("Expected ASI = %f, got %f", expected, output)
	}
}

func TestGMP_CalcVelocitySpectrumIntensity(t *testing.T) {
	expected := 105.78267
	gmp.CalcVelocitySpectrumIntensity(testSpectra)
	output := np.Round(gmp.VelocitySpectrumIntensity, 5)
	if output != expected {
		t.Errorf("Expected ASI = %f, got %f", expected, output)
	}
}

func TestGMP_CalcA95(t *testing.T) {
	expected := 0.12115
	gmp.CalcA95(testMotion)
	output := np.Round(gmp.A95, 5)
	if output != expected {
		t.Errorf("Expected A95 = %f, got %f", expected, output)
	}
}

func TestGMP_CalcPredominantPeriod(t *testing.T) {
	expected := 0.7
	gmp.CalcPredominantPeriod(testSpectra)
	output := np.Round(gmp.PredominantPeriod, 5)
	if output != expected {
		t.Errorf("Expected Predominant Period = %f, got %f", expected, output)
	}
}

func TestGMP_CalcMeanPeriod(t *testing.T) {
	expected := 1.05355
	gmp.CalcMeanPeriod(testMotion)
	output := np.Round(gmp.MeanPeriod, 5)
	if output != expected {
		t.Errorf("Expected Mean Period = %f, got %f", expected, output)
	}
}

func TestGMP_CalcUniformDuration(t *testing.T) {
	expected := 14.2
	gmp.CalcUniformDuration(testMotion)
	output := np.Round(gmp.UniformDuration, 5)
	if output != expected {
		t.Errorf("Expected Uniform Duration = %f, got %f", expected, output)
	}
}

func TestGMP_CalcBracketedDuration(t *testing.T) {
	expected := 36.3
	gmp.CalcBracketedDuration(testMotion)
	output := np.Round(gmp.BracketedDuration, 5)
	if output != expected {
		t.Errorf("Expected Bracketed Duration = %f, got %f", expected, output)
	}
}

func TestGMP_CalcAriasIntensity(t *testing.T) {
	expected := 0.34798
	gmp.CalcAriasIntensity(testMotion)
	output := np.Round(gmp.AriasIntensity, 5)
	if output != expected {
		t.Errorf("Expected Arias Intensity = %f, got %f", expected, output)
	}
}

func TestGMP_CalcSignificantDuration(t *testing.T) {
	expected := 12.05
	gmp.CalcSignificantDuration(testMotion)
	output := np.Round(gmp.SignificantDuration, 5)
	if output != expected {
		t.Errorf("Expected Significant Duration = %f, got %f", expected, output)
	}
}

func TestGMP_CalcRMSAcceleration(t *testing.T) {
	expected := 0.00942
	gmp.CalcRMSAcceleration(testMotion)
	output := np.Round(gmp.RmsAcceleration, 5)
	if output != expected {
		t.Errorf("Expected RMS Acceleration = %f, got %f", expected, output)
	}
}

func TestGMP_CalcRMSVelocity(t *testing.T) {
	expected := 2.535
	testMotion.FromAcceleration()
	gmp.CalcRMSVelocity(testMotion)
	output := np.Round(gmp.RmsVelocity, 3)
	if output != expected {
		t.Errorf("Expected RMS Velocity = %f, got %f", expected, output)
	}
}

func TestGMP_CalcRMSDisplacement(t *testing.T) {
	expected := 4.713
	testMotion.FromAcceleration()
	gmp.CalcRMSDisplacement(testMotion)
	output := np.Round(gmp.RmsDisplacement, 3)
	if output != expected {
		t.Errorf("Expected RMS Displacement = %f, got %f", expected, output)
	}
}

func TestGMP_CalcCharacteristicIntensity(t *testing.T) {
	expected := 0.01458
	gmp.CalcCharacteristicIntensity(testMotion)
	output := np.Round(gmp.CharacteristicIntensity, 5)
	if output != expected {
		t.Errorf("Expected Characteristic Intensity = %f, got %f", expected, output)
	}
}

func TestGMP_CalcSpecificEnergyDensity(t *testing.T) {
	expected := 1636.71118
	testMotion.FromAcceleration()
	gmp.CalcSpecificEnergyDensity(testMotion)
	output := np.Round(gmp.SpecificEnergyDensity, 5)
	if output != expected {
		t.Errorf("Expected Specific Energy Density = %f, got %f", expected, output)
	}
}

func TestGMP_CalcCumulativeAbsoluteVelocity(t *testing.T) {
	expected := 626.62117
	gmp.CalcCumulativeAbsoluteVelocity(testMotion)
	output := np.Round(gmp.CumulativeAbsoluteVelocity, 5)
	if output != expected {
		t.Errorf("Expected Cumulative Absolute Velocity = %f, got %f", expected, output)
	}
}
