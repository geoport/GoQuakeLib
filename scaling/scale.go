package Scaling

import (
	"errors"
	np "github.com/geoport/numpy4go/vectors"
)

func isPeriodInRange(periods []float64, period float64) bool {
	return np.Contains(np.Round(periods, 3).([]float64), np.Round(period, 3).(float64))
}

// ScaleByMSE returns the scale factor for the response spectrum to match the target spectrum based on minimum mean squared error.
func ScaleByMSE(
	periods, targetSpectrum []float64, spectralAccelerations [][]float64, minScalingPeriod,
	maxScalingPeriod float64,
) ([]float64, error) {
	var err error

	scaleFactors := make([]float64, len(spectralAccelerations))

	if !isPeriodInRange(periods, minScalingPeriod) || !isPeriodInRange(periods, maxScalingPeriod) {
		err = errors.New("periods must be in the list of periods")
		return scaleFactors, err
	}

	indexT1 := np.ArgMin(np.Abs(np.SumWith(periods, -minScalingPeriod)))
	indexT2 := np.ArgMin(np.Abs(np.SumWith(periods, -maxScalingPeriod)))

	for i, responseSpectrum := range spectralAccelerations {
		responseRange := responseSpectrum[indexT1 : indexT2+1]
		targetRange := targetSpectrum[indexT1 : indexT2+1]
		scaleFactors[i] = np.Sum(np.MultiplyBy(targetRange, responseRange)) / np.Sum(np.Pow(responseRange, 2))
	}

	return scaleFactors, err
}

// ScaleByPeriodRange returns the scale factor for the response spectrum to match the target spectrum
// so that all the spectral values between specified periods will be above target spectrum.
func ScaleByPeriodRange(
	periods, targetSpectrum []float64, spectralAccelerations [][]float64, minScalingPeriod,
	maxScalingPeriod float64,
) ([]float64, error) {
	var err error
	scaleFactors := make([]float64, len(spectralAccelerations))

	if !isPeriodInRange(periods, minScalingPeriod) || !isPeriodInRange(periods, maxScalingPeriod) {
		err = errors.New("periods must be in the list of periods")
		return scaleFactors, err
	}

	indexT1 := np.ArgMin(np.Abs(np.SumWith(periods, -minScalingPeriod)))
	indexT2 := np.ArgMin(np.Abs(np.SumWith(periods, -maxScalingPeriod)))

	for i, responseSpectrum := range spectralAccelerations {
		responseRange := responseSpectrum[indexT1 : indexT2+1]
		targetRange := targetSpectrum[indexT1 : indexT2+1]
		scaleFactors[i] = np.Max(np.DividedBy(targetRange, responseRange))
	}

	return scaleFactors, err
}

// ScaleBySinglePeriod returns the scale factor for the response spectrum to match the target spectrum
// so that spectral value at specified period scaled to the target spectrum.
func ScaleBySinglePeriod(
	periods, targetSpectrum []float64, spectralAccelerations [][]float64,
	scalingPeriod float64,
) ([]float64, error) {
	var err error
	scaleFactors := make([]float64, len(spectralAccelerations))

	if !isPeriodInRange(periods, scalingPeriod) {
		err = errors.New("scaling period must be in the list of periods")
		return scaleFactors, err
	}

	index := np.ArgMin(np.Abs(np.SumWith(periods, -scalingPeriod)))

	for i, responseSpectrum := range spectralAccelerations {
		scaleFactors[i] = targetSpectrum[index] / responseSpectrum[index]
	}

	return scaleFactors, err
}
