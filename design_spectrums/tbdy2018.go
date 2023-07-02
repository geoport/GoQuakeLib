package design_spectrums

import (
	np "github.com/geoport/numpy4go/vectors"
	"math"
)

func InsertPeriods(periods []float64, startPeriod, endPeriod float64) []float64 {
	for _, period := range []float64{startPeriod, endPeriod} {
		if !np.Contains(periods, period) {
			periodIndex, _ := np.Where(
				periods, func(x float64) bool {
					return x < period
				},
			)
			periods = np.Insert(periods, periodIndex[len(periodIndex)-1]+1, period)
		}
	}

	return periods
}

func getPeriods(maxPeriod, periodStep float64) []float64 {
	periods := np.Arange(periodStep, maxPeriod, periodStep)
	periods[len(periods)-1] = maxPeriod
	return periods
}

// GetSpectrumByTBDY returns design spectrum and periods according to TBDY 2018
func GetSpectrumByTBDY(periodStep, maxPeriod, SDS, SD1 float64, insertPeriods bool) ([]float64, []float64) {
	periods := getPeriods(maxPeriod, periodStep)
	TA := 0.2 * SD1 / SDS
	TB := SD1 / SDS
	TL := 6.
	if insertPeriods {
		periods = InsertPeriods(periods, TA, TB)
	}
	spectrum := make([]float64, len(periods))
	for i, period := range periods {
		if period <= TA {
			spectrum[i] = (0.4 + 0.6*period/TA) * SDS
		} else if TA < period && period <= TB {
			spectrum[i] = SDS
		} else if TB < period && period <= TL {
			spectrum[i] = SD1 / period
		} else {
			spectrum[i] = SD1 * TL / math.Pow(period, 2)
		}
	}
	return spectrum, periods
}
