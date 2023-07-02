package response_spectra

import (
	np "github.com/geoport/numpy4go/vectors"
	"math"
)

type ResponseSpectraData struct {
	SpectralAccelerations []float64
	SpectralVelocities    []float64
	SpectralDisplacements []float64
	PseudoAccelerations   []float64
	PseudoVelocities      []float64
	Periods               []float64
}

func GetTimeSeries(
	constants map[string][]float64, omega2 []float64, numSteps, numPeriods int, accelerations []float64, dt float64,
) map[string][][]float64 {
	xd := np.Zeros(numSteps-1, numPeriods)
	xv := np.Zeros(numSteps-1, numPeriods)
	xa := np.Zeros(numSteps-1, numPeriods)

	for i := 0; i < numSteps-1; i++ {
		dug := accelerations[i+1] - accelerations[i]
		z1 := np.MultiplyBy(constants["f2"], dug)
		z2 := np.MultiplyBy(constants["f2"], accelerations[i])
		z3 := np.MultiplyBy(constants["f1"], dug)
		z4 := np.MultiplyBy(z1, 1/dt)

		var a, b []float64
		if i == 0 {
			b = np.SumWith(z2, np.MultiplyBy(z3, -1))
			a = np.SumWith(np.MultiplyBy(constants["f5"], b), np.MultiplyBy(constants["f4"], z4))
		} else {
			b = np.SumWith(xd[i-1], np.SumWith(z2, np.MultiplyBy(z3, -1)))
			a = np.SumWith(
				np.SumWith(
					np.MultiplyBy(constants["f5"], b), np.MultiplyBy(constants["f4"], z4),
				),
				np.MultiplyBy(constants["f4"], xv[i-1]),
			)
		}

		z321 := np.SumWith(z3, np.MultiplyBy(np.SumWith(z1, z2), -1))
		ag1 := np.MultiplyBy(a, constants["g1"])
		ah1 := np.MultiplyBy(a, constants["h1"])
		bg2 := np.MultiplyBy(b, constants["g2"])
		bh2 := np.MultiplyBy(b, constants["h2"])

		xd[i] = np.SumWith(ag1, np.SumWith(bg2, z321))
		xv[i] = np.SumWith(ah1, np.MultiplyBy(np.SumWith(bh2, z4), -1))
		o2xd := np.MultiplyBy(omega2, xd[i])
		f6xv := np.MultiplyBy(constants["f6"], xv[i])
		xa[i] = np.MultiplyBy(np.SumWith(f6xv, o2xd), -1)
	}
	timeSeriesData := map[string][][]float64{
		"xd": xd,
		"xv": xv,
		"xa": xa,
	}
	return timeSeriesData
}

func ResponseSpectra(accelerations []float64, dt float64, periods []float64, damping float64) *ResponseSpectraData {
	if periods[0] == 0 {
		periods[0] = 1e-6
	}
	omega := np.DividedBy(np.Repeat(2*math.Pi, len(periods)), periods)
	omega2 := np.Pow(omega, 2)
	omega3 := np.Pow(omega, 3)
	omegaD := np.MultiplyBy(omega, math.Sqrt(1-damping*damping))
	var constants = map[string][]float64{}
	constants["f1"] = np.DividedBy(np.Repeat(2*damping/dt, len(omega)), omega3)
	constants["f2"] = np.DividedBy(np.Ones(len(omega)), omega2)
	constants["f3"] = np.MultiplyBy(omega, damping)
	constants["f4"] = np.DividedBy(np.Ones(len(omega)), omegaD)
	constants["f5"] = np.MultiplyBy(constants["f3"], constants["f4"])
	constants["f6"] = np.MultiplyBy(constants["f3"], 2)
	constants["e"] = np.Apply(np.MultiplyBy(constants["f3"], -dt), math.Exp)
	constants["s"] = np.Apply(np.MultiplyBy(omegaD, dt), math.Sin)
	constants["c"] = np.Apply(np.MultiplyBy(omegaD, dt), math.Cos)
	constants["g1"] = np.MultiplyBy(constants["e"], constants["s"])
	constants["g2"] = np.MultiplyBy(constants["e"], constants["c"])

	oDg1 := np.MultiplyBy(omegaD, constants["g1"])
	oDg2 := np.MultiplyBy(omegaD, constants["g2"])
	f3g1 := np.MultiplyBy(constants["f3"], constants["g1"])
	f3g2 := np.MultiplyBy(constants["f3"], constants["g2"])
	constants["h1"] = np.SumWith(oDg2, np.MultiplyBy(f3g1, -1))
	constants["h2"] = np.SumWith(oDg1, f3g2)

	timeSeriesData := GetTimeSeries(constants, omega2, len(accelerations), len(periods), accelerations, dt)
	var spectraData = ResponseSpectraData{
		Periods:               periods,
		SpectralAccelerations: np.Max2D(np.Abs2D(timeSeriesData["xa"]), 0),                     // g
		SpectralVelocities:    np.MultiplyBy(np.Max2D(np.Abs2D(timeSeriesData["xv"]), 0), 981), // cm/s
		SpectralDisplacements: np.MultiplyBy(np.Max2D(np.Abs2D(timeSeriesData["xd"]), 0), 981), // cm
	}
	spectraData.PseudoVelocities = np.MultiplyBy(omega, spectraData.SpectralDisplacements) // cm/s
	spectraData.PseudoAccelerations = np.MultiplyBy(
		np.MultiplyBy(omega2, spectraData.SpectralDisplacements), 1.0/981,
	) // g

	return &spectraData
}
