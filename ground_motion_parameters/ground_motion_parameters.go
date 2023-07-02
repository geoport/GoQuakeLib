package ground_motion_parameters

import (
	fs "github.com/geoport/GoQuakeLib/fourier_spectrum"
	Filtering "github.com/geoport/GoQuakeLib/processing/filtering"
	GoQuake "github.com/geoport/GoQuakeLib/response_spectra"
	ts "github.com/geoport/GoQuakeLib/time_series"
	"math"
	"sort"

	np "github.com/geoport/numpy4go/vectors"
)

type GMPData struct {
	Pga                           float64
	PgaTime                       float64
	Pgv                           float64
	PgvTime                       float64
	Pgd                           float64
	PgdTime                       float64
	HousnerIntensity              float64
	SustainedMaxAcceleration      float64
	SustainedMaxVelocity          float64
	EffectiveDesignAcceleration   float64
	AccelerationSpectrumIntensity float64
	VelocitySpectrumIntensity     float64
	A95                           float64
	PredominantPeriod             float64
	MeanPeriod                    float64
	UniformDuration               float64
	BracketedDuration             float64
	SignificantDuration           float64
	EffectiveDuration             float64
	AriasIntensity                float64
	AriasIntensityArray           []float64
	RmsAcceleration               float64
	RmsVelocity                   float64
	RmsDisplacement               float64
	CharacteristicIntensity       float64
	SpecificEnergyDensity         float64
	SpecificEnergyDensityArray    []float64
	CumulativeAbsoluteVelocity    float64
}

func (gmp *GMPData) CalcPGA(motion ts.MotionData) {
	pga := np.Max(np.Abs(motion.Accelerations))
	pgaIndex := np.ArgMax(np.Abs(motion.Accelerations))
	gmp.Pga = pga
	gmp.PgaTime = motion.Times[pgaIndex]
}

func (gmp *GMPData) CalcPGV(motion ts.MotionData) {
	pgv := np.Max(np.Abs(motion.Velocities))
	pgvIndex := np.ArgMax(np.Abs(motion.Velocities))
	gmp.Pgv = pgv
	gmp.PgvTime = motion.Times[pgvIndex]
}

func (gmp *GMPData) CalcPGD(motion ts.MotionData) {
	pgd := np.Max(np.Abs(motion.Displacements))
	pgdIndex := np.ArgMax(np.Abs(motion.Displacements))
	gmp.Pgd = pgd
	gmp.PgdTime = motion.Times[pgdIndex]
}

func (gmp *GMPData) CalcHousnerIntensity(spectraData *GoQuake.ResponseSpectraData) {
	periods := spectraData.Periods
	dx := periods[2] - periods[1]
	filterFunc := func(x float64) bool { return x >= 0.1 && x <= 2.5 }
	indexes, _ := np.Where(np.Round(periods, 2).([]float64), filterFunc)
	psv := spectraData.PseudoVelocities[indexes[0]:indexes[len(indexes)-1]]
	gmp.HousnerIntensity = np.Cumtrapz(psv, dx, 0)[len(psv)-1]
}

func (gmp *GMPData) CalcSustainedMaxAcceleration(motion ts.MotionData) {
	accelerations := np.Abs(motion.Accelerations)
	sort.Float64s(accelerations)
	gmp.SustainedMaxAcceleration = accelerations[len(accelerations)-3]
}

func (gmp *GMPData) CalcSustainedMaxVelocity(motion ts.MotionData) {
	velocities := np.Abs(motion.Velocities)
	sort.Float64s(velocities)
	gmp.SustainedMaxVelocity = velocities[len(velocities)-3]
}

func (gmp *GMPData) CalcEffectiveDesignAcceleration(motion ts.MotionData) {
	accelerations := np.Abs(motion.Accelerations)
	dt := motion.TimeStep
	cornerFreqs := []float64{9}
	err, filteredAccelerations := Filtering.FilterSignal(accelerations, cornerFreqs, 1, "lowpass", "butterworth", dt)
	if err != nil {
		panic(err)
	}
	gmp.EffectiveDesignAcceleration = np.Max(np.Abs(filteredAccelerations))
}

func (gmp *GMPData) CalcAccelerationSpectrumIntensity(spectraData *GoQuake.ResponseSpectraData) {
	periods := spectraData.Periods
	dx := periods[2] - periods[1]
	filterFunc := func(x float64) bool { return x >= 0.1 && x <= 0.5 }
	indexes, _ := np.Where(np.Round(periods, 2).([]float64), filterFunc)
	sa := spectraData.SpectralAccelerations[indexes[0]:indexes[len(indexes)-1]]
	gmp.AccelerationSpectrumIntensity = np.Cumtrapz(sa, dx, 0)[len(sa)-1]
}

func (gmp *GMPData) CalcVelocitySpectrumIntensity(spectraData *GoQuake.ResponseSpectraData) {
	periods := spectraData.Periods
	dx := periods[2] - periods[1]
	filterFunc := func(x float64) bool { return x >= 0.1 && x <= 2.5 }
	indexes, _ := np.Where(np.Round(periods, 2).([]float64), filterFunc)
	sv := spectraData.SpectralVelocities[indexes[0]:indexes[len(indexes)-1]]
	gmp.VelocitySpectrumIntensity = np.Cumtrapz(sv, dx, 0)[len(sv)-1]
}

func (gmp *GMPData) CalcA95(motion ts.MotionData) {
	accelerations := motion.Accelerations
	dt := motion.TimeStep
	// Es is arias intensity
	Es := np.Cumtrapz(np.Pow(accelerations, 2), dt, 0)[len(accelerations)-1]
	callback := func(A2 float64) float64 {
		clippedA2 := np.Clip(np.Pow(accelerations, 2), 0, A2)
		// Ex is the area between the line at A95^2 and pga^2
		Ex := Es - np.Cumtrapz(clippedA2, dt, 0)[len(accelerations)-1]

		return Ex/Es - 0.05
	}
	optimizer := func() float64 {
		boundary1 := 0.0
		boundary2 := np.Max(np.Pow(accelerations, 2))
		A2 := (boundary1 + boundary2) / 2

		n := 0

		// keep iterating till the difference is less than 0.001 or number of steps is greater than 100
		for math.Abs(callback(A2)) > 0.01 || n < 100 {
			n++
			if boundary1 == boundary2 && boundary1 == A2 && n > 10 {
				return 0
			}
			if callback(A2) > 0 {
				boundary1 = A2
			} else {
				boundary2 = A2
			}
			A2 = (boundary1 + boundary2) / 2
		}
		return math.Sqrt(A2)
	}
	A95 := optimizer()
	gmp.A95 = A95
}

func (gmp *GMPData) CalcPredominantPeriod(spectraData *GoQuake.ResponseSpectraData) {
	sa := spectraData.SpectralAccelerations
	maxSAIndex := np.ArgMax(sa)

	gmp.PredominantPeriod = spectraData.Periods[maxSAIndex]
}

func (gmp *GMPData) CalcMeanPeriod(motion ts.MotionData) {
	f, fa, _ := fs.FourierSpectrum(motion.Accelerations, motion.TimeStep)
	indexes, _ := np.Where(f, func(x float64) bool { return x >= 0.25 && x <= 20 })

	var A, B float64
	for _, i := range indexes {
		A += math.Pow(fa[i], 2) / f[i]
		B += math.Pow(fa[i], 2)
	}
	gmp.MeanPeriod = A / B
}

func (gmp *GMPData) CalcUniformDuration(motion ts.MotionData) {
	accelerations := motion.Accelerations
	A0 := np.Max(np.Abs(motion.Accelerations)) * 0.05
	higherIndexes, _ := np.Where(np.Pow(accelerations, 2), func(x float64) bool { return x > A0*A0 })
	gmp.UniformDuration = motion.TimeStep * float64(len(higherIndexes))
}

func (gmp *GMPData) CalcBracketedDuration(motion ts.MotionData) {
	accelerations := motion.Accelerations
	A0 := np.Max(np.Abs(motion.Accelerations)) * 0.05
	indexes, _ := np.Where(np.Pow(accelerations, 2), func(x float64) bool { return x >= A0*A0 })
	gmp.BracketedDuration = motion.TimeStep + motion.Times[indexes[len(indexes)-1]] - motion.Times[indexes[0]]
}

func (gmp *GMPData) CalcAriasIntensity(motion ts.MotionData) {
	g := 9.81 // m/s^2
	accelerations := motion.Accelerations
	dt := motion.TimeStep
	acc2 := np.Pow(np.MultiplyBy(accelerations, g), 2)
	Ia := np.MultiplyBy(np.Cumtrapz(acc2, dt, 0), math.Pi*0.5/g)
	gmp.AriasIntensity = Ia[len(Ia)-1]
	gmp.AriasIntensityArray = Ia
}

func (gmp *GMPData) CalcSignificantDuration(motion ts.MotionData) {
	p1 := 5.
	p2 := 95.
	if gmp.AriasIntensity == 0 {
		gmp.CalcAriasIntensity(motion)
	}
	Ia := gmp.AriasIntensityArray
	IaNormalized := np.MultiplyBy(Ia, 100/np.Max(Ia))

	indexes, _ := np.Where(np.Round(IaNormalized, 2).([]float64), func(x float64) bool { return x >= p1 && x <= p2 })
	gmp.SignificantDuration = motion.Times[indexes[len(indexes)-1]] - motion.Times[indexes[0]]
}

func (gmp *GMPData) CalcRMSAcceleration(motion ts.MotionData) {
	accelerations := motion.Accelerations
	acc2 := np.Pow(accelerations, 2)
	acc2Cum := np.Cumtrapz(acc2, motion.TimeStep, 0)
	gmp.RmsAcceleration = math.Sqrt(acc2Cum[len(acc2Cum)-1] / motion.Times[len(motion.Times)-1])
}

func (gmp *GMPData) CalcRMSVelocity(motion ts.MotionData) {
	velocities := motion.Velocities
	vel2 := np.Pow(velocities, 2)
	vel2Cum := np.Cumtrapz(vel2, motion.TimeStep, 0)
	gmp.RmsVelocity = math.Sqrt(vel2Cum[len(vel2Cum)-1] / motion.Times[len(motion.Times)-1])
}

func (gmp *GMPData) CalcRMSDisplacement(motion ts.MotionData) {
	displacements := motion.Displacements
	disp2 := np.Pow(displacements, 2)
	disp2Cum := np.Cumtrapz(disp2, motion.TimeStep, 0)
	gmp.RmsDisplacement = math.Sqrt(disp2Cum[len(disp2Cum)-1] / motion.Times[len(motion.Times)-1])
}

func (gmp *GMPData) CalcCharacteristicIntensity(motion ts.MotionData) {
	if gmp.RmsAcceleration == 0 {
		gmp.CalcRMSAcceleration(motion)
	}
	RMSAccel := gmp.RmsAcceleration
	gmp.CharacteristicIntensity = math.Pow(RMSAccel, 1.5) * math.Sqrt(motion.Times[len(motion.Times)-1])
}

func (gmp *GMPData) CalcSpecificEnergyDensity(motion ts.MotionData) {
	velocities := motion.Velocities
	sedArray := np.Cumtrapz(np.Pow(velocities, 2), motion.TimeStep, 0)
	gmp.SpecificEnergyDensity = sedArray[len(sedArray)-1]
	gmp.SpecificEnergyDensityArray = sedArray
}

func (gmp *GMPData) CalcCumulativeAbsoluteVelocity(motion ts.MotionData) {
	g := 9.81 // m/s^2
	accelerations := motion.Accelerations
	CAV := np.Cumtrapz(np.Abs(accelerations), motion.TimeStep, 0)[len(accelerations)-1] * 100 * g
	gmp.CumulativeAbsoluteVelocity = CAV
}

func (gmp *GMPData) CalcGMP(motion ts.MotionData, spectra *GoQuake.ResponseSpectraData) *GMPData {
	gmp.CalcAriasIntensity(motion)
	gmp.CalcPGA(motion)
	gmp.CalcPGV(motion)
	gmp.CalcPGD(motion)
	gmp.CalcHousnerIntensity(spectra)
	gmp.CalcSustainedMaxAcceleration(motion)
	gmp.CalcSustainedMaxVelocity(motion)
	gmp.CalcEffectiveDesignAcceleration(motion)
	gmp.CalcAccelerationSpectrumIntensity(spectra)
	gmp.CalcVelocitySpectrumIntensity(spectra)
	gmp.CalcA95(motion)
	gmp.CalcPredominantPeriod(spectra)
	gmp.CalcMeanPeriod(motion)
	gmp.CalcUniformDuration(motion)
	gmp.CalcBracketedDuration(motion)
	gmp.CalcSignificantDuration(motion)
	gmp.CalcRMSAcceleration(motion)
	gmp.CalcRMSVelocity(motion)
	gmp.CalcRMSDisplacement(motion)
	gmp.CalcCharacteristicIntensity(motion)
	gmp.CalcSpecificEnergyDensity(motion)
	gmp.CalcCumulativeAbsoluteVelocity(motion)

	return gmp
}
