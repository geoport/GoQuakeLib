package GoQuakeLib

import (
	np "github.com/geoport/numpy4go/vectors"
)

type MotionData struct {
	Accelerations []float64
	Velocities    []float64
	Displacements []float64
	Times         []float64
	TimeStep      float64
	AccUnit       string
	VelUnit       string
	DispUnit      string
}

func (md *MotionData) FromAcceleration() ([]float64, []float64, []float64) {
	unitDict := map[string]float64{
		"g":      1,
		"m/s2":   1 / 9.81,
		"cm/s2":  1 / 981,
		"mm/s2":  1 / 9810,
		"inc/s2": 0.0025900792,
		"ft/s2":  1 / 32.17404855643,
	}
	if unitDict[md.AccUnit] == 0 {
		panic("Unsupported unit")
	}
	if len(md.Accelerations) == 0 {
		panic("No acceleration data")
	}
	if md.TimeStep == 0 {
		panic("Time step is zero")
	}
	md.Accelerations = np.MultiplyBy(md.Accelerations, unitDict[md.AccUnit]) //g

	velocities := np.Cumtrapz(md.Accelerations, md.TimeStep, md.Accelerations[0])
	md.Velocities = np.MultiplyBy(velocities, 981) // cm/s

	displacements := np.Cumtrapz(md.Velocities, md.TimeStep, md.Velocities[0])
	md.Displacements = displacements // cm

	md.Times = np.Arange(0, float64(len(md.Accelerations))*md.TimeStep, md.TimeStep)

	return md.Accelerations, md.Velocities, md.Displacements
}

func (md *MotionData) FromVelocity() ([]float64, []float64, []float64) {
	unitDict := map[string]float64{
		"m/s":   100,
		"cm/s":  1,
		"mm/s":  0.1,
		"inc/s": 2.54,
		"ft/s":  30.48,
	}
	if unitDict[md.VelUnit] == 0 {
		panic("Unsupported unit")
	}
	if len(md.Velocities) == 0 {
		panic("No velocity data")
	}
	if md.TimeStep == 0 {
		panic("Time step is zero")
	}

	md.Times = np.Arange(0, float64(len(md.Velocities))*md.TimeStep, md.TimeStep)
	md.Velocities = np.MultiplyBy(md.Velocities, unitDict[md.VelUnit]) //cm/s
	accelerations := np.DividedBy(np.Diff(md.Velocities), np.Diff(md.Times))
	accelerations = np.MultiplyBy(np.Insert(accelerations, 0, md.Velocities[0]), 1.0/981)
	md.Accelerations = accelerations // g

	displacements := np.Cumtrapz(md.Velocities, md.TimeStep, md.Velocities[0])
	md.Displacements = displacements // cm

	return md.Accelerations, md.Velocities, md.Displacements
}

func (md *MotionData) FromDisplacement() ([]float64, []float64, []float64) {
	unitDict := map[string]float64{
		"m":   100,
		"cm":  1,
		"mm":  0.1,
		"inc": 2.54,
		"ft":  30.48,
	}
	if unitDict[md.DispUnit] == 0 {
		panic("Unsupported unit")
	}
	if len(md.Displacements) == 0 {
		panic("No displacement data")
	}
	if md.TimeStep == 0 {
		panic("Time step is zero")
	}

	md.Times = np.Arange(0, float64(len(md.Displacements))*md.TimeStep, md.TimeStep)
	md.Displacements = np.MultiplyBy(md.Displacements, unitDict[md.DispUnit]) // cm

	velocities := np.DividedBy(np.Diff(md.Displacements), np.Diff(md.Times))
	velocities = np.Insert(velocities, 0, md.Displacements[0])
	md.Velocities = velocities // cm/s

	accelerations := np.DividedBy(np.Diff(md.Velocities), np.Diff(md.Times))
	accelerations = np.MultiplyBy(np.Insert(accelerations, 0, md.Velocities[0]), 1.0/981) // g
	md.Accelerations = accelerations

	return md.Accelerations, md.Velocities, md.Displacements
}
