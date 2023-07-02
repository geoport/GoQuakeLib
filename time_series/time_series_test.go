package GoQuakeLib

import (
	td "github.com/geoport/GoQuakeLib/TestData"
	np "github.com/geoport/numpy4go/vectors"
	"reflect"
	"testing"
)

var testMotion = MotionData{
	Accelerations: td.TestMotion["Accelerations"].([]float64),
	Velocities:    td.TestMotion["Velocities"].([]float64),
	Displacements: td.TestMotion["Displacements"].([]float64),
	TimeStep:      td.TestMotion["TimeStep"].(float64),
	AccUnit:       td.TestMotion["AccUnit"].(string),
	VelUnit:       td.TestMotion["VelUnit"].(string),
	DispUnit:      td.TestMotion["DispUnit"].(string),
}

func TestMotionData_FromAcceleration(t *testing.T) {
	motion := testMotion
	motion.AccUnit = "m/s2"
	_, vel, disp := motion.FromAcceleration()
	if !reflect.DeepEqual(np.Round(vel, 4), np.Round(motion.Velocities, 4)) {
		t.Errorf("Velocities Expected %v, got %v", motion.Velocities[:10], vel[:10])
	}
	if !reflect.DeepEqual(np.Round(disp, 4), np.Round(motion.Displacements, 4)) {
		t.Errorf("Displacements Expected %v, got %v", motion.Displacements[:10], disp[:10])
	}
}

func TestMotionData_FromVelocity(t *testing.T) {
	motion := testMotion

	acc, _, disp := motion.FromVelocity()
	maxAcc := np.Max(acc)
	maxDisp := np.Max(disp)
	if np.Round(maxAcc, 3) != 0.015 {
		t.Errorf("PGA Expected %f, got %f", 0.015, maxAcc)
	}
	if np.Round(maxDisp, 3) != 1.035 {
		t.Errorf("PGD Expected %f, got %f", 1.035, maxDisp)
	}
}

func TestMotionData_FromDisplacement(t *testing.T) {
	motion := testMotion

	acc, vel, _ := motion.FromDisplacement()
	maxAcc := np.Max(acc)
	maxVel := np.Max(vel)
	if np.Round(maxAcc, 3) != 0.015 {
		t.Errorf("PGA Expected %f, got %f", 0.015, maxAcc)
	}
	if np.Round(maxVel, 3) != 2.423 {
		t.Errorf("PGV Expected %f, got %f", 2.423, maxVel)
	}
}
