package Filtering

import (
	np "github.com/geoport/numpy4go/vectors"
	"testing"
)

func TestButtap(t *testing.T) {
	z, p, k := buttap(4)
	if len(z) != 0 {
		t.Errorf("Expected no zeros, got %v", z)
	}
	if len(p) != 4 {
		t.Errorf("Expected 4 poles, got %v", len(p))
	}
	if k != 1.0 {
		t.Errorf("Expected gain of 1.0, got %v", k)
	}
}

func TestCheb1ap(t *testing.T) {
	z, p, k := cheb1ap(2)
	if len(z) != 0 {
		t.Errorf("Expected no zeros, got %v", z)
	}
	if len(p) != 2 {
		t.Errorf("Expected 2 poles, got %v", len(p))
	}
	if np.Round(k, 2) != 1.43 {
		t.Errorf("Expected gain of 1.43, got %v", k)
	}

}
