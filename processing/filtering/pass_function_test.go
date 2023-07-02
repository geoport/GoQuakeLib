package Filtering

import (
	"github.com/geoport/numpy4go/vectors"
	"testing"
)

func TestLp2lpZpk(t *testing.T) {
	z, p, k := buttap(4)
	zLp, pLp, kLp := lp2lpZpk(z, p, k, []float64{2})
	if len(zLp) != 0 {
		t.Errorf("Expected no zeros, got %v", len(zLp))
	}
	if len(pLp) != 4 {
		t.Errorf("Expected 3 poles, got %v", len(pLp))
	}
	if kLp != 16 {
		t.Errorf("Expected gain of 1.0, got %v", kLp)
	}
}

func TestHp2lpZpk(t *testing.T) {
	z, p, k := buttap(4)
	zHp, pHp, kHp := lp2hpZpk(z, p, k, []float64{2})
	if len(zHp) != 4 {
		t.Errorf("Expected 4 zeros, got %v", len(zHp))
	}
	if len(pHp) != 4 {
		t.Errorf("Expected 4 poles, got %v", len(pHp))
	}
	if kHp != 1 {
		t.Errorf("Expected gain of 1.0, got %v", kHp)
	}
}

func TestBp2lpZpk(t *testing.T) {
	z, p, k := buttap(4)
	zBp, pBp, kBp := lp2bpZpk(z, p, k, 2, 2)
	if len(zBp) != 4 {
		t.Errorf("Expected 4 zeros, got %v", len(zBp))
	}
	if len(pBp) != 8 {
		t.Errorf("Expected 8 poles, got %v", len(pBp))
	}
	if kBp != 16 {
		t.Errorf("Expected gain of 16.0, got %v", kBp)
	}
}

func TestBs2lpZpk(t *testing.T) {
	z, p, k := buttap(4)
	zBs, pBs, kBs := lp2bsZpk(z, p, k, 2, 2)
	if len(zBs) != 8 {
		t.Errorf("Expected 8 zeros, got %v", len(zBs))
	}
	if len(pBs) != 8 {
		t.Errorf("Expected 8 poles, got %v", len(pBs))
	}
	if kBs != 1 {
		t.Errorf("Expected gain of 1, got %v", kBs)
	}
}

func TestBilinearZpk(t *testing.T) {
	z, p, k := buttap(4)
	zZ, pZ, kZ := bilinearZpk(z, p, k, 10)
	if len(zZ) != 4 {
		t.Errorf("Expected 4 zeros, got %v", len(zZ))
	}
	if len(pZ) != 4 {
		t.Errorf("Expected 4 poles, got %v", len(pZ))
	}
	if vectors.Round(kZ, 7) != 5.5e-6 {
		t.Errorf("Expected gain of 5.5e-6, got %v", vectors.Round(kZ, 7))
	}
}

func TestZpk2Tf(t *testing.T) {
	z, p, k := buttap(4)
	zZ, pZ, kZ := bilinearZpk(z, p, k, 10)
	b, a := zpk2Tf(zZ, pZ, kZ)

	if len(b) != 5 || vectors.Round(b[1], 6) != 2.2e-05 {
		t.Errorf("Unexpected b")
	}
	if len(a) != 5 || vectors.Round(a[2], 2) != 5.25 {
		t.Errorf("Unexpected a")
	}
}
