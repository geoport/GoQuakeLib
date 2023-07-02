package Filtering

import (
	"github.com/geoport/numpy4go/vectors"
	"testing"
)

func TestPoly(t *testing.T) {
	roots := []complex128{
		complex(0.958, 0.088),
		complex(0.911, 0.035),
		complex(0.911, -0.035),
		complex(0.958, -0.088),
	}
	coefficients := Poly(roots)

	if len(coefficients) != 5 || vectors.Round(coefficients[2], 2) != 5.25 {
		t.Errorf("Expected coefficients of [ 1.,-3.73,5.25,-3.28,0.77], got %v", coefficients)
	}
}
