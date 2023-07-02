package Filtering

import (
	np "github.com/geoport/numpy4go/vectors"
	"math"
	"math/cmplx"
)

func buttap(N int) ([]complex128, []complex128, float64) {
	if N <= 0 {
		panic("Filter order must be a nonnegative integer")
	}

	var z []complex128
	m := np.Arange(-float64(N)+1, float64(N), 2)

	p := make([]complex128, N)
	for i, mi := range m {
		theta := complex(0., 1.) * complex(math.Pi, 0) * complex(mi, 0) / complex(2.*float64(N), 0)
		p[i] = -cmplx.Exp(theta)
	}

	k := 1.0

	return z, p, k
}

func cheb1ap(N int) ([]complex128, []complex128, float64) {
	if N <= 0 {
		panic("Filter order must be a nonnegative integer")
	}

	var z []complex128
	rp := 0.5

	eps := math.Sqrt(math.Pow(10, 0.1*rp) - 1)
	mu := math.Asinh(1/eps) / float64(N)

	m := np.Arange(-float64(N)+1, float64(N), 2)
	theta := np.MultiplyBy(m, math.Pi/float64(N*2))

	p := make([]complex128, len(m))
	pNeg := make([]complex128, len(m))

	for i, thetai := range theta {
		p[i] = -cmplx.Sinh(complex(mu, 0) + complex(thetai, 0)*complex(0, 1))
		pNeg[i] = p[i] * complex(-1, 0)
	}

	k := real(complexVectorProductSum(pNeg))

	if N%2 == 0 {
		k = k / math.Sqrt(1+eps*eps)
	}

	return z, p, k
}
