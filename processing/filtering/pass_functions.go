package Filtering

import (
	"github.com/geoport/numpy4go/vectors"
	"math"
)

func lp2lpZpk(z []complex128, p []complex128, k float64, wo []float64) ([]complex128, []complex128, float64) {
	z = atleast1d(z).([]complex128)
	p = atleast1d(p).([]complex128)
	degree := relativeDegree(z, p)

	zLp := make([]complex128, len(z))
	pLp := make([]complex128, len(p))
	kLp := k * math.Pow(wo[0], float64(degree))

	for i, zi := range z {
		zLp[i] = zi * complex(wo[0], 0)
	}

	for i, pi := range p {
		pLp[i] = pi * complex(wo[0], 0)
	}

	return zLp, pLp, kLp
}

func lp2hpZpk(z []complex128, p []complex128, k float64, wo []float64) ([]complex128, []complex128, float64) {
	z = atleast1d(z).([]complex128)
	p = atleast1d(p).([]complex128)
	degree := relativeDegree(z, p)

	zHp := make([]complex128, len(z))
	zNeg := make([]complex128, len(z))
	pHp := make([]complex128, len(p))
	pNeg := make([]complex128, len(p))

	for i, zi := range z {
		zHp[i] = complex(wo[0], 0) / zi
		zNeg[i] = complex(-1, 0) * zi
	}
	for i := 0; i < degree; i++ {
		zHp = append(zHp, 0)
	}

	for i, pi := range p {
		pHp[i] = complex(wo[0], 0) / pi
		pNeg[i] = complex(-1, 0) * pi
	}
	kHp := k * real(complexVectorProductSum(zNeg)/complexVectorProductSum(pNeg))

	return zHp, pHp, kHp
}

func lp2bpZpk(z []complex128, p []complex128, k float64, wo float64, bw float64) ([]complex128, []complex128, float64) {
	z = atleast1d(z).([]complex128)
	p = atleast1d(p).([]complex128)
	degree := relativeDegree(z, p)

	zLp := make([]complex128, len(z))
	pLp := make([]complex128, len(p))

	for i, zi := range z {
		zLp[i] = zi * complex(bw/2, 0)
	}

	for i, pi := range p {
		pLp[i] = pi * complex(bw/2, 0)
	}

	zBp := make([]complex128, len(z)*2)
	pBp := make([]complex128, len(p)*2)

	for i, zi := range zLp {
		diff2 := getRootDiff(zi, wo)
		zBp[i] = zi + diff2
		zBp[i+len(zLp)] = zi - diff2
	}

	for i, pi := range pLp {
		diff2 := getRootDiff(pi, wo)
		pBp[i] = pi + diff2
		pBp[i+len(pLp)] = pi - diff2
	}

	for i := 0; i < degree; i++ {
		zBp = append(zBp, 0)
	}

	kBp := k * math.Pow(bw, float64(degree))

	return zBp, pBp, kBp
}

func lp2bsZpk(z []complex128, p []complex128, k float64, wo float64, bw float64) ([]complex128, []complex128, float64) {
	z = atleast1d(z).([]complex128)
	p = atleast1d(p).([]complex128)
	degree := relativeDegree(z, p)

	zHp := make([]complex128, len(z))
	pHp := make([]complex128, len(p))
	zNeg := make([]complex128, len(z))
	pNeg := make([]complex128, len(p))

	for i, zi := range z {
		zHp[i] = complex(bw/2, 0) / zi
		zNeg[i] = complex(-1, 0) * zi
	}

	for i, pi := range p {
		pHp[i] = complex(bw/2, 0) / pi
		pNeg[i] = complex(-1, 0) * pi
	}

	zBs := make([]complex128, len(z)*2)
	pBs := make([]complex128, len(p)*2)

	for i, zi := range zHp {
		diff2 := getRootDiff(zi, wo)
		zBs[i] = zi + diff2
		zBs[i+len(zHp)] = zi - diff2
	}

	for i, pi := range pHp {
		diff2 := getRootDiff(pi, wo)
		pBs[i] = pi + diff2
		pBs[i+len(pHp)] = pi - diff2
	}

	for i := 0; i < degree; i++ {
		zBs = append(zBs, complex(0, wo))
	}
	for i := 0; i < degree; i++ {
		zBs = append(zBs, complex(0, -wo))
	}

	kBs := k * real(complexVectorProductSum(zNeg)/complexVectorProductSum(pNeg))

	return zBs, pBs, kBs
}

func bilinearZpk(z, p []complex128, k float64, fs float64) ([]complex128, []complex128, float64) {
	z = atleast1d(z).([]complex128)
	p = atleast1d(p).([]complex128)
	degree := relativeDegree(z, p)

	fs2 := complex(2*fs, 0)

	zZ := make([]complex128, len(z))
	pZ := make([]complex128, len(p))
	fs2z := make([]complex128, len(z))
	fs2p := make([]complex128, len(p))

	for i, zi := range z {
		zZ[i] = (fs2 + zi) / (fs2 - zi)
		fs2z[i] = fs2 - zi
	}
	for i, pi := range p {
		pZ[i] = (fs2 + pi) / (fs2 - pi)
		fs2p[i] = fs2 - pi
	}

	for i := 0; i < degree; i++ {
		zZ = append(zZ, -1)
	}

	kZ := k * real(complexVectorProductSum(fs2z)/complexVectorProductSum(fs2p))

	return zZ, pZ, kZ
}

func zpk2Tf(z, p []complex128, k float64) ([]float64, []float64) {
	z = atleast1d(z).([]complex128)
	b := vectors.MultiplyBy(Poly(z), k)
	a := atleast1d(Poly(p)).([]float64)

	return b, a
}
