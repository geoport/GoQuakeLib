package Filtering

import (
	"math/cmplx"
	"reflect"
)

func Poly(roots []complex128) []float64 {
	degree := len(roots)
	coefficients := make([]complex128, degree+1)
	coefficients[0] = complex(1.0, 0.0)

	for _, root := range roots {
		for j := degree; j >= 1; j-- {
			coefficients[j] = coefficients[j] - root*coefficients[j-1]
		}
	}

	return getReal(coefficients)
}

func atleast1d(arr interface{}) interface{} {
	arrValue := reflect.ValueOf(arr)

	// Check if the input is already at least 1-dimensional
	if arrValue.Kind() == reflect.Slice {
		return arr
	}

	// Convert scalar values to a slice of length 1
	sliceType := reflect.SliceOf(arrValue.Type())
	slice := reflect.MakeSlice(sliceType, 1, 1)
	slice.Index(0).Set(arrValue)

	return slice.Interface()
}

func complexVectorProductSum(vector []complex128) complex128 {
	result := complex(1, 0) // Initialize the result as 1 + 0i

	for _, value := range vector {
		result *= value
	}

	return result
}

func getRootDiff(a complex128, wo float64) complex128 {
	wo2 := cmplx.Pow(complex(wo, 0), 2)
	a2 := cmplx.Pow(a, 2)
	diff2 := cmplx.Pow(a2-wo2, 0.5)
	return diff2
}

func relativeDegree(z []complex128, p []complex128) int {
	return len(p) - len(z)
}

func getReal(a []complex128) []float64 {
	reals := make([]float64, len(a))
	for i, ai := range a {
		reals[i] = real(ai)
	}
	return reals
}

func setCutoffFrequencies(cornerFrequencies []float64, timeStep float64) []float64 {
	var cutoff []float64
	nyq := 0.5 / timeStep

	for _, freq := range cornerFrequencies {
		cutoff = append(cutoff, freq/nyq)
	}
	return cutoff
}
