package Filtering

import (
	"errors"
	"github.com/geoport/numpy4go/vectors"
	"math"
)

func iirFilter(N int, Wn []float64, btype, ftype string) ([]float64, []float64) {
	var z, p []complex128
	var k float64

	if ftype == "butterworth" {
		z, p, k = buttap(N)
	} else if ftype == "cheby1" {
		z, p, k = cheb1ap(N)
	} else {
		panic("Filter type not supported")
	}

	fs := 2.
	WnPi := vectors.MultiplyBy(Wn, math.Pi/fs)
	tanWn := vectors.Apply(WnPi, math.Tan)
	warped := vectors.MultiplyBy(tanWn, fs*2)

	if btype == "lowpass" {
		z, p, k = lp2lpZpk(z, p, k, warped)
	} else if btype == "highpass" {
		z, p, k = lp2hpZpk(z, p, k, warped)
	} else {
		bw := warped[1] - warped[0]
		wo := math.Sqrt(warped[0] * warped[1])
		if btype == "bandpass" {
			z, p, k = lp2bpZpk(z, p, k, wo, bw)
		} else if btype == "bandstop" {
			z, p, k = lp2bsZpk(z, p, k, wo, bw)
		} else {
			panic("Filter type not supported")
		}
	}

	z, p, k = bilinearZpk(z, p, k, fs)

	return zpk2Tf(z, p, k)
}

func linearFilter(x []float64, b []float64, a []float64) []float64 {
	N := len(a) - 1 // Degree of the denominator
	// Initialize state variables
	d := make([]float64, N)

	// Apply the filter
	var y []float64
	for _, xn := range x {
		yn := b[0]*xn + d[0]
		for i := 1; i < N; i++ {
			d[i-1] = b[i]*xn - a[i]*yn + d[i]
		}
		d[N-1] = b[N]*xn - a[N]*yn
		y = append(y, yn)
	}

	return y
}

func checkInput(
	signal []float64, cornerFreqs []float64, filterOrder int, btype string, ffunc string, timeStep float64,
) error {
	if len(signal) == 0 {
		return errors.New("signal is empty")
	}

	if len(cornerFreqs) == 0 {
		return errors.New("corner frequencies are empty")
	}

	if filterOrder < 1 {
		return errors.New("filter order must be a positive integer")
	}

	if vectors.Contains([]string{"lowpass", "highpass", "bandpass", "bandstop"}, btype) == false {
		return errors.New("filter type not supported")
	}

	if ffunc != "butterworth" && ffunc != "cheby1" {
		return errors.New("filter function not supported")
	}

	if timeStep <= 0 {
		return errors.New("time step must be a positive number")
	}

	Wn := setCutoffFrequencies(cornerFreqs, timeStep)

	for _, wn := range Wn {
		if wn <= 0 || wn >= 1 {
			return errors.New("cutoff frequencies must be between 0 and 1")
		}
	}

	if (btype == "bandpass" || btype == "bandstop") && len(cornerFreqs) != 2 {
		return errors.New("corner frequencies must be a vector of length 2")
	}

	if (btype == "lowpass" || btype == "highpass") && len(cornerFreqs) != 1 {
		return errors.New("corner frequencies must be a vector of length 1")
	}

	if len(cornerFreqs) == 2 && cornerFreqs[0] >= cornerFreqs[1] {
		return errors.New("corner frequencies must be in ascending order")
	}

	return nil
}

// FilterSignal applies a digital filter to a given signal.
//
// The function takes the following parameters:
//   - signal: A slice of float64 values representing the input signal.
//   - cornerFreqs: A slice of float64 values specifying the corner frequencies of the filter.
//   - filterOrder: An integer representing the order of the filter. (>= 1)
//   - btype: A string indicating the type of the filter ("lowpass", "highpass", "bandpass", or "bandstop").
//   - ffunc: A string specifying the filter function to use ("butterworth" or "cheby1").
//   - timeStep: A float64 value representing the time step between consecutive samples of the signal.
//
// The function returns an error and a filtered signal as output. If an error occurs during the filtering process,
// the error parameter will be non-nil, and the filtered signal will be nil. Otherwise, the error parameter will be nil,
// and the filtered signal will contain the result of the filtering operation.
//
// Example usage:
//
//	signal := []float64{1.2, 3.5, 2.1, 4.7, 2.9}
//	cornerFreqs := []float64{10.0, 20.0}
//	filterOrder := 4
//	btype := "lowpass"
//	ffunc := "butterworth"
//	timeStep := 0.01
//
//	err, filteredSignal := FilterSignal(signal, cornerFreqs, filterOrder, btype, ffunc, timeStep)
//	if err != nil {
//	    fmt.Println("Error:", err)
//	    return
//	}
//	fmt.Println("Filtered signal:", filteredSignal)
//
// The function internally checks the validity of the input parameters, performs the necessary calculations,
// and applies the specified filter to the input signal. The result is a filtered signal that reduces or eliminates
// certain frequency components based on the chosen filter type and parameters.
//
// The filter functions supported by this implementation are "butterworth" and "cheby1".
// The filter type options are "lowpass", "highpass", "bandpass", and "bandstop".
//
// It's important to note that this function assumes a single-channel signal and does not support multichannel signals.
// Also, the code depends on additional helper functions: checkInput, setCutoffFrequencies, iirFilter, and linearFilter,
// which handle input validation, cutoff frequency calculation, IIR filter coefficient computation, and linear filtering,
// respectively.
//
// The function is part of a larger program or package focused on digital signal processing and filtering operations.
// For a complete understanding of the implementation, it's recommended to review the source code of the dependent functions.
//
// For more information on digital signal processing and filter design, refer to relevant literature and resources.
func FilterSignal(
	signal []float64, cornerFreqs []float64, filterOrder int, btype string, ffunc string, timeStep float64,
) (error, []float64) {
	err := checkInput(signal, cornerFreqs, filterOrder, btype, ffunc, timeStep)
	if err != nil {
		return err, nil
	}
	Wn := setCutoffFrequencies(cornerFreqs, timeStep)

	var b, a []float64

	if ffunc == "butterworth" {
		b, a = iirFilter(filterOrder, Wn, btype, "butterworth")
	} else if ffunc == "cheby1" {
		b, a = iirFilter(filterOrder, Wn, btype, "cheby1")
	} else {
		return errors.New("unsupported filter function"), nil
	}

	filteredSignal := linearFilter(signal, b, a)

	return nil, filteredSignal
}
