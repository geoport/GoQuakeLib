package processing

import (
	"errors"
	np "github.com/geoport/numpy4go/vectors"
)

// BaselineCorrection performs baseline correction on a given signal.
//
// The function takes the following parameters:
//   - signal: A slice of float64 values representing the input signal.
//   - time: A slice of float64 values representing the corresponding time values for each sample in the signal.
//   - order: An integer representing the order of the polynomial to fit for baseline correction.
//
// The function returns a slice of float64 values representing the corrected signal after baseline subtraction.
//
// Example usage:
//
//	signal := []float64{1.2, 3.5, 2.1, 4.7, 2.9}
//	time := []float64{0.1, 0.2, 0.3, 0.4, 0.5}
//	order := 2
//
//	correctedSignal := BaselineCorrection(signal, time, order)
//	fmt.Println("Corrected signal:", correctedSignal)
//
// The function utilizes a polynomial fitting approach to estimate the baseline of the input signal.
// It fits a polynomial of the specified order to the given signal and time values and subtracts the fitted polynomial
// from the original signal to obtain the baseline-corrected signal.
//
// The code depends on an external library or module called "numpy4go" (likely a numerical computation library),
// which provides the functions Polyval and Subtract for polynomial evaluation and subtraction, respectively.
// Please ensure that this library is imported and available for the code to execute successfully.
//
// It's important to note that the baseline correction assumes that the baseline variations in the signal can be approximated
// by a polynomial function of the specified order. The success of the correction heavily relies on the appropriateness
// of the chosen order for the given signal.
//
// The function is part of a larger program or package focused on signal processing or data analysis tasks,
// and it provides a convenient way to perform baseline correction on signals before further analysis or processing.
//
// For more information on baseline correction techniques and polynomial fitting, refer to relevant literature and resources.
func BaselineCorrection(signal []float64, time []float64, order int) (error, []float64) {
	if len(signal) != len(time) {
		return errors.New("Signal and time vectors must be of equal length"), nil
	}
	if len(signal) == 0 {
		return errors.New("Signal vector is empty"), nil
	}
	if order < 1 {
		return errors.New("Order must be a positive integer"), nil
	}

	pred := np.Polyval(time, signal, order, time).([]float64)
	correctedData := np.Substract(signal, pred)
	return nil, correctedData
}
