package fourier_spectrum

import (
	"math"

	"github.com/eripe970/go-dsp-utils"
	np "github.com/geoport/numpy4go/vectors"
)

func FourierSpectrum(data []float64, timeStep float64) ([]float64, []float64, []float64) {
	signal := dsp.Signal{SampleRate: 1 / timeStep, Signal: data}
	spectrum, _ := signal.FrequencySpectrum()
	frequency := spectrum.Frequencies

	fourierAmplitudes := spectrum.Spectrum
	T := float64(len(data)) * timeStep
	aRMS2 := np.DividedBy(np.Cumtrapz(np.Pow(data, 2), timeStep, 0), T)
	aRMS := math.Pow(aRMS2[len(aRMS2)-1], 0.5)
	powerAmplitudes := np.DividedBy(np.Pow(fourierAmplitudes, 2), math.Pi*T*math.Pow(aRMS, 2))

	return frequency, fourierAmplitudes, powerAmplitudes
}
