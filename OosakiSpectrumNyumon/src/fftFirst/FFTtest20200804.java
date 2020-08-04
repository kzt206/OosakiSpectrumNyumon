package fftFirst;

public class FFTtest20200804 {
	public static void main(String... args) {
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		FFT fft = new FFT(16, testData, 2);
		
		fft.fir_fft(16, testData, 2, 0, -1);
		
	}
}
