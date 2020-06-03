package fftFirst;

public class FFTtest {
	public static void main(String... args) {
		FFT fft = new FFT();
		
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		
		
		fft.fft(16,testData,2.,16,-1);
	}
}
