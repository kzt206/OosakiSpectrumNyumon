package fftFirst;

public class FFTtest {
	public static void main(String... args) {
		
		
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		FFT fft = new FFT(16,testData,2);
		
//		fft.cfft(16,testData,2.,16,-1);

		fft.ffft(16,testData,2.,16,-1);

	
	}
}
