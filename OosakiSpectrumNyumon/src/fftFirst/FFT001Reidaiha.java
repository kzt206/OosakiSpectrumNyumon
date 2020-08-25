package fftFirst;

public class FFT001Reidaiha {
	public static void main(String... args) {
		reidaiha001();
	}

	public static void reidaiha001() {
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		//2020/8/25
		int nOfData = 16;
		int samplingFrequency = 2;
		FFT reidaihaFft = new FFT(nOfData, testData, samplingFrequency);
		
		reidaihaFft.fft(nOfData, testData, samplingFrequency, nOfData, -1);
		
		reidaihaFft.fas(nOfData, testData, samplingFrequency, nOfData);
		
		
		//		// 20200805
//		double dt = 0.5;
//		int n = 16;
//		double lower = 0.4;
//		double upper = 0.7;
//		double alpha = 0.1;
	}
}
