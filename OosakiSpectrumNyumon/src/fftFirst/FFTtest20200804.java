package fftFirst;


public class FFTtest20200804 {
	public static void main(String... args) {
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };
		double temp;
		
		FFT fft = new FFT(16, testData, 2);
		
		double[][] x = fft.fir_fft(16, testData, 2, 0, -1);
	
		for(int i = 0;i<testData.length;i++) {
			temp = Math.sqrt(Math.pow(x[i][0], 2.) + Math.pow(x[i][1], 2.));
			System.out.printf("i:%4d, xr:%11.6f, xi:%11.6f, amp:%11.6f\n",i,x[i][0],x[i][1],temp);
		}
		
	}
}
