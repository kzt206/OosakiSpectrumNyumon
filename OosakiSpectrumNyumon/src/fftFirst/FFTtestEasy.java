package fftFirst;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;



public class FFTtestEasy {
	public static void main(String... args) {

		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		FFT fft = new FFT(16, testData, 2);

//		fft.cfft(16,testData,2.,16,-1);

		fft.ffft(16, testData, 2., 16, -1);

		System.out.println();
//		System.out.println(fft.sinc(3.1415));

		System.out.println();
		
		//20200802
		double dt = 0.5;
	    int n = 16;
		double lower = 0.4;
		double upper = 0.7;
		double alpha = 0.1;
		
		//元波形のフーリエスペクトル
		double temp = Math.log(n)/Math.log(2);
		int nn = (int) Math.pow(2.,Math.ceil(temp));
		int nfold = nn/2+1;
		double[] xr0 = new double[nn];
		double[] xi0 = new double[nn];
		//データの読み込み
		for(int i=0;i<n;i++) {
			xr0[i] = testData[i];
			xi0[i] = 0.;
		}
		for(int i=n;i<nn;i++) {
			xr0[i] = 0.;
			xi0[i] = 0.;
		}
		//データに対するFFTの計算
		int ind = -1;
		FFT fft2 = new FFT(nn, xr0, 2); // 2 = 1/dt
		double[][] coef2 = fft.ffft(nn, xr0, 2, nn, ind);
		System.out.println();
		for(int i=0;i<nfold;i++) {
			System.out.println(coef2[i][0] + " " + coef2[i][1]);
		}
		
		//フィルター
		double fe1 = lower*dt;
		double fe2 = upper*dt;
		double delta = alpha *dt;
		int j = (int)(3.1/delta+0.5) - 1;
		if(j%2 == 1) {
			j = j+1;
		}
		temp = Math.log(n+j)/Math.log(2);
		nn = (int) Math.pow(2.,Math.ceil(temp));
		nfold = nn/2 + 1;
		
		double[] xr = new double[nn];
		double[] xi = new double[nn];
		double[] yr = new double[nn];
		double[] yi = new double[nn];
		double[] br = new double[nn];
		double[] bi = new double[nn];
		double[] b = new double[j+1];
		double[] w = new double[j+1];
		//データの再読み込み
		for(int i=0;i<n;i++) {
			xr[i] = testData[i];
			xi[i] = 0.;
		}
		for(int i=n;i<nn;i++) {
			xr[i] = 0.;
			xi[i] = 0.;
		}
		//データに対するFFTの計算
		ind = -1;
		FFT fft3 = new FFT(nn, xr, 2); // 2 = 1/dt
		double[][] coef3 = fft.ffft(nn, xr, 2, nn, ind);
		System.out.println();
		w = fft3.hannigWindow(w, j+1);
		System.out.println();
		for(int i =0;i<j+1;i++) {
			System.out.println(w[i]);
		}
		fft3.bpf(fe1, fe2, j, b, w);
		
		
		// end 20200802
		
		
//		//20200725
//		double[] waveData = new double[1024];
//		File file = new File("5993wave.txt");
//		try(BufferedReader br  = new BufferedReader(new FileReader(file))){
//			String text;
//			int i = 0;
//			while((text = br.readLine()) != null) {
//				System.out.println(text);
//				if(i<1024) {
//					waveData[i] = Double.parseDouble(text);
//				}
//				i++;
//			}
//		}catch (Exception e) {
//			e.printStackTrace();
//		}finally {
//		
//		}
//		
//		int i=0;
//		for(double d :waveData) {
//			System.out.print(i + " ");
//			System.out.println(d);
//			i++;
//		}
//		
//		FFT fft2 = new FFT(1024, waveData, 100);
//		fft2.ffft(1024, waveData, 100., 1024, -1);
//		
//		
//		//20200725 end
		
	}
}
