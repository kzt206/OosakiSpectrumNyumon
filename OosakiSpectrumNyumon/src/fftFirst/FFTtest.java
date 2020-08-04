package fftFirst;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;


public class FFTtest {
	public static void main(String... args) {

		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

//		FFT fft = new FFT(16, testData, 2);
//		fft.cfft(16,testData,2.,16,-1);
//		fft.ffft(16, testData, 2., 16, -1);


		//20200725
		double[] waveData = new double[1024];
		File file = new File("5993wave.txt");
		try(BufferedReader br  = new BufferedReader(new FileReader(file))){
			String text;
			int i = 0;
			while((text = br.readLine()) != null) {
//				System.out.println(text);
				if(i<1024) {
					waveData[i] = Double.parseDouble(text);
				}
				i++;
			}
		}catch (Exception e) {
			e.printStackTrace();
		}finally {
		
		}
		
//		int i=0;
//		for(double d :waveData) {
//			System.out.print(i + " ");
//			System.out.println(d);
//			i++;
//		}
		
		FFT fft2 = new FFT(1024, waveData, 100);
		fft2.ffft(1024, waveData, 100., 1024, -1);
		

		
		//20200725 end
		
		
	}
}
