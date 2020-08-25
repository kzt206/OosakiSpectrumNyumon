package fftFirst;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class FFT002ElCentro {
	public static void main(String... args) {

		// Input data El Centro Wave NS test data
		int nOfData = 1024;
		int samplingFrequency = 50;
		double dt = 0.02; // Sampling 50Hz?
		double[] time = new double[nOfData];
		double[] wave = new double[nOfData];

		File file = new File("86-El_Centro.csv");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			int nLine = 0;
			int i = 0;
			while ((line = br.readLine()) != null) {
				if (nLine >= 3 && i < nOfData) {
					String[] data = line.split(",", 0);
					time[i] = Double.parseDouble(data[0]);
					wave[i] = Double.parseDouble(data[1]);
					// waveData[i] = Double.parseDouble(text);
					for (String elem : data) {
						System.out.print(elem + ",");

					}
					System.out.println();
					i++;
				}
				nLine++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {

		}

		FFT reidaihaFft = new FFT(nOfData, wave, samplingFrequency);

		reidaihaFft.fft(nOfData, wave, samplingFrequency, nOfData, -1);
		double[][] fas = reidaihaFft.fas(nOfData, wave, samplingFrequency, nOfData);

		
		
	}
}
