package fftFirst;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class FFTtest20200811 {
	public static void main(String... args) {
//		fftTest();

//		akatani5993();

		// Input data akatani test data
		int n = 1024;
		double dt = 0.01;
		double lower = 1;
		double upper = 7;
		double alpha = 1;
		
		
		double[] waveData = new double[n];
		File file = new File("5993wave.txt");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String text;
			int i = 0;
			while ((text = br.readLine()) != null) {
				if (i < n) {
					waveData[i] = Double.parseDouble(text);
				}
				i++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {

		}

		double[] filterData = filteringBPF(waveData, n, dt,lower,upper,alpha);
		
		
		
		// output data
		File outfile = new File("5993_filtered_wave4.txt");
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(outfile))) {
			String text;
			for (int i = 0; i < n; i++) {
				text = String.valueOf(filterData[i]);
				bw.write(text);
				bw.newLine();
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {

		}
		
	}

	public static void akatani5993(String... args) {

		System.out.println("akatani");

		// 20200811

		// Input data
		int n = 1024;

		double[] waveData = new double[n];
		File file = new File("5993wave.txt");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String text;
			int i = 0;
			while ((text = br.readLine()) != null) {
				if (i < n) {
					waveData[i] = Double.parseDouble(text);
				}
				i++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {

		}

//		int i=0;
//		for(double d :waveData) {
//			System.out.print(i + " ");
//			System.out.println(d);
//			i++;
//		}
//		FFT fft2 = new FFT(1024, waveData, 100);
//		fft2.ffft(1024, waveData, 100., 1024, -1);

		// Shori start
		double dt = 0.01;
		double lower = 1;
		double upper = 7;
		double alpha = 1;

		// 元波形のフーリエスペクトル
		double temp = Math.log(n) / Math.log(2);
		int nn = (int) Math.pow(2., Math.ceil(temp));
		int nfold = nn / 2 + 1;
		double[] xr0 = new double[nn];
		double[] xi0 = new double[nn];
		// データの読み込み
		for (int i = 0; i < n; i++) {
			xr0[i] = waveData[i];
			xi0[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr0[i] = 0.;
			xi0[i] = 0.;
		}
		// データに対するFFTの計算
		int ind = -1;
		FFT fft2 = new FFT(nn, xr0, 2); // 2 = 1/dt
		double[][] coef2 = fft2.fir_fft(nn, xr0, xi0, 2, nn, ind);

		// フィルター
		double fe1 = lower * dt;
		double fe2 = upper * dt;
		double delta = alpha * dt;
		int j = (int) (3.1 / delta + 0.5) - 1;
		if (j % 2 == 1) {
			j = j + 1;
		}
		temp = Math.log(n + j) / Math.log(2);
		nn = (int) Math.pow(2., Math.ceil(temp));
		nfold = nn / 2 + 1;

		double[] xr = new double[nn];
		double[] xi = new double[nn];
		double[] yr = new double[nn];
		double[] yi = new double[nn];
		double[] br = new double[nn];
		double[] bi = new double[nn];
		double[] b = new double[j + 2];
		double[] w = new double[j + 1];
		// データの再読み込み
		for (int i = 0; i < n; i++) {
			xr[i] = waveData[i];
			xi[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr[i] = 0.;
			xi[i] = 0.;
		}
		// データに対するFFTの計算
		ind = -1;
		FFT fft3 = new FFT(nn, xr, 2); // 2 = 1/dt
		double[][] coef3 = fft3.fir_fft(nn, xr, xi, 2, nn, ind);

		w = fft3.hannigWindow(w, j + 1);

		b = fft3.bpf(fe1, fe2, j, b, w);

		for (int i = 0; i < j; i++) {
			br[i] = b[i];
			bi[i] = 0.;
		}
		for (int i = j + 1; i < nn; i++) {
			br[i] = 0.;
			bi[i] = 0.;
		}

		// フィルターに対するFFTの計算
		ind = -1;
		FFT fft4 = new FFT(nn, br, 2); // 2 = 1/dt

		double[][] coef4 = fft4.fir_fft(nn, br, bi, 2, nn, ind);

		for (int i = 0; i < nn; i++) {
			yr[i] = coef3[i][0] * coef4[i][0] - coef3[i][1] * coef4[i][1];
			yi[i] = coef3[i][1] * coef4[i][0] + coef3[i][0] * coef4[i][1];
		}

		// フーリエ逆変換
		ind = 1;
		double[][] coef5 = fft4.fir_fft(nn, yr, yi, 2, nn, ind);
		double temp2;

		// 結果を出力（調整波形）
		System.out.println("<<<<  filtered waveform  >>>>");
		for (int i = 0; i < n; i++) {
			temp2 = coef5[j / 2 - 1 + i][0] / nn;
			System.out.printf("%6.3f %14.10f\n", i * dt, temp2); // OK
		}

		// Shori end

		// output data
		File outfile = new File("5993_filtered_wave3.txt");
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(outfile))) {
			String text;
			for (int i = 0; i < n; i++) {
				temp2 = coef5[j / 2 - 1 + i][0] / nn;
				text = String.valueOf(temp2);
				bw.write(text);
				bw.newLine();
//				System.out.println("wirte");
			}
//			int i = 0;
//			while((text = bw.readLine()) != null) {
//				if(i<1024) {
//					waveData[i] = Double.parseDouble(text);
//				}
//				i++;
//			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {

		}
	}

	public static void fftTest(String... args) {

		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		// 20200805-2
		double dt = 0.5;
		int n = 16;
		double lower = 0.4;
		double upper = 0.7;
		double alpha = 0.1;

		// 元波形のフーリエスペクトル
		double temp = Math.log(n) / Math.log(2);
		int nn = (int) Math.pow(2., Math.ceil(temp));
		int nfold = nn / 2 + 1;
		double[] xr0 = new double[nn];
		double[] xi0 = new double[nn];
		// データの読み込み
		for (int i = 0; i < n; i++) {
			xr0[i] = testData[i];
			xi0[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr0[i] = 0.;
			xi0[i] = 0.;
		}
		// データに対するFFTの計算
		int ind = -1;
		FFT fft2 = new FFT(nn, xr0, 2); // 2 = 1/dt
		double[][] coef2 = fft2.fir_fft(nn, xr0, xi0, 2, nn, ind);

		// フィルター
		double fe1 = lower * dt;
		double fe2 = upper * dt;
		double delta = alpha * dt;
		int j = (int) (3.1 / delta + 0.5) - 1;
		if (j % 2 == 1) {
			j = j + 1;
		}
		temp = Math.log(n + j) / Math.log(2);
		nn = (int) Math.pow(2., Math.ceil(temp));
		nfold = nn / 2 + 1;

		double[] xr = new double[nn];
		double[] xi = new double[nn];
		double[] yr = new double[nn];
		double[] yi = new double[nn];
		double[] br = new double[nn];
		double[] bi = new double[nn];
		double[] b = new double[j + 2];
		double[] w = new double[j + 1];
		// データの再読み込み
		for (int i = 0; i < n; i++) {
			xr[i] = testData[i];
			xi[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr[i] = 0.;
			xi[i] = 0.;
		}
		// データに対するFFTの計算
		ind = -1;
		FFT fft3 = new FFT(nn, xr, 2); // 2 = 1/dt
		double[][] coef3 = fft3.fir_fft(nn, xr, xi, 2, nn, ind);

		w = fft3.hannigWindow(w, j + 1);

		b = fft3.bpf(fe1, fe2, j, b, w);

		for (int i = 0; i < j; i++) {
			br[i] = b[i];
			bi[i] = 0.;
		}
		for (int i = j + 1; i < nn; i++) {
			br[i] = 0.;
			bi[i] = 0.;
		}

		// フィルターに対するFFTの計算
		ind = -1;
		FFT fft4 = new FFT(nn, br, 2); // 2 = 1/dt

		double[][] coef4 = fft4.fir_fft(nn, br, bi, 2, nn, ind);

		for (int i = 0; i < nn; i++) {
			yr[i] = coef3[i][0] * coef4[i][0] - coef3[i][1] * coef4[i][1];
			yi[i] = coef3[i][1] * coef4[i][0] + coef3[i][0] * coef4[i][1];
		}

		// フーリエ逆変換
		ind = 1;
		double[][] coef5 = fft4.fir_fft(nn, yr, yi, 2, nn, ind);

		// 結果を出力（調整波形）
		System.out.println("<<<<  filtered waveform  >>>>");
		for (int i = 0; i < n; i++) {
			System.out.printf("%10.6f\n", coef5[j / 2 - 1 + i][0] / nn); // OK
		}

		// end 20200805-2
	}

	public static double[] filteringBPF(double[] waveData, int n, double dt, double lower, double upper, double alpha) {

		double[] filteredData = new double[waveData.length];
		System.out.println("filter start");

		// 元波形のフーリエスペクトル
		double temp = Math.log(n) / Math.log(2);
		int nn = (int) Math.pow(2., Math.ceil(temp));
		int nfold = nn / 2 + 1;
		double[] xr0 = new double[nn];
		double[] xi0 = new double[nn];
		// データの読み込み
		for (int i = 0; i < n; i++) {
			xr0[i] = waveData[i];
			xi0[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr0[i] = 0.;
			xi0[i] = 0.;
		}
		// データに対するFFTの計算
		int ind = -1;
		FFT fft2 = new FFT(nn, xr0, 2); // 2 = 1/dt
		double[][] coef2 = fft2.fir_fft(nn, xr0, xi0, 2, nn, ind);

		// フィルター
		double fe1 = lower * dt;
		double fe2 = upper * dt;
		double delta = alpha * dt;
		int j = (int) (3.1 / delta + 0.5) - 1;
		if (j % 2 == 1) {
			j = j + 1;
		}
		temp = Math.log(n + j) / Math.log(2);
		nn = (int) Math.pow(2., Math.ceil(temp));
		nfold = nn / 2 + 1;

		double[] xr = new double[nn];
		double[] xi = new double[nn];
		double[] yr = new double[nn];
		double[] yi = new double[nn];
		double[] br = new double[nn];
		double[] bi = new double[nn];
		double[] b = new double[j + 2];
		double[] w = new double[j + 1];
		// データの再読み込み
		for (int i = 0; i < n; i++) {
			xr[i] = waveData[i];
			xi[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr[i] = 0.;
			xi[i] = 0.;
		}
		// データに対するFFTの計算
		ind = -1;
		FFT fft3 = new FFT(nn, xr, 2); // 2 = 1/dt
		double[][] coef3 = fft3.fir_fft(nn, xr, xi, 2, nn, ind);

		w = fft3.hannigWindow(w, j + 1);

		b = fft3.bpf(fe1, fe2, j, b, w);

		for (int i = 0; i < j; i++) {
			br[i] = b[i];
			bi[i] = 0.;
		}
		for (int i = j + 1; i < nn; i++) {
			br[i] = 0.;
			bi[i] = 0.;
		}

		// フィルターに対するFFTの計算
		ind = -1;
		FFT fft4 = new FFT(nn, br, 2); // 2 = 1/dt

		double[][] coef4 = fft4.fir_fft(nn, br, bi, 2, nn, ind);

		for (int i = 0; i < nn; i++) {
			yr[i] = coef3[i][0] * coef4[i][0] - coef3[i][1] * coef4[i][1];
			yi[i] = coef3[i][1] * coef4[i][0] + coef3[i][0] * coef4[i][1];
		}

		// フーリエ逆変換
		ind = 1;
		double[][] coef5 = fft4.fir_fft(nn, yr, yi, 2, nn, ind);
		double temp2;

		// 結果を出力（調整波形）
		System.out.println("<<<<  filtered waveform  >>>>");
		for (int i = 0; i < n; i++) {
			temp2 = coef5[j / 2 - 1 + i][0] / nn;
			filteredData[i] = temp2;
			System.out.printf("%6.3f %14.10f\n", i * dt, temp2); // OK
		}

		return filteredData;

	}
}
