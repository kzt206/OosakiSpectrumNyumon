package fftFirst;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;


public class ElCentro {
	public static void main(String... args) {

		// Input data El Centro Wave NS test data
		int n = 1024;
		double dt = 0.02; // Sampling 50Hz?
		double lower = 1;
		double upper = 7;
		double alpha = 1;
		double[] time = new double[n];
		double[] wave = new double[n];

		double[] waveData = new double[n];
		File file = new File("86-El_Centro.csv");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			int nLine = 0;
			int i = 0;
			while ((line = br.readLine()) != null) {
				if (nLine >= 3) {
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

		System.out.println();
//		double[] filterData = filteringBPF(waveData, n, dt, lower, upper, alpha);

//		// output data
//		File outfile = new File("El_Centro_filtered.txt");
//		try (BufferedWriter bw = new BufferedWriter(new FileWriter(outfile))) {
//			String text;
//			for (int i = 0; i < n; i++) {
//				text = String.valueOf(filterData[i]);
//				bw.write(text);
//				bw.newLine();
//			}
//		} catch (Exception e) {
//			e.printStackTrace();
//		} finally {
//
//		}
	}

	
	// OOSAKI p.235 FPAC(Fourier spectrum,Power spectrum,Auto Correlation)
	// ind 100 retun F,ind 010 return G, ind 001 return R
	public static double[] fpac(int n, double[] x, int nd1, double dt, int ind, int nd2) {
		double[] x2 = new double[nd1];
		double[] f = new double[nd2];  // Fourier spectrum
		double[] g = new double[nd2];  // power spectrum
		double[] r = new double[nd2];  // Auto Correlation

		// Initialization
		
		Complex[] c = new Complex[8192];
		for(int i = 0;i<n;i++) {
			c[i].setComplex(x[i], 0.);
		}
		int nt = 2;
		if(nt<n) {
			nt *= 2;
		}
		if(nt != n) {
			for(int i = n;i<nt;i++) {
				c[i].setComplex(0.,0.);
			}
		}
		int nfold = nt/2 + 1;
		double t = nt*dt;
		double df = 1./t;
		
		
		// Fourier transform
		FFT fft1 = new FFT( n, x, 50);
		
		Complex[] c2 = fft1.fast(nt, c, 8192, -1);
		
		// Fourier spectrum
		for(int i = 0; i<nfold ; i++) {
			f[i] = c2[i].abs()*dt;
		}
		
		// Power spectrum
		g[0] = f[0]*f[0]/t;
		for(int i = 1; i<nfold-1 ; i++) {
			g[i] = 2. * f[i]*f[i] /t;
		}
		
		// Auto Correlation
		Complex[] c3 = new Complex[c2.length];
		for(int i = 1;i <nt;i++) {
			c3[i] = c2[i].multiply(c2[i].conjg());
		}
		Complex[] c4 = fft1.fast(nt, c3, 8192, 1);
		double r0 = c[0].real();
		for(int i=0;i<nfold;i++) {
			r[i] = c[i].real()/r0;
		}
		
		//
		switch (ind) {
		case 100: {
			return f;
		}
		case 010: {
			return g;
		}
		case 001: {
			return r;
		}
		}

		return f;

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
