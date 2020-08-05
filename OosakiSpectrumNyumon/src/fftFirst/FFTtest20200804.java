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
		
		fftTest();
	}
	
	public static void fftTest(String... args) {
		
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };
		
		//20200805
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
		double[][] coef2 = fft2.ffft(nn, xr0, 2, nn, ind);
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
		double[][] coef3 = fft3.ffft(nn, xr, 2, nn, ind);
		System.out.println();
		w = fft3.hannigWindow(w, j+1);
		System.out.println("\n<<<<  w  >>>>");
		for(int i =0;i<j+1;i++) {
			System.out.println(w[i]);
		}
		b = fft3.bpf(fe1, fe2, j, b, w);
		System.out.println("\n<<<<  b  >>>>");
		for(int i =0;i<j+1;i++) {
			System.out.println(b[i]);
		}
		
		for(int i=0;i<j;i++) {
			br[i] = b[i];
			bi[i] = 0.;
		}
		for(int i=j+1;i<nn;i++) {
			br[i] = 0.;
			bi[i] = 0.;
		}
		
		//フィルターに対するFFTの計算
		ind = -1;
		FFT fft4 = new FFT(nn, br, 2); // 2 = 1/dt
		System.out.println("\n\n<<fft4>>");
		double[][] coef4 = fft4.ffft(nn, br, 2, nn, ind);
		for(int i =0;i<nn;i++) {
			yr[i] = coef3[i][0] * coef4[i][0] - coef3[i][1] * coef4[i][1];
			yi[i] = coef3[i][1] * coef4[i][0] + coef3[i][0] * coef4[i][1];
		}
		System.out.println("\n<< 7,8 >>");
		for(int i=0;i<nfold;i++) {
			System.out.printf("%8.6f, %10.8f\n", (double)i/nn/dt,Math.sqrt(yr[i]*yr[i] + yi[i]*yi[i]));
		}
		//フーリエ逆変換
		
		
		
		// end 20200805
	}
	
}
