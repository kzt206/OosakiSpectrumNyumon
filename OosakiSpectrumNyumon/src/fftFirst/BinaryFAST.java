package fftFirst;

public class BinaryFAST {
	public static void main(String... args) {

//		String[] data = { "000", "001", "010", "011", "100", "101", "110", "111" };
		double[] data = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };
		Complex[] complex = new Complex[data.length];
		for (int i = 0; i < data.length; i++) {
			complex[i] = new Complex(data[i], 0.);
		}

		int N = 16;
		int i = 1;
		int j = 1;
		int m = N / 2;
		double deltaT = 0.5;

		// binary
		for (i = 1; i < N + 1; i++) {
//			int j = 1;
			if (i >= j) {
				// goto 110 
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
			do {
				if (j <= m) { // 120
					// j = j + m;
					break;
					// goto 130
				} else {
					j = j - m;
					m = m / 2;
				}
			} while (m >= 2);
			j = j + m; // 130

		}

		// print out
		for (double s : data) {
			System.out.println(s);
		}

		System.out.println();

//		System.out.println(Math.cos(3.1415 * 5));
//		Complex testc=new Complex(5.,1);
//		System.out.println(testc.cexp());
//		
//		System.out.println();

		int IND = -1; // -1:Fourie Transform 1:Inverse
		int kmax = 1;
		while (kmax < N) {
			int istep = kmax * 2;
			for (int k = 1; k < kmax + 1; k++) {
				Complex theta = new Complex(0., Math.PI * IND * (k - 1) / kmax);
				for (int ii = k; ii <= N; ii += istep) {
					int jj = ii + kmax;
					Complex tmp = complex[jj - 1].multiply(theta.cexp());
					complex[jj - 1] = complex[ii - 1].diff(tmp);
					complex[ii - 1] = complex[ii - 1].add(tmp);

				}

			}
			kmax = istep;
		}

		// print out
		for (int ii = 0; ii < complex.length; ii++) {
			complex[ii] = complex[ii].divide(N);
			System.out.printf("%7.3f, %7.3f, %7.3f \n", complex[ii].real(), complex[ii].image(), complex[ii].abs()); // 有限複素フーリエ係数
		}

		System.out.println();
		// print out
		for (int ii = 0; ii < complex.length / 2 + 1; ii++) {
			double a = 2 * complex[ii].real();
			double b = -2 * complex[ii].image();
			double amp = 2 * complex[ii].abs();
			double phase = Math.atan(-1*b/a)/Math.PI*180.;
			double power;
			if(ii == 0 || ii == complex.length/2) {
				power = amp*amp/4;
			}else {
				power = amp*amp/2;
			}
			double fas = amp *N*deltaT/2;
			System.out.printf("f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3f, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii / (N * deltaT), a,b,amp,phase,fas,power); // 有限フーリエ係数
		}

	}
}
