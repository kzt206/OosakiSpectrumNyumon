package fftFirst;

public class BinaryTest {
	public static void main(String... args) {
		
		String[] data = { "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
				"1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"};
		
		int N = 16;
		int i = 1;
		int j = 1;
		int m = N / 4;
		
		// binary 
		for (i = 1; i < N + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				String temp = data[j - 1];
				data[j - 1] = data[i - 1];
				data[i - 1] = temp;
			}

			m = N / 2; // 110
			do {
				if (j <= m) { // 120
					//j = j + m;
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
		for(String s : data) {
			System.out.println(s);
		}
		

	}
}
