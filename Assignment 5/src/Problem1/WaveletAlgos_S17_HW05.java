package Problem1;

public class WaveletAlgos_S17_HW05 {

	static void displaySignal(double[] sig){
		for(double d: sig)
			System.out.println(d);
		System.out.println();
	}

	static void test01() {
		double[] signal = {1, 2, 3, 4};
		signal = CDF24.ordDWTForNumIters(signal, 1);
		displaySignal(signal);
		CDF24.ordInvDWTForNumIters(signal, 1);
	}

	static void test02() {
		double[] signal = {1, 2, 3, 4, 5, 6, 7, 8};
		signal = CDF24.ordDWTForNumIters(signal, 1);
		CDF24.ordInvDWTForNumIters(signal, 1);
	}

	static void test03() {
		double[] signal = {1, 2, 3, 4, 5, 6, 7, 8};
		signal = CDF24.ordDWTForNumIters(signal, 2);
		displaySignal(signal);
		signal = CDF24.ordInvDWTForNumIters(signal, 2);
		displaySignal(signal);
	}

	public static void main(String[] args) {
		test03();
	}
}