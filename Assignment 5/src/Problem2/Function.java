package Problem2;

public class Function {
	public double v(double x) {
		return 0;
	}

	public double[] generateRangeInterval(double[] domain_values) {
		return null;
	}
	
	public static void addRandomNoiseToSignal(double[] signal) {
		for(int i = 0; i < signal.length; i++) {
			signal[i] += Math.random()/2.0;
		}
	}
}