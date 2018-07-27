package Problem1;

import java.util.Arrays;

public class CDF24 {
	static final double SQRT_OF_3 = Math.sqrt(3);
	protected final static double denom = 4 * Math.sqrt(2);

	// forward transform smoothing coefficients
	protected final static double h0 = (1 + SQRT_OF_3)/denom;
	protected final static double h1 = (3 + SQRT_OF_3)/denom; 
	protected final static double h2 = (3 - SQRT_OF_3)/denom; 
	protected final static double h3 = (1 - SQRT_OF_3)/denom;

	// forward transform wavelet coefficients
	protected final static double g0 =  h3;
	protected final static double g1 = -h2;
	protected final static double g2 =  h1;
	protected final static double g3 = -h0;

	// Inverse transform coefficients for smoothed values
	protected final static double Ih0 = h2;
	protected final static double Ih1 = g2;
	protected final static double Ih2 = h0;
	protected final static double Ih3 = g0;

	// Inverse transform for wavelet values
	protected final static double Ig0 = h3;
	protected final static double Ig1 = g3;
	protected final static double Ig2 = h1;
	protected final static double Ig3 = g1;


	public static double[] ordDWTForNumIters(double[] signal, int num_iters) {
		int length = signal.length;
		if (length < 4 || !Utils.isPowerOf2(length)) {
			System.out.println("CDF can't be done: Signal length is either less than 4 or not power of 2");
			return null;
		}

		int i, j, half;
		double tmp[] = new double[length];
		while(length >= 4)  {
			half = length/2;
			for(i = 0, j = 0; j < length-3; i += 1, j += 2) {
				tmp[i] = h0*signal[j] + h1*signal[j+1] + h2*signal[j+2] + h3*signal[j+3];
				tmp[half+i] = g0*signal[j] + g1*signal[j+1] + g2*signal[j+2] + g3*signal[j+3];
			}

			tmp[i] = h0*signal[length-2] + h1*signal[length-1] + h2*signal[0] + h3*signal[1];
			tmp[half+i] = g0*signal[length-2] + g1*signal[length-1] + g2*signal[0] + g3*signal[1];

			signal = Arrays.copyOf(tmp, tmp.length);
			//tmp = null;
			length /= 2;
		}
		return signal;
	}

	public static double[] ordInvDWTForNumIters(double[] transformedSignal, int num_iters) {
		int length = transformedSignal.length;
		if (length < 4 || !Utils.isPowerOf2(length)) {
			System.out.println("CDF inverse can't be done: Signal length is either less than 4 or not power of 2");
			return null;
		}

		int transform_length = length / (1 << (num_iters-1));

		while (transform_length <= length) {
			int mid = transform_length/2;

			double[] temp = new double[transform_length];

			temp[0] = Ih0*transformedSignal[mid-1] + Ih1*transformedSignal[transform_length-1] + Ih2*transformedSignal[0] + Ih3*transformedSignal[mid];
			temp[1] = Ig0*transformedSignal[mid-1] + Ig1*transformedSignal[transform_length-1] + Ig2*transformedSignal[0] + Ig3*transformedSignal[mid];

			int i = 0, j = 2;

			while(i < mid-1) {               
				temp[j] = Ih0*transformedSignal[i] + Ih1*transformedSignal[mid+i] + Ih2*transformedSignal[i+1] + Ih3*transformedSignal[mid+i+1];
				temp[j+1] = Ig0*transformedSignal[i] + Ig1*transformedSignal[mid+i] + Ig2*transformedSignal[i+1] + Ig3*transformedSignal[mid+i+1];
				i += 1; j += 2;
			}

			transformedSignal = Arrays.copyOf(temp, temp.length);
			transform_length *= 2;
		}
		
		return transformedSignal;
	}
}