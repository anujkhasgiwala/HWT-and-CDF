package Problem2;

import Problem1.CDF24;
import Problem1.Utils;

public class Ripples_F_ex_4_4_p34 extends Function {
	public static enum DWT {CDF24, HWT};
	public static enum COEFF {S, D};

	static final int D8_START_1024 = 128;
	static final int D8_END_1024 = 255;

	static final int D7_START_1024 = 64;
	static final int D7_END_1024 = 127;

	static final int D6_START_1024 = 32;
	static final int D6_END_1024 = 63;

	static final int S6_START_1024 = 0;
	static final int S6_END_1024 = 31;
	
	static double[] sDomain = partition(0, 511, 1);
	static double[] sRangeFig_4_12_p34 = new double[512];
	static Ripples_F_ex_4_4_p34 sRipples_F_p34 = new Ripples_F_ex_4_4_p34();
	static double[] sRange  = new double[512];
	
	public final static double FSNORM = Math.sqrt(2);
    public final static double FDNORM = 1/FSNORM;

	@Override
	public double v(double t) {
		if (0.0 <= t && t < 0.25) {
			return Math.sin(4*Math.PI*t);
		} else if (0.25 <= t && t < 0.75) {
			return 1 + Math.sin(4*Math.PI*t);
		} else if (0.75 <= t && t <= 1) {
			return Math.sin(4*Math.PI*t);
		} else {
			return 0;
		}
	}

	static void multires_fig_4_13_cdf_p34(String message, int range_start, int range_end) {
		multires_fig_4_13_p34_aux(message, DWT.CDF24, 6, range_start, range_end);
	}

	static void multires_fig_4_13_hwt_p34(String message, int range_start, int range_end) {
		multires_fig_4_13_p34_aux(message, DWT.HWT, 6, range_start, range_end);
	}

	static void multires_fig_4_13_p34_aux(String message, DWT dwt, int num_iters, int range_start, int range_end) {
		for(int i = 0; i < 1024; i++)  {
			//sRangeFig_4_12_p34[i] = sRipples_F_p34.v(sDomainFig_4_12_p34[i]);
		}

		for(int i = 1; i < 1024; i += 32) {
			sRangeFig_4_12_p34[i] += 2;
		}

		//final int NUM_ITERS = 6;
		forwardDWTForNumIters(sRangeFig_4_12_p34, dwt, num_iters, range_start, range_end);

		//CDF44.orderedDWTForNumIters(sRangeFig_4_12_p33, 6, false);

		double[] signal = new double[sRangeFig_4_12_p34.length];
		for(int i = 0; i < 1024; i++) {
			if ( i >= range_start && i <= range_end ) {
				signal[i] = sRangeFig_4_12_p34[i];
			}
			else {
				signal[i] = 0;
			}
		}

		System.out.println("=========================");
		System.out.println(message);
		//display_signal(signal);
		System.out.println("Inversed Signal");
		System.out.println("=========================");

		inverseDWTForNumIters(signal, dwt, num_iters);
		//CDF44.orderedInverseDWTForNumIters(signal, 6, false);

		//display_signal(signal);
		System.out.println("=========================");
	}

	static void fig_4_13_D8_CDF44_p34() {
		multires_fig_4_13_cdf_p34("Fig. 4.13, 06-06-07-D8-09-010, CDF(4,4), p. 33", D8_START_1024, D8_END_1024);
	}

	static void fig_4_13_D8_HWT_p34() {
		multires_fig_4_13_hwt_p34("Fig. 4.13, 06-06-07-D8-09-010, HWT, p. 33", D8_START_1024, D8_END_1024);
	}

	static void fig_4_13_D7_CDF44_p34() {
		multires_fig_4_13_cdf_p34("Fig. 4.13, 06-06-D7-08-09-010, CDF(4,4), p. 33", D7_START_1024, D7_END_1024);
	}

	static void fig_4_13_D7_HWT_p34() {
		multires_fig_4_13_hwt_p34("Fig. 4.13, 06-06-D7-08-09-010, HWT, p. 33", D7_START_1024, D7_END_1024);
	}

	static void fig_4_13_D6_CDF44_p34() {
		multires_fig_4_13_cdf_p34("Fig. 4.13, 06-D6-07-08-09-010, CDF(4,4), p. 33", D6_START_1024, D6_END_1024);
	}

	static void fig_4_13_D6_HWT_p34() {
		multires_fig_4_13_hwt_p34("Fig. 4.13, 06-D6-07-08-09-010, HWT, p. 33", D6_START_1024, D6_END_1024);
	}

	static void fig_4_13_S6_CDF44_p34() {
		multires_fig_4_13_cdf_p34("Fig. 4.13, S6-06-07-08-09-010, CDF(4,4), p. 33", S6_START_1024, S6_END_1024);
	}

	static void fig_4_13_S6_HWT_p34() {
		multires_fig_4_13_hwt_p34("Fig. 4.13, S6-06-07-08-09-010, HWT, p. 33", S6_START_1024, S6_END_1024);
	}
	
	public static double[] partition(double from, double upto, double step) {
        if (upto <= from) return null;
        int n = (int)((upto - from)/step) + 1;
        double[] interval = new double[n];
        int i;
        double curr;
        for(i = 0, curr = from; i < n; i++, curr += step) {
            interval[i] = curr/upto;
        }
        return interval;
    }
	
	public static void forwardDWTForNumIters(double[] signal, DWT dwt, int num_iters, int range_start, int range_end) {
        if ( dwt == DWT.CDF24 ) {
            CDF24.ordDWTForNumIters(signal, num_iters);
        }
        else if ( dwt == DWT.HWT ) {
            orderedNormalizedFastHaarWaveletTransformForNumIters(signal, num_iters);
        }
        else {
            throw new IllegalArgumentException("Illegal dwt value: " + dwt);
        }
    }
	
	public static void inverseDWTForNumIters(double[] signal, DWT dwt, int num_iters) {
        if ( dwt == DWT.CDF24 ) {
            CDF24.ordInvDWTForNumIters(signal, num_iters);
        }
        else if ( dwt == DWT.HWT ) {
            orderedNormalizedFastInverseHaarWaveletTransformForNumIters(signal, num_iters);
        }
        else {
            throw new IllegalArgumentException("Illegal dwt value: " + dwt);
        }
    }
	
	public static void orderedNormalizedFastHaarWaveletTransformForNumIters(double[] sample, int num_iters) {
        final int n = sample.length;
        if ( !Utils.isPowerOf2(n) ) return;
        final int NUM_SWEEPS = (int) (Math.log(n) / Math.log(2.0));
        if ( num_iters > NUM_SWEEPS ) return;
        double acoeff, ccoeff;
        if ( NUM_SWEEPS == 1 ) {
            acoeff = FSNORM * (sample[0] + sample[1])/2.0;
            ccoeff = FDNORM * (sample[0] - sample[1]);
            sample[0] = acoeff;
            sample[1] = ccoeff;
            return;
        }
        double[] acoeffs;
        double[] ccoeffs;
        for (int SWEEP_NUM = 1; SWEEP_NUM <= num_iters; SWEEP_NUM++) {
            int size = (int) Math.pow(2.0, (double) (NUM_SWEEPS - SWEEP_NUM)); 
            acoeffs = new double[size]; // where we place a-coefficients
            ccoeffs = new double[size]; // where we place c-coefficients
            int ai = 0; // index over acoeffs
            int ci = 0; // index over ccoeffs
            int end = ((int) Math.pow(2.0, (double) (NUM_SWEEPS - SWEEP_NUM + 1))) - 1;
            for (int i = 0; i <= end; i += 2) {
                acoeffs[ai++] = FSNORM * (sample[i] + sample[i + 1])/2.0;
                ccoeffs[ci++] = FDNORM * (sample[i] - sample[i + 1]);
            }
            
            for (int i = 0; i < size; i++) {
                sample[i] = acoeffs[i];
                sample[i + size] = ccoeffs[i];
            }
        }
    }
	
	public static void orderedNormalizedFastInverseHaarWaveletTransformForNumIters(double[] sample, int num_iters) {
        int n = sample.length;
        if (n < 2 || !Utils.isPowerOf2(n)) {
            return;
        }
        n = (int) (Math.log(n) / Math.log(2.0));
        if (num_iters > n) {
            return;
        }
        double a0 = 0;
        double a1 = 0;
        double[] restored_vals = null;
        int GAP = (int)(Math.pow(2.0, n-num_iters));
        for (int L = 1; L <= num_iters; L++) {
            restored_vals = null;
            restored_vals = new double[2 * GAP]; // restored values at level L
            for (int i = 0; i < GAP; i++) {
                double d = FSNORM * sample[GAP + i];
                double s = FDNORM * sample[i];
                a0 = s + d/2;
                a1 = s - d/2;
                restored_vals[2 * i] = a0;
                restored_vals[2 * i + 1] = a1;
            }
            System.arraycopy(restored_vals, 0, sample, 0, 2 * GAP);
            GAP *= 2;
        }
    }
	
	public static void main(String[] args) {
		generateRange(512);
		addRandomNoiseToSignal(sRange);
        computeSTDsForDWTMultiresCoeffs(sRange, DWT.HWT, 3, true);
	}
	
	public static void generateRange(int size) {
        Ripples_F_ex_4_4_p34 curve1 = new Ripples_F_ex_4_4_p34();
        int i = 0;
        for(double x: sDomain) {
            sRange[i++] = curve1.v(x);
        }
    }
	
	public static void computeSTDsForDWTMultiresCoeffs(double[] signal, DWT dwt, int num_iters, boolean display_signal_flag) {
        forwardDWTForNumIters(signal, dwt, num_iters);
        double[] reconstructed_signal = null; 
        final int sig_len = signal.length;
        int start = 0; int end = 0;
        
        for(int scale = 1; scale <= num_iters; scale++) {
            start = getDCoeffsStart(sig_len, scale);
            end   = getDCoeffsEnd(sig_len, scale);
            reconstructed_signal = new double[signal.length];
            multiresSignalReconstruct(signal, reconstructed_signal, dwt, COEFF.D, num_iters, scale);
            System.out.println("start = " + start + "; end = " + end);
            System.out.println("STD of D[" + (num_iters - scale) + "] = " + Utils.computeCorrectedSTDInRange(signal, start, end));
            if ( display_signal_flag ) Utils.displaySignalRange(signal, start, end);
            multiresSignalReconstruct(signal, reconstructed_signal, dwt, COEFF.S, num_iters, scale);
            start = 0;
            end   = getSCoeffsEnd(sig_len, scale);
            System.out.println("start = " + 0 + "; end = " + end);
            System.out.println("STD of S[" + (num_iters - scale) + "] = " + Utils.computeCorrectedSTDInRange(signal, start, end));
            if ( display_signal_flag) Utils.displaySignalRange(signal, start, end);
            System.out.println("*********************");
        }
        System.out.println();
    }
	
	public static void forwardDWTForNumIters(double[] signal, DWT dwt, int num_iters) {
        if ( dwt == DWT.CDF24 ) {
            CDF24.ordDWTForNumIters(signal, num_iters);
        }
        else if ( dwt == DWT.HWT ) {
            orderedNormalizedFastHaarWaveletTransformForNumIters(signal, num_iters);
        }
        else {
            throw new IllegalArgumentException("Illegal dwt value: " + dwt);
        }
    }
	
	public static void multiresSignalReconstruct(double[] signal_transform, double[] reconstructed_signal, DWT dwt, COEFF coeff, int num_scales, int scale_num) {
        final int signal_size = signal_transform.length;
        int range_start = 0;
        int range_end = 0;
        
        switch ( coeff ) {
            case S: 
                range_end = getSCoeffsEnd(signal_size, scale_num);
                break;
            case D:
                range_start = getDCoeffsStart(signal_size, scale_num);
                range_end = getDCoeffsEnd(signal_size, scale_num);
                break;
            default:
                throw new IllegalArgumentException("multiresSignalReconstruct: unknown coeff type");
        }
       
        System.arraycopy(signal_transform, 0, reconstructed_signal, 0, signal_size);
        for(int i = 0; i < signal_size; i++) {
            if ( i < range_start || i > range_end ) {
                reconstructed_signal[i] = 0;
            }
        }
        
        inverseDWTForNumIters(reconstructed_signal, dwt, num_scales);
    }
	
	static double[] generate_signal_with_noise_ex_4_4_p34(int len, int x, double v) {
        double[] signal = new double[len];
        double[] signal_domain = partition(0.0, 1.0, 1.0/len);
 
        Ripples_F_ex_4_4_p34 sf = new Ripples_F_ex_4_4_p34();
        
        for(int i = 0; i < len; i++) {
            signal[i] = sf.v(signal_domain[i]);
        }
        
        // add noise at x
        signal[x] = v;
        signal_domain = null;
        
        return signal;
    }
    
	static void display_signal_with_noise_ex_4_4_p34_512() {
        double[] signal = generate_signal_with_noise_ex_4_4_p34(512, 200, 4.0);
        displaySignal(signal);
        signal = null;
    }
	
	public static int getDCoeffsStart(int signal_len, int curr_scale) {
        return (int)(signal_len/Math.pow(2.0, curr_scale));
    }
	
	public static int getDCoeffsEnd(int signal_len, int curr_scale) {
        return 2*getDCoeffsStart(signal_len, curr_scale) - 1;
    }
	
	public static int getSCoeffsEnd(int signal_len, int curr_scale) {
        return (int)(signal_len/Math.pow(2.0, curr_scale))-1;
    }
	
	static void displaySignal(double[] signal) {
        for(double d: signal) {
            System.out.println(d);
        }
    }
}