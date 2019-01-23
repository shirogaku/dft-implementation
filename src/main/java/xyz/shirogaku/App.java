package xyz.shirogaku;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

/**
 * Hello world!
 *
 */

/**
 * DFT
 * Euler's Identity = e^(iθ) = cos(θ) + isin(θ)
 * 
 *       N-1
 * X(k) = Σ  x(n)e^-(ikn2π/N)
 *       n=0
 * 
 * IDFT
 *           N-1
 * x(n) = 1/N Σ  X(k)e^(ikn2π/N)
 *           k=0
 * 
 * N = Sample data length
 * k = Frequency-axis sample position
 * n = Time-axis sample position
 * i = Imaginary number
 * 
 * Thanks to https://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
 */

public class App {
    public static void dft(double[] inreal, double[] inimag,
                            double[] outreal, double[] outimag) {
        if(!(inreal.length == inimag.length && inimag.length == outreal.length && outreal.length == outimag.length))
            throw new IllegalArgumentException("Real and imaginary length are not same.");

        int N = inreal.length;
        for(int k = 0; k < N; k++) {                                                                        // Loop through all frequency bin
            double sumreal = 0, sumimag = 0;
            for(int n = 0; n < N; n++) {                                                                    // Loop through all sample data and sum all value
                double angle = 2 * Math.PI * n * k / N;
                sumreal +=  inreal[n] * Math.cos(angle) + inimag[n] * Math.sin(angle);                      // Re(x(n)) * cos(2π(nk/N)) + Im(x(n)) * sin(2π(nk/N))
                                                                                                            //                          +
                sumimag += -inreal[n] * Math.sin(angle) + inimag[n] * Math.cos(angle);                      //-Re(x(n)) * sin(2π(nk/N)) + Im(x(n)) * cos(2π(nk/N))
            }
            outreal[k] = sumreal;
            outimag[k] = sumimag;
        }
    }

    public static void idft(double[] inreal, double[] inimag, double[] value) {
        int N = inreal.length;
        for(int n = 0; n < N; n++) {
            double sum = 0;
            for(int k = 0; k < N; k++) {
                double angle = 2 * Math.PI * n * k / N;
                sum += (inreal[k] * Math.cos(angle)) - (inimag[k] * Math.sin(angle));
            }
            value[n] = sum / N;
        }
    }

    public static double getFreqPerBin(double samplingrate, long N) {
        return samplingrate / N;                                                                               // Bin * Samplingrate / N
    }

    public static double getVectorValue(double real, double imag) {
        return Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));                                               //√(Re^2) + (Im^2)
    }

    public static double getPhaseOfBin(double real, double imag) {
        return Math.atan2(imag, real);
    }

    public static void main( String[] args ) throws IOException {
        SinData data = new SinData(1024, 0.1);
        double[] iimag = new double[data.getValue().length];
        double[] oreal = new double[data.getValue().length];
        double[] oimag = new double[data.getValue().length];
        double[] vector = new double[data.getValue().length];
        double[] phase = new double[data.getValue().length];
        double[] idft = new double[data.getValue().length];

        for(int i = 0;i < iimag.length; i++) {
            iimag[i] = 0;
            idft[i] = 0;
        }
        dft(data.getValue(), iimag, oreal, oimag);
        
        for(int i = 0; i < vector.length; i++) {
           vector[i] = getVectorValue(oreal[i], oimag[i]);
           phase[i] = getPhaseOfBin(oreal[i], oimag[i]);
        }
        idft(oreal, oimag, idft);
        
        try(BufferedWriter bw = Files.newBufferedWriter(Paths.get("./dft.csv"), StandardCharsets.UTF_8, StandardOpenOption.CREATE)) {
            bw.write("original,oreal,oimag,vector,phase,idft");
            bw.newLine();
            for(int i = 0; i < oreal.length; i++) {
                bw.write(String.format("%f,%f,%f,%f,%f,%f",data.getValue()[i] , oreal[i], oimag[i], vector[i], phase[i], idft[i]));
                bw.newLine();
            }
        }
    }
}
