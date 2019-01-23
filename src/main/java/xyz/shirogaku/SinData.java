package xyz.shirogaku;

public class SinData {
    private int mSize = 1024;
    private double mIncrement = 0.1; 
    private double[] mValue;

    public SinData() {
        this.generateData();
    };

    public SinData(int size, double increment) {
        mSize = size;
        mIncrement = increment;
        this.generateData();
    }

    public double[] getValue() {
        return mValue;
    }

    private void generateData() {
        mValue = new double[mSize];
        double sin_value = 0.01;
        for(int i = 0; i < mSize; i++) {
            mValue[i] = Math.sin(sin_value);
            sin_value =+ mIncrement;
        }
    }
}