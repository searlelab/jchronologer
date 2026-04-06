package org.searlelab.jchronologer.dlib;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.InflaterOutputStream;

/**
 * Encoding helpers for SQLite DLIB blobs.
 */
public final class DlibCodec {

    private DlibCodec() {
    }

    public static byte[] encode(double[] values) {
        ByteBuffer buffer = ByteBuffer.allocate(values.length * Double.BYTES).order(ByteOrder.LITTLE_ENDIAN);
        for (double value : values) {
            buffer.putDouble(value);
        }
        return buffer.array();
    }

    public static byte[] encode(float[] values) {
        ByteBuffer buffer = ByteBuffer.allocate(values.length * Float.BYTES).order(ByteOrder.LITTLE_ENDIAN);
        for (float value : values) {
            buffer.putFloat(value);
        }
        return buffer.array();
    }

    public static byte[] encode(boolean[] values) {
        byte[] encoded = new byte[values.length];
        for (int i = 0; i < values.length; i++) {
            encoded[i] = values[i] ? (byte) 1 : (byte) 0;
        }
        return encoded;
    }

    public static boolean[] unitBooleanArray(int size) {
        boolean[] values = new boolean[size];
        for (int i = 0; i < size; i++) {
            values[i] = true;
        }
        return values;
    }

    public static byte[] compress(byte[] bytes) {
        try {
            ByteArrayOutputStream buffer = new ByteArrayOutputStream();
            try (DeflaterOutputStream output = new DeflaterOutputStream(buffer)) {
                output.write(bytes);
            }
            return buffer.toByteArray();
        } catch (IOException e) {
            throw new IllegalStateException("Failed to compress DLIB payload.", e);
        }
    }

    public static byte[] decompress(byte[] bytes) {
        try {
            ByteArrayOutputStream buffer = new ByteArrayOutputStream();
            try (InflaterOutputStream output = new InflaterOutputStream(buffer)) {
                output.write(bytes);
            }
            return buffer.toByteArray();
        } catch (IOException e) {
            throw new IllegalStateException("Failed to decompress DLIB payload.", e);
        }
    }

    public static double[] decodeDoubles(byte[] bytes) {
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        double[] values = new double[bytes.length / Double.BYTES];
        for (int i = 0; i < values.length; i++) {
            values[i] = buffer.getDouble();
        }
        return values;
    }

    public static float[] decodeFloats(byte[] bytes) {
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[bytes.length / Float.BYTES];
        for (int i = 0; i < values.length; i++) {
            values[i] = buffer.getFloat();
        }
        return values;
    }
}
