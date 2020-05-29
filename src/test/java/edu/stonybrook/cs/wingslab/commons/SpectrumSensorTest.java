package edu.stonybrook.cs.wingslab.commons;

import org.junit.Assert;
import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

class SpectrumSensorTest {

    @Test
    void testToString() {
        SpectrumSensor ss1 = new SpectrumSensor( new RX(new Element(new Point(1, 1),
                10)), 1.0, 1);
        System.out.println(ss1);
    }

    @Test
    void sensorGeneratorRandom() {
        SpectrumSensor.sensorGeneratorRandom(400, new Square(1000), 0.732, 1.0, 15);
    }

    @Test
    void sensorGeneratorUniform() {
        SpectrumSensor.sensorGeneratorUniform(930, new Square(1000), 0.732, 1.0, 15);
    }

    @Test
    public void testSensorReader() throws IOException {
        SpectrumSensor[] sensors = SpectrumSensor.SensorReader("resources/sensors/square1000/900/sensors.txt");
        Assert.assertTrue(sensors.length == 900);
    }
}