package edu.stonybrook.cs.wingslab.commons;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import smile.interpolation.KrigingInterpolation;
import smile.interpolation.variogram.ExponentialVariogram;

import java.util.Comparator;
import java.util.PriorityQueue;

/**2-D Ordinary Kriging interpolation method.
 * A list of locations would be passed along with values.
 * @author Mohammad Ghaderi, <mghaderibane@cs.stonybrook.edu>
 * @since 1.1
 */
public class OrdinaryKriging {
//    private class SensorDistance {// used to find nearest objects(heap)
//        int id;
//        double distance;
//
//        SensorDistance(int id, double distance, double power) {
//            this.id = id;
//            this.distance = distance;
//        }
//
//        double getDistance() { return this.distance; }
//
//        int getId() { return this.id; }
//    }
//
//    private final boolean detrended;
//    private final boolean selectPartitionCenter;
//    // array of sensors
//    private final SpectrumSensor[] sss;
//    // array of PUs + Active SUs
//    private final TX[] txs;
//    // number of close sensors to be selected
//    private final int numSssSelected;
//    // OK PATH-LOSS exponent
//    private final double alpha;
//    private final LogDistancePM logDistancePM;
//    private double[] txsSuPL;
//    private final int cellSize;
//    private final Element requestingSu;


//    public double[] getTxsSuPL() { return txsSuPL; }
    private final double[][] locations;
    private final double[] values;
    private final double x;
    private final double y;

    /**
     * OrdinaryKriging constructor given a set of locations and values to be fitted and a location to do interpolation
     * @param locations set of locations <x,y>
     * @param values values to be trained
     * @param x x coordinate of interpolating location
     * @param y y coordinate of interpolating location*/
    public OrdinaryKriging(double[][] locations, double[] values, double x, double y){
        super();
        if (locations.length != values.length)
            throw new IllegalArgumentException("locations and values do not have the same length!");
        this.locations = locations;
        this.values = values;
        this.x = x;
        this.y = y;
    }

    public double interpolate(){
        int n = this.locations.length;
        ParametricUnivariateFunction variogramFunc = variogramFunction("Exponential");
        CurveFitter<ParametricUnivariateFunction> curveFitter = new CurveFitter<>(new LevenbergMarquardtOptimizer());
        for (int i = 0; i < n; i++)
            curveFitter.addObservedPoint(distance(locations[i][0], locations[i][1], this.x, this.y), values[i]);
        double[] initialGuess = {Math.random(), Math.random()};
        System.out.println("fitting");
        double[] params = curveFitter.fit(variogramFunc, initialGuess); //params[0]: sill, params[1]: range

        if (params[0] == 0)
            return 0.0;
        System.out.println("interpolating");
        KrigingInterpolation krigingInterpolation = new KrigingInterpolation(this.locations, this.values,
                new ExponentialVariogram(Math.abs(params[1]), Math.abs(params[0])), null);

        return krigingInterpolation.interpolate(this.x, this.y);
    }

    private static double distance(double x1, double y1, double x2, double y2){
        return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
    }

    private static ParametricUnivariateFunction variogramFunction(String type){
        if (type.equals("Exponential"))
            return new ParametricUnivariateFunction() {
                @Override
                public double value(double v, double... doubles) {
                    //doubles[0]: sill, doubles[1]: range
                    return doubles[0] * (1 - Math.exp(-v / doubles[1])); // sill * (1-exp(-v/range))
                }

                @Override
                public double[] gradient(double v, double... doubles) {
                    return new double[]{
                            1.0 * (1 - Math.exp(-v / doubles[1])), //df/d(sill) = (1 - exp(-v / range))
                    -doubles[0] * Math.exp(-v / doubles[1]) * v / (doubles[1] * doubles[1])};
                    //df/d(range) = - sill * exp(-v / range) *  v / range^2 = - v * sill * exp(-v/range) / range^2
                }
            };
        return null;
    }
//    /**OrdinaryKriging constructor
//     * @param sss list of spectrum sensors
//     * @param txs list of transmitters
//     * @param numSssSelected number of selected sensors to do interpolation
//     * @param detrended is detrended enabled
//     * @param selectPartitionCenter is center of partitions going to be selected
//     * @param alpha path-loss exponent
//     * @param cellSize cell size of the grid
//     * */
//    public OrdinaryKriging(TX[] txs, SpectrumSensor[] sss, Element requestingSu,
//                           double alpha, int cellSize, int numSssSelected,
//                           boolean detrended, boolean selectPartitionCenter){
//        super();
//        this.sss = sss;
//        this.txs = txs;
//        this.alpha = alpha;
//        this.numSssSelected = numSssSelected;
//        this.detrended = detrended;
//        this.selectPartitionCenter = selectPartitionCenter;
//        this.logDistancePM = new LogDistancePM(this.alpha);
//        this.cellSize = cellSize;
//        this.requestingSu = requestingSu;
//        this.txsSuPL = new double[txs.length];
//
//        // interpolate for each tx
//        for (int txIdx = 0; txIdx < txs.length; txIdx++){
//
//        }
//    }
//
//    /**OrdinaryKriging constructor with default values
//     * @param sss list of spectrum sensors
//     * @param txs list of transmitters
//     * @param alpha path-loss exponent
//     * @param cellSize cell size of the grid
//     * */
//    public OrdinaryKriging(TX[] txs, SpectrumSensor[] sss, Element requestingSu, double alpha, int cellSize){
//        this(txs, sss, requestingSu, alpha, cellSize, 5, true, true);
//    }
//
//    // distance component based on log-distance and alpha
//    private double distanceComponent(double distance){
//        return this.logDistancePM.pathLoss(distance);
//    }
//
//    // find num nearest objects using heap
//    private static SensorDistance[] findNearestElements(SensorDistance[] elements, int num){
//        PriorityQueue<SensorDistance> pq = new PriorityQueue<>(num,
//                Comparator.comparing(SensorDistance::getDistance).reversed());
//        for (SensorDistance element : elements) {
//            if (pq.size() < num)
//                pq.add(element);
//            else if (element.distance < pq.peek().distance) {
//                pq.poll();
//                pq.add(element);
//            }
//        }
//        SensorDistance[] nearest = new SensorDistance[pq.size()];
//        int cnt = 0;
//        for (SensorDistance elementDistance : pq)
//            nearest[cnt++] = elementDistance;
//        return nearest;
//    }
}
