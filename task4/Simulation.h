//
// Created by Hesham Othman on 17.11.19.
//

#ifndef XRAYS_SIMULATIONS_SIMULATION_H
#define XRAYS_SIMULATIONS_SIMULATION_H


#include <Spectrum.h>
#include "../task3/Particle.h"
#include "MedImage.h"

class cSimulation {
private:
    // all positions in meter and voltages in Kv
    double tubePos[array_size];
    double tubeVoltage;
    double spherePos[array_size];
    double sphereRadious;
    double pixelSize;
    cSpectrum xRayTube;
    cSpectrum attCoeff;
    cSpectrum energy;
    double current_time_product;
public:
    void setCurrentTimeProduct(double currentTimeProduct);

private:
    void getPixelDirection(double x, double y, double pDouble[3]);
public:
    cSimulation();
    ~cSimulation();
    void setTubePosition(double *p);
    void setSpherePosition(double *p);
    void setTubeVoltage(double tubeVoltageKv);
    void setSphereRadious(double sphereRadious);
    void setPixelSize(double pixelSize);
    void prepare();
    void simulate(cMedImage<double> &image);

    double getT(double center, double sphereRadious, double *t);
};


#endif //XRAYS_SIMULATIONS_SIMULATION_H
