//
// Created by Hesham Othman on 17.11.19.
//

#include <Material.h>
#include "Simulation.h"
#include "MedImage.h"
#include "MedImageIO.h"
#include "../task5/random.h"

cSimulation::cSimulation() {
    setTubePosition(new double[3] {0,0,1.2});
    setTubeVoltage(75);
    setSpherePosition(new double[3] {0,0,0.2});
    setSphereRadious(0.15);
    setPixelSize(.01);
    // in mAs
    setCurrentTimeProduct(10);

};

cSimulation::~cSimulation() {

}
void cSimulation::setTubeVoltage(double tubeVoltage) {
    cSimulation::tubeVoltage = tubeVoltage;
}
void cSimulation::setSphereRadious(double sphereRadious) {
    cSimulation::sphereRadious = sphereRadious;
}
void cSimulation::setPixelSize(double pixelSize) {
    cSimulation::pixelSize = pixelSize;
}

void cSimulation::setTubePosition(double p[array_size]) {
    tubePos[0] = p[0];
    tubePos[1] = p[1];
    tubePos[2] = p[2];
}

void cSimulation::setSpherePosition(double p[array_size]) {
    spherePos[0] = p[0];
    spherePos[1] = p[1];
    spherePos[2] = p[2];
}

void cSimulation::prepare() {
    double minEnergy;
    std::string spectrumName;

    xRayTube.readSpectrum("/Users/Hesham/CLionProjects/xrays-simulations/task1/SRO33100ROT350.dat", tubeVoltage, minEnergy , spectrumName);
    energy.resize(xRayTube.size());
    // multiply the spec by the corresponding energy
    // ambiguous
    double step = (tubeVoltage - minEnergy) / xRayTube.size();
    for(signed i = 0; i < xRayTube.size() ; i++ ) {
        energy[i] = (minEnergy + step*i);
//        xRayTube *= energy[i];
    }
    // the Z component of the tubePos is the distance between the tube and the detector. ????
    double f = ((pixelSize*pixelSize) / (tubePos[2] * tubePos[2])) * 1000000;
    f *= current_time_product;
    f *= step;
    xRayTube *= f;

    cMaterial water;
    std::vector<double> spec;
    water.addElement(1, 2);
    water.addElement(8, 1);
    water.setDensity(1.0);
    water.setName("water particle");
    water.getAttSpec(spec, minEnergy, tubeVoltage, xRayTube.size());
    attCoeff = spec;
    // attCoeff is 101  ????
//    attCoeff.resize(100);
}

void cSimulation::setCurrentTimeProduct(double currentTimeProduct) {
    current_time_product = currentTimeProduct;
}

void cSimulation::simulate(cMedImage<double> &image) {
    cParticle particle;
    double pixelDir[3];
//    double pixelPos[2] = {};
    double mindistToSphereCenter;
    double t = 0;
    double sum;
    double stdev;
//    double newSum;
    cRandKiss kiss((uint32_t)time(NULL));
    cRandNormal normal(kiss);
    normal.prepare();

    // draw image
    for (signed i = 0 ; i < image.nRow; i++) {
        for (signed j = 0 ; j < image.nCol; j++) {
            particle.setPosition(tubePos);
            getPixelDirection(pixelSize* (j - (image.nCol-1)/2.0), pixelSize* (i - (image.nRow-1)/2.0), pixelDir);
            particle.setDirection(pixelDir);
            mindistToSphereCenter = particle.minDist2Point(spherePos);
            if (mindistToSphereCenter < sphereRadious ) {
                t = getT(mindistToSphereCenter,sphereRadious,  &t);
            } else {
                t = 0;
            }
            cSpectrum I;
            I = attCoeff;
            I *= -t;
            I.exp();
            I *= xRayTube;

            I *= energy;
            sum = I.sum();
            stdev = sqrt(sum);
            sum += stdev*normal.rand();

            image[i][j] = ::log(sum);
        }
    }
    // store the image as a bitmap file
}

void cSimulation::getPixelDirection(double x, double y, double *pixelDir) {
    double length;
     pixelDir[0] = x-tubePos[0];
     pixelDir[1] = y-tubePos[1];
     pixelDir[2] = 0-tubePos[2];
     length = sqrt(pixelDir[0] *pixelDir[0] + pixelDir[1]*pixelDir[1] + pixelDir[2]*pixelDir[2]);
     pixelDir[0] /= length;
     pixelDir[1] /= length;
     pixelDir[2] /= length;
}

double cSimulation::getT(double distToC, double sphR, double *t) {
    return 2 * sqrt((sphR*sphR) - (distToC*distToC)) ;
}


