//
// Created by Hesham Othman on 17.11.19.
//

//#ifndef XRAYS_SIMULATIONS_TESTS_H
//#define XRAYS_SIMULATIONS_TESTS_H

#include "MedImage.h"
#include "MedImageIO.h"
#include "Simulation.h"

void testMedImage();
void testSimulation();
void testMedImage() {
    cMedImage<double> image;
    unsigned size = 99;

    // create image initialized with zero
    image.create(size,size,0);

    // draw image
    for (unsigned i = 0; i < size; i++) {
        image[i][i] = image[i][size - i - 1] = 1;
    }
        // set pixel value for black and white
    image.pixelValueForBlack = -1;
    image.pixelValueForWhite = 2 ;

    // store the image as a bitmap file
    cMedImageIO<double>::writeBmp(image,"deleteMe.bmp");

}
void testSimulation() {
    cMedImage<double> image;
    image.create(50,50, 0);
//    image.nCol = 50;
//    image.nRow = 50;
    cSimulation simulation;
    simulation.setTubeVoltage(40);
    simulation.setCurrentTimeProduct(1);
    simulation.setSphereRadious(0.0025);
    simulation.setPixelSize(0.00015);
    simulation.prepare();
    simulation.simulate(image);
    image.autoWindowing();
//    cMedImageIO<double>::writeBmp(image,"first.bmp");
    cMedImageIO<double>::writeBmp(image,"second.bmp");

}
//#endif //XRAYS_SIMULATIONS_TESTS_H
