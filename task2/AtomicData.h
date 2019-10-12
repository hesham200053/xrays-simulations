//
// Created by Hesham Othman on 07.10.19.
// class to read and store the cross section data to provide it for any photon energy by
// log-log interpolation.
//

#ifndef SPECTRUM_ATOMICDATA_H
#define SPECTRUM_ATOMICDATA_H

#include <vector>

typedef struct {
    float x;
    float y;
} tPoint2d;

class cAtomicData {
private:
    // should the type def put here ?
    static std::vector<tPoint2d> tcs[100];
    static float A[100];
    static bool prepared;
public:
    cAtomicData();
//    ~cAtomicData();
    void prepare();
    double getStdAtomicWeight(unsigned Z);
    double getTotalCrossSection(unsigned Z, double energy);
};

#endif //SPECTRUM_ATOMICDATA_H
