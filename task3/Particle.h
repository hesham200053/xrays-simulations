//
// Created by Hesham Othman on 21.10.19.
//

#ifndef XRAYS_SIMULATIONS_PARTICLE_H
#define XRAYS_SIMULATIONS_PARTICLE_H


#include <vector>

const int array_size = 3;
class cParticle {
friend void testParticle();
private:
    double p[array_size];
    double u[array_size];
    void setElements(double src[array_size], const double *dest, double x = 0, double y = 0, double z = 0) const;
    double magnitude(const double *point) const;
 public:
    cParticle();
    cParticle(const cParticle &particle);

    void getPosition(double p[array_size]) const;
    void getPosition(double &x, double &y, double &z) const;

    void getDirection(double d[array_size]) const;
    void getDirection(double &x, double &y, double &z) const;

    void setPosition(double p[array_size]);
    void setPosition(double x, double y, double z);

    void setDirection(double d[array_size]);
    void setDirection(double x, double y, double z);

    double minDistToOrigin() const;
    double minDist2Point(double x, double y, double z) const;
    double minDist2Point(double d[array_size]) const;

    double dot(const double *arr1, const double *arr2) const;
    double sum(const double *arr) const;
};


#endif //XRAYS_SIMULATIONS_PARTICLE_H
