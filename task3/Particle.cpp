//
// Created by Hesham Othman on 21.10.19.
//

#include <cmath>
#include "Particle.h"
#include <iostream>

cParticle::cParticle() {

}
cParticle::cParticle(const cParticle &particle) {
    for(int i = 0; i < array_size; ++i)
    {
        this->p[i] = particle.p[i];
        this->u[i] = particle.u[i];
    }
}


void cParticle::getPosition(double &x, double &y, double &z) const {
    x = this->p[0];
    y = this->p[1];
    z = this->p[2];
}

void cParticle::getDirection(double &x, double &y, double &z) const {
    x = this->u[0];
    y = this->u[1];
    z = this->u[2];
}

void cParticle::setPosition(double x, double y, double z) {
    setElements(this->p, NULL, x,y,z);
}

void cParticle::setDirection(double x, double y, double z) {
    setElements(this->u,NULL,x,y,z);
}

void cParticle::getPosition(double *p) const {
    setElements(p, this->p);

}
void cParticle::setDirection(double *d) {
    setElements(this->u, d);
}

void cParticle::getDirection(double *d) const {
    setElements(d, this->u);
}

void cParticle::setPosition(double *p) {
    setElements(this->p, p);
}

void cParticle::setElements(double *src, const double *dest, double x, double y, double z) const {
    if (dest != NULL) {
        for(int i = 0; i < array_size ; ++i) {
            src[i] = dest[i];
        }
    } else {
        src[0] = x;
        src[1] = y;
        src[2] = z;
    }

}
double cParticle::sum(const double *arr) const {
    double result = 0;
    for(int i = 0; i< array_size ; i++){
        result+=arr[i];
    }
    return result;
}

double cParticle::dot(const double *arr1, const double *arr2) const {
    double result = 0;
    for(int i = 0; i< array_size ; i++){
        result+= arr1[i]*arr2[i];
    }
    return result;
}

double cParticle::magnitude(const double *point) const {
    return std::sqrt(point[0] * point[0] +
                     point[1] * point[1] +
                     point[2] * point[2] );
};
double cParticle::minDistToOrigin() const {
//    you may need plane equation in case of a general use case
//    double planeEquation[array_size] = {};
//    getDirection(planeEquation);
    // for origin planes, use 0 for planes passing by the origin
    double d = 0  ;
    double lambda = (d - dot(u, p)) /  dot(u, u);
    double intersectionPoint[array_size] = {};
    for(int i = 0; i < array_size; ++i) {
        intersectionPoint[i] = lambda * u[i] + p[i];
    }
    return magnitude(intersectionPoint);
}

double cParticle::minDist2Point(double x, double y, double z) const {
    double temp[array_size] = {x,y,z};
    return minDist2Point(temp);
}

//CD = |ABxAC| / |AB|
double cParticle::minDist2Point(double *point) const {
//    double d = u[0] * (p[0]-point[0]) + u[1] * (p[1]-point[1]) + u[2] * (p[2]-point[2]) ;

//    double lambda = (d - dot(u, p)) /  dot(u, u);
//    double intersectionPoint[array_size] = {};
//    for(int i = 0; i < array_size; ++i) {
//        intersectionPoint[i] = lambda * u[i] + p[i];
//    }
//    double distance[array_size] = {intersectionPoint[0] - point[0],
//                                   intersectionPoint[1] - point[1],
//                                   intersectionPoint[2] - point[2]};
//    return magnitude(distance);

    double d = u[0] * (p[0]-point[0]) + u[1] * (p[1]-point[1]) + u[2] * (p[2]-point[2]);
    double x = p[0] - point[0] - u[0] * d;
    double y = p[1] - point[1] - u[1] * d;
    double z = p[2] - point[2] - u[2] * d;
    d = sqrt(x*x + y*y + z*z);
    return d;

}






