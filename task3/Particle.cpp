//
// Created by Hesham Othman on 21.10.19.
//

#include <cmath>
#include "Particle.h"

cParticle::cParticle() {

}
cParticle::cParticle(const cParticle &particle) {
    for(int i = 0; i < array_size; ++i)
    {
        this->position[i] = particle.position[i];
        this->direction[i] = particle.direction[i];
    }
}


void cParticle::getPosition(double &x, double &y, double &z) const {
    x = this->position[0];
    y = this->position[1];
    z = this->position[2];
}

void cParticle::getDirection(double &x, double &y, double &z) const {
    x = this->direction[0];
    y = this->direction[1];
    z = this->direction[2];
}

void cParticle::setPosition(double x, double y, double z) {
    setElements(this->position, NULL, x,y,z);
}

void cParticle::setDirection(double x, double y, double z) {
    setElements(this->direction,NULL,x,y,z);
}

void cParticle::getPosition(double *p) const {
    setElements(p, this->position);

}
void cParticle::setDirection(double *d) {
    setElements(this->direction, d);
}

void cParticle::getDirection(double *d) const {
    setElements(d, this->direction);
}

void cParticle::setPosition(double *p) {
    setElements(this->position, p);
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

double cParticle::minDistToOrigin() const {
    double planeEquation[array_size] = {1,1,1};
    getDirection(planeEquation);
    // for origin planes, use 0 for planes passing by the origin
    double d = 6  ;
    double lambda = (d - dot(planeEquation, this->position)) /  dot(planeEquation, this->direction);
    double intersectionPoint[array_size] = {};
    for(int i = 0; i < array_size; ++i) {
        intersectionPoint[i] = lambda*this->direction[i] + this->position[i];
    }

    return std::sqrt(std::pow(intersectionPoint[0],2) +
                     std::pow(intersectionPoint[1],2) +
                     std::pow(intersectionPoint[2],2));
}






