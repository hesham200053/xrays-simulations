//
// Created by Hesham Othman on 21.10.19.
//

#ifndef XRAYS_SIMULATIONS_TESTS_H
#define XRAYS_SIMULATIONS_TESTS_H

#include "Particle.h"

void testParticle();
void testMinDestToOrigin();

void testParticle() {

}

void testMinDestToOrigin() {
    cParticle particle;
    double position[array_size] = {1,0,1};
    double direction[array_size] = {3,-2,1};
    particle.setPosition(position);
    particle.setDirection(direction);

    double restult = particle.minDistToOrigin();
    cout << restult << endl;
}


#endif //XRAYS_SIMULATIONS_TESTS_H

