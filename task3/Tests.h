//
// Created by Hesham Othman on 21.10.19.
//

#ifndef XRAYS_SIMULATIONS_TESTS_H
#define XRAYS_SIMULATIONS_TESTS_H

#include "Particle.h"

void testParticle();
void testMinDestToOrigin();

void testParticle() {
    const unsigned n = 100;	// number of tests per method

    // initialize random generator to perform different tests for each run
    srand((unsigned)time(NULL));

    // print heading
    cout << "+----------------------+" << endl;
    cout << "| test class cParticle |" << endl;
    cout << "+----------------------+" << endl;

    //========================================================
    cout << "test setPosition() and getPosition()" << endl;
    for (unsigned k = 0; k < n; k++) {
        double pos[3];		// random position of particle
        cParticle a, b;		// objects to perform tests on
        double p[3];		// porsition read from class
        double x, y, z;		// porsition read from class

        // random position
        pos[0] = rand()*1.0 / RAND_MAX;
        pos[1] = rand()*1.0 / RAND_MAX;
        pos[2] = rand()*1.0 / RAND_MAX;

        // perform tests
        a.setPosition(pos[0], pos[1], pos[2]);
        b.setPosition(pos);
        a.getPosition(p);
        a.getPosition(x, y, z);

        // check results
        if (a.p[0] != pos[0] || a.p[1] != pos[1] || a.p[2] != pos[2])
            throw runtime_error("Bad method setPosition(x, y, z).");
        if (b.p[0] != pos[0] || b.p[1] != pos[1] || b.p[2] != pos[2])
            throw runtime_error("Bad method setPosition(pos[3]).");
        if (x != pos[0] || y != pos[1] || z != pos[2])
            throw runtime_error("Bad method getPosition(&x, &y, &z).");
        if (p[0] != pos[0] || p[1] != pos[1] || p[2] != pos[2])
            throw runtime_error("Bad method getPosition(pos[3]).");
    }
    cout << "  finished" << endl;

    //========================================================
    cout << "test setDirection() and getDirection()" << endl;
    for (unsigned k = 0; k < n; k++) {
        double dir[3];		// random direction of particle
        cParticle a, b;		// objects to perform tests on
        double d[3];		// direction read from class
        double x, y, z;		// direction read from class

        // random direction
        dir[0] = rand()*1.0 / RAND_MAX;
        dir[1] = rand()*1.0 / RAND_MAX;
        dir[2] = rand()*1.0 / RAND_MAX;

        // perform tests
        a.setDirection(dir[0], dir[1], dir[2]);
        b.setDirection(dir);
        a.getDirection(d);
        a.getDirection(x, y, z);

        // check results
        if (a.u[0] != dir[0] || a.u[1] != dir[1] || a.u[2] != dir[2])
            throw runtime_error("Bad method setDirection(x, y, z).");
        if (b.u[0] != dir[0] || b.u[1] != dir[1] || b.u[2] != dir[2])
            throw runtime_error("Bad method setDirection(dir[3]).");
        if (x != dir[0] || y != dir[1] || z != dir[2])
            throw runtime_error("Bad method getDirection(&x, &y, &z).");
        if (d[0] != dir[0] || d[1] != dir[1] || d[2] != dir[2])
            throw runtime_error("Bad method getDirection(dir[3]).");
    }
    cout << "  finished" << endl;

    //========================================================
    cout << "test copy constructor" << endl;
    for (unsigned k = 0; k < n; k++) {
        double pos[3];	// random position of particle
        double dir[3];	// random direction of particle
        cParticle a;	// object to perform tests on

        // random position and direction
        pos[0] = rand()*1.0 / RAND_MAX;
        pos[1] = rand()*1.0 / RAND_MAX;
        pos[2] = rand()*1.0 / RAND_MAX;
        dir[0] = rand()*1.0 / RAND_MAX;
        dir[1] = rand()*1.0 / RAND_MAX;
        dir[2] = rand()*1.0 / RAND_MAX;
        a.setPosition(pos);
        a.setDirection(dir);

        // perform test
        cParticle b(a);

        // check result
        if (b.p[0] != pos[0] || b.p[1] != pos[1] || b.p[2] != pos[2])
            throw runtime_error("Bad position.");
        if (b.u[0] != dir[0] || b.u[1] != dir[1] || b.u[2] != dir[2])
            throw runtime_error("Bad direction.");
    }
    cout << "  finished" << endl;

    //========================================================
    cout << "test minDist2Origin()" << endl;
    for (unsigned k = 0; k < n; k++) {
        double pos[3];	// random position of particle
        double dir[3];	// random direction of particle
        double polar;	// random polar direction
        double azimuth;	// random azimuth
        // direction
        cParticle a;	// object to perform tests on
        double x, y, z;	// temporary values for distance calculation
        double d;		// calculated distance

        // random position and direction
        pos[0] = rand()*1.0 / RAND_MAX;
        pos[1] = rand()*1.0 / RAND_MAX;
        pos[2] = rand()*1.0 / RAND_MAX;
        polar = acos(1-rand()*2.0/RAND_MAX);
        azimuth = rand()*2.0*M_PI/RAND_MAX;
        dir[0] = sin(azimuth)*sin(polar);
        dir[1] = cos(azimuth)*sin(polar);
        dir[2] = cos(polar);

        // set parameters and evaluate distance
        a.setPosition(pos);
        a.setDirection(dir);
        d = dir[0] * pos[0] + dir[1] * pos[1] + dir[2] * pos[2];
        //d /= dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
        x = pos[0] - dir[0] * d;
        y = pos[1] - dir[1] * d;
        z = pos[2] - dir[2] * d;
        d = sqrt(x*x + y * y + z * z);

        // perform test and check result
        if (fabs(a.minDistToOrigin() - d) > 1e-10)
            throw runtime_error("Bad distance to origin.");
    }
    cout << "  finished" << endl;

    //========================================================
    cout << "test minDist2Point()" << endl;
    for (unsigned k = 0; k < n; k++) {
        double pos[3];	// random position of particle
        double dir[3];	// random direction of particle
        double polar;	// random polar direction
        double azimuth;	// random azumuthal direction
        double point[3];	// random point to evaluate minimum distance to
        cParticle a;	// object to perform tests on
        double x, y, z;	// temporary values for distance calculation
        double d;		// calculated distance

        // random position, direction and point
        pos[0] = rand()*1.0 / RAND_MAX;
        pos[1] = rand()*1.0 / RAND_MAX;
        pos[2] = rand()*1.0 / RAND_MAX;
        polar = acos(1 - rand()*2.0 / RAND_MAX);
        azimuth = rand()*2.0*M_PI / RAND_MAX;
        dir[0] = sin(azimuth)*sin(polar);
        dir[1] = cos(azimuth)*sin(polar);
        dir[2] = cos(polar);
        point[0] = rand()*1.0 / RAND_MAX;
        point[1] = rand()*1.0 / RAND_MAX;
        point[2] = rand()*1.0 / RAND_MAX;

        // set parameters and evaluate distance
        a.setPosition(pos);
        a.setDirection(dir);
        d = dir[0] * (pos[0]-point[0]) + dir[1] * (pos[1]-point[1]) + dir[2] * (pos[2]-point[2]);
        x = pos[0] - point[0] - dir[0] * d;
        y = pos[1] - point[1] - dir[1] * d;
        z = pos[2] - point[2] - dir[2] * d;
        d = sqrt(x*x + y*y + z*z);

        // perform tests and check results
//        if (fabs(a.minDist2Point(point) - d) > 1e-10)
//            throw runtime_error("Bad distance to origin.");
//        if (fabs(a.minDist2Point(point[0], point[1], point[2]) - d) > 1e-10)
//            throw runtime_error("Bad distance to origin.");
    }
    cout << "  finished" << endl;
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

