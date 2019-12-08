//
// Created by Hesham Othman on 07.10.19.
//

#include "Material.h"
#include "AtomicData.h"
//#include "../task1/Spectrum.h"


cMaterial::cMaterial() {
    clear();
}

cMaterial::cMaterial(const cMaterial &newMaterial) {
    // what is the diff between those two
    // this->Z = newMaterial.Z;
    Z = newMaterial.Z;
    fraction = newMaterial.fraction;
    density = newMaterial.density;
    name = newMaterial.name;
}

void cMaterial::clear() {
    Z = {};
    fraction = {};
    density = 1;
    name = "new material";
}

void cMaterial::setName(const std::string &newName) {
    name = newName;
}

const std::string &cMaterial::getName() {
    return name;
}

void cMaterial::addElement(unsigned newZ, double newFraction) {
    Z.push_back(newZ);
    fraction.push_back(newFraction);
}

unsigned cMaterial::getNoOfElements() {
    return Z.size(); // or fraction.size() ?
}

unsigned cMaterial::getAtomicNumber(unsigned i) {
    return Z[i];
}

double cMaterial::getFraction(unsigned i) {
    return fraction[i];
}

void cMaterial::setDensity(double newDensity) {
    density = newDensity;
}

const double &cMaterial::getDensity() {
    return density;
}

// the units here are pretty important
double cMaterial::getAttCoeff(double energy) {
    cAtomicData atomicData;
    atomicData.prepare();
    double sum_sigK_fk, sum_ak_fk;
    int i = 0;
    for (auto &v: Z) {
        sum_sigK_fk += atomicData.getTotalCrossSection(v, energy) * fraction[i];
        sum_ak_fk += atomicData.getStdAtomicWeight(v) * fraction[i++];
    }
    return density * Na * (sum_sigK_fk/sum_ak_fk) * pow(10,-22) ;
}

// preparing the attenuation spectrum after a specific material
void cMaterial::getAttSpec(std::vector<double> &spec, double minEnergy, double tubeVoltage, unsigned energySteps) {
//    double step = ::ceil((tubeVoltage - minEnergy) / energySteps);
    double step = (tubeVoltage - minEnergy) / energySteps;
    double energy = minEnergy;
    while ( energy + step/2 < tubeVoltage ) {
        spec.push_back(getAttCoeff(energy));
        energy += step;
    }
}


double cMaterial::getMeanFreePath(double energy) {
    return 1 / getAttCoeff(energy);
}



