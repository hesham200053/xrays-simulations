//
// Created by Hesham Othman on 07.10.19.
//

#ifndef SPECTRUM_MATERIAL_H
#define SPECTRUM_MATERIAL_H


#include <vector>
#include <string>
#include <cmath>

static const double Na = 6.022140857e23;

class cMaterial {
private:
    std::vector<unsigned> Z;            // each element std atomic number
    std::vector<unsigned> fraction;     // no. particles for each element
    double density;                     // in g/m^3
    std::string name;
public:
     cMaterial();
     cMaterial(const cMaterial &newMaterial);
//     ~cMaterial();
     void clear();
     void setName(const std::string &newName);
     const std::string &getName();
     void addElement(unsigned newZ, double newFraction);
     unsigned getNoOfElements();
     unsigned getAtomicNumber(unsigned i);
     double getFraction(unsigned i);
     void setDensity(double newDensity);
     const double &getDensity();
     double  getAttCoeff(double energy);
     void getAttSpec( std::vector<double> &spec,
                      double minEnergy, double tubeVoltage,
                      unsigned energySteps);
     double getMeanFreePath(double energy);

};


#endif //SPECTRUM_MATERIAL_H
