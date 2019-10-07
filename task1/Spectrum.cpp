//
// Created by Hesham Othman on 02.10.19.
//

using namespace std;
#include <cmath>
#include <iostream>
#include "Spectrum.h"
#include "fstream"

cSpectrum::cSpectrum() {
    *this = 0;
}
cSpectrum::cSpectrum(const cSpectrum &newSpec) {
    spec = newSpec.spec;
}
unsigned cSpectrum::size() {
    return (unsigned) spec.size();
}
void cSpectrum::resize(unsigned newSize) {
    spec.resize(newSize);
}

const vector<double> cSpectrum::get() const {
    //return spec;
    return this->spec;
}
// what does the ampersand do here ?
double& cSpectrum:: operator[](unsigned index) {
    if (index >= spec.size())
        throw runtime_error("Index_out_of_range_in_class_cSpectrum");
    return spec[index];
}
cSpectrum &cSpectrum::operator=(double value) {
    for (unsigned i = 0; i<spec.size() ; i++)
        spec[i] = value ;
    // what is the asterisk here ?
    return *this ;
}
cSpectrum &cSpectrum::operator=(const vector<double > &newSpec) {
    spec = newSpec;
    return *this;
}
cSpectrum &cSpectrum::operator=(const cSpectrum &newSpec) {
    spec = newSpec.spec;
    return *this;
}
cSpectrum &cSpectrum::operator+(const cSpectrum &summand) {
    if (spec.size() != summand.spec.size())
        throw runtime_error("adding spectra of different sizes not possible");
    for (unsigned i = 0; i < spec.size() ; i++) {
        spec[i] += summand.spec[i];
    }
    return *this;
}
// multiplying vectors element wise
cSpectrum &cSpectrum::operator*(const cSpectrum &factor) {
    if (spec.size() != factor.spec.size())
        throw runtime_error("multiplying spectra of different sizes not possible");
    for (unsigned i = 0; i < spec.size() ; i++) {
        spec[i] *= factor.spec[i];
    }
    return *this;
}
// kind like a scaler for the whole vector ?
cSpectrum &cSpectrum::operator*(double factor) {
    for (unsigned i = 0; i < spec.size() ; i++) {
        spec[i] *= factor;
    }
    return *this;
}
// ??
cSpectrum &cSpectrum::exp() {
    for (unsigned i = 0; i < spec.size() ; i++) {
        spec[i] = ::exp(spec[i]);
    }
    return *this;
}

double cSpectrum::sum() {
    double result;
    for (unsigned i = 0; i < spec.size() ; i++) {
        result += spec[i];
    }
    return result;
}

// when you pass values by ambersand, it makes you able to control the variables and can access them
void cSpectrum::readSpectrum(const string &fname, double tubeVoltage,
                             double &minEnergy, string &spectrumName) {
    ifstream inp;
    bool valueFound = false;
    // open file
    inp.open(fname, ios::binary);
    if(!inp.is_open()) throw runtime_error(" file couldn't be opened");

    // 1- title
    char title[15];
    inp.read(title, 14);
    title[14] = 0;
    // cout << "title: " << title << endl;

    // 2- length of the name Nn
    uint32_t Nn;
    inp.read((char*) &Nn, sizeof(Nn));
    // cout << "Nn: " << Nn << endl;

    // 3- spectrum name
    char *charOfName = new char[Nn+1];
    inp.read(charOfName, Nn);
    charOfName[Nn] = 0;
    spectrumName = charOfName;
    // cout << "Name characters: " << charOfName << endl;

    // 4- spectra number
    uint32_t Ns;
    inp.read((char*) &Ns, sizeof(Ns));
    // cout << "Ns: " << Ns << endl;

    for (int i = 0 ; i < Ns ; i++) {
        float tVoltage;
        inp.read((char*) &tVoltage, sizeof(tVoltage));

        // minimum Energy
        float mEnergy;
        inp.read((char*) &mEnergy, sizeof(mEnergy));
        // cout << tVoltage << ", " << mEnergy << endl;
        minEnergy = mEnergy;

        // table length
        uint32_t Nt;
        inp.read((char*) &Nt, sizeof(Nt));

        float tableEntry;
        if (fabs(tVoltage - tubeVoltage) < 1e-3) {
            valueFound = true;
            cout << spectrumName << " spectra values for minimum energy of " << minEnergy <<
            " and a tube voltage of " << tubeVoltage << endl;
            for (int j = 0; j < Nt; j++) {
                inp.read((char*) &tableEntry, sizeof(tableEntry));
                // for debugging reasons 
                // cout << j << ": " << tableEntry << endl;
            }
        }  else {
            //skip this entry
            inp.seekg(Nt * sizeof(tableEntry), ios_base::cur);
        }
    }

    if (!valueFound) throw runtime_error("value_couldn't_be_found");
}
