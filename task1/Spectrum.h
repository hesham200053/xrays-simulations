//
// Created by Hesham Othman on 02.10.19.
//

#ifndef SPECTRUM_SPECTRUM_H
#define SPECTRUM_SPECTRUM_H

#include <vector>

class cSpectrum {
private:
    std::vector<double> spec;
public:
    // copy constructor
    cSpectrum();
    cSpectrum(const cSpectrum &newSpec);
    unsigned size();
    void resize(unsigned newSize);
    const std::vector<double> get() const;
    // what does the ampersand do here ?
    double &operator[](unsigned index);
    cSpectrum &operator=(double value);
    cSpectrum &operator=(const std::vector<double > &newSpec);
    cSpectrum &operator=(const cSpectrum &newSpec);
    cSpectrum &operator+(const cSpectrum &summand);
    // multiplying vectors element wise
    cSpectrum &operator*(const cSpectrum &factor);
    // kind like a scalar for the whole vector ?
    cSpectrum &operator*(double factor);
    cSpectrum &exp();
    double sum();
    void readSpectrum(const std::string &fname, double tubeVoltage,
                      double &minEnergy, std::string &spectrumName);

};

#endif //SPECTRUM_SPECTRUM_H
