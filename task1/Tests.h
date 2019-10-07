//
// Created by Hesham Othman on 07.10.19.
//
#include <iostream>
#include "Spectrum.h"
#include <cmath>

using namespace std;

void testSpectrum();
void testMultiplyingWholeSpecs();
void testAddingSpec();
void testMultiplyingSpecWithScalar();
void testExp();
void testSum();

void testTask1() {
        cout << "Test_task1_cSpectrum" << endl;
        testSpectrum();
        testAddingSpec();
        testMultiplyingWholeSpecs();
        testMultiplyingSpecWithScalar();
        testExp();
        testSum();
        cout << "Tests_task1_successful_run" << endl;

        cSpectrum tmp;
        double minEnergy;
        double tubeVoltage = 75;
        string spectrumName;
        tmp.readSpectrum("/Users/Hesham/CLionProjects/Spectrum/task1/SRO33100ROT350.dat", tubeVoltage, minEnergy , spectrumName);
}

void testSpectrum() {
    cSpectrum spec;
    cout << "check_size()_and_resize()" << endl;
    if (spec.size() != 0)
        throw runtime_error("Unexpected_size_in_spectrum");
    spec.resize(100);
    if (spec.size() != 100)
        throw runtime_error("Unexpected_size_in_spectrum_after_setting_the_value");
    spec[40] = 90;
    if(spec[40] == 90)
        cout << "operate [] works" << endl;
    cout << "__finished" << endl;
}

void testAddingSpec() {
    cSpectrum spec1;
    cSpectrum spec2;
    cSpectrum result;
    spec1 = {1, 3, 5.3, 23.2, 5};
    spec2 = {4, 6, 1.2, 41.3, 4};
    result = {5, 9, 6.5, 64.5, 9};
    cSpectrum resultAddition =  spec1 + spec2;
    if ( resultAddition.size() != result.size() )
        throw runtime_error(" adding does not work properly, size problem");
    for (unsigned i = 0; i < result.size() ; i++) {
//        printf("%.20lf %.20lf\n", resultAddition[i], result[i]);
        if (fabs(resultAddition[i]/result[i] - 1) > 1e-12)
            throw runtime_error("precision error");
    }
    cout << "+ works" << endl;
}

void testMultiplyingWholeSpecs() {
    cSpectrum spec1;
    cSpectrum spec2;
    cSpectrum resultExpected;
    spec1 = {1, 3, 5.3, 23.2, 5};
    spec2 = {4, 6, 1.2, 41.3, 4.3};

    resultExpected = {4, 18, 5.3*1.2, 23.2*41.3, 5*4.3};
    cSpectrum resultMul =  spec1 * spec2;
    for (unsigned i = 0; i < resultMul.size() ; i++) {
        if (fabs(resultMul[i]) != fabs(resultExpected[i]))
            throw runtime_error(" multiplying does not work properly");
    }
    cout << "* works" << endl;
}

void testMultiplyingSpecWithScalar() {
    cSpectrum spec1;
    cSpectrum resultExpected;
    spec1 = {1, 3, 5.3, 23.2, 5};

    resultExpected = {3, 9, 3*5.3, 3*23.2, 15};
    cSpectrum resultMulScalar =  spec1 * 3;

    for (unsigned i = 0; i < resultMulScalar.size() ; i++) {
        if (fabs(resultMulScalar[i]) != fabs(resultExpected[i]))
            throw runtime_error(" scalar multiplying does not work properly");
    }
    cout << "* scalar works" << endl;
}

void testExp() {
    cSpectrum spec1;
    spec1 = {1, 3, 5.3, 3.2, 5};

    cSpectrum resultExp = spec1.exp();
    if (resultExp[1] != ::exp(3))
        throw runtime_error(" exp does not work properly");
    cout << "e works" << endl;
}

void testSum() {
    cSpectrum spec1;
    spec1 = {1, 3, 5, 3, 5};
    if (spec1.sum() != 17 )
        throw runtime_error(" sum does not work properly");
    cout << "sum works" << endl;
}