//
// Created by Hesham Othman on 07.10.19.
//

#ifndef SPECTRUM_TESTS_H
#define SPECTRUM_TESTS_H

#include <fstream>
#include "Material.h"
#include "AtomicData.h"

void testMaterial();
void testH2O();

void testTask2() {

    cout << "Test_task2_cMaterial" << endl;
    testMaterial();
    cAtomicData cAtomicData;
    cAtomicData.prepare();
    // values should be in KeV
    double sigma_Aluminum  = cAtomicData.getTotalCrossSection(13, 27.9);
    cout << sigma_Aluminum << endl;

    cMaterial water;
    cSpectrum tmp;
    double minEnergy;
    double tubeVoltage = 75;
    string spectrumName;
    tmp.readSpectrum("/Users/Hesham/CLionProjects/xrays-simulations/task1/SRO33100ROT350.dat", tubeVoltage, minEnergy , spectrumName);


    water.addElement(1, 2);
    water.addElement(8, 1);
    water.setDensity(1.0);
    water.setName("water particle");

    // values should be in KeV
    double totalAttCoeff = water.getAttCoeff(50);
    cout << "Total attenuation coefficient: " << totalAttCoeff << " cm e-3 .b"<< endl;

    cSpectrum attSpecClass;
    vector<double> attSpec;
    water.getAttSpec(attSpec, 10, 75, tmp.size());
    attSpecClass = attSpec;

    cSpectrum I = attSpecClass;
    I *= -0.01;
    I.exp();
    I *= tmp;
    cout << "sum behind water: " << I.sum() << endl;
    cout << "attSpec size: " <<  I.size() << endl;

//    for (auto &v: attSpec)
//        cout << v << endl;
    cout << "Tests_task2_successful_run" << endl;
}
void testMaterial() {
    cMaterial material;
    cAtomicData atomicData;
}

//=============================================================================
void testH2O() {
    cSpectrum spec;
    double voltage = 75;
    string specName;
    double minEnergy;
    cMaterial H2O;

    cout << "+--------------------------------------+" << endl;
    cout << "| test class cMaterial with pure water |" << endl;
    cout << "+--------------------------------------+" << endl;

    // read spectrum
    spec.readSpectrum("/Users/Hesham/CLionProjects/xrays-simulations/task1/SRO33100ROT350.dat", voltage, minEnergy, specName);
//    cout << "spectrum name: " << specName << endl;
//    cout << "minimum energy: " << minEnergy << endl;

    // prepare material
    H2O.addElement(1, 2);
    H2O.addElement(8, 1);
    H2O.setDensity(1.0);
    H2O.setName("pure water");

    // get attenuation spectrum
    cSpectrum attSpec;
    vector<double> tmp;
    H2O.getAttSpec(tmp, minEnergy, voltage, spec.size());
    attSpec = tmp;
    // right answer = 1.26E+06
    // x-ray spectrum behind 1 cm water
    cSpectrum I = attSpec;
    I *= -0.01;
    I.exp();
    I *= spec;
    cout << "sum behind water: " << I.sum() << endl;
    cout << "attSpec size: " <<  attSpec.size() << endl;

    // write spectrum
    ofstream out("/Users/Hesham/CLionProjects/xrays-simulations/task1/specBehindWater.csv");
    if (!out.is_open()) throw runtime_error("Could not open file 'specBehindWater.csv' for writing.");
    out << "energy,intensity" << endl;
    out << "(keV),(#/keV mAs µsr)" << endl;
    for (unsigned k = 0; k < I.size(); k++) {
        out << minEnergy + k * (voltage - minEnergy) / I.size() << ",";
        out << I[k] << endl;
    }
    out << voltage << ",0" << endl;
    out.close();
    cout << "spectrum written in file 'specBehindWater.csv'." << endl;
}


#endif //SPECTRUM_TESTS_H
