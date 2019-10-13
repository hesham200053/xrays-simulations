#include <iostream>
#include "task1/Tests.h"
#include "task2/Tests.h"

using namespace std;

int main() {
    try {
        // task 1
        // testTask1();
        // task 2
        cout << "Test_task2_cMaterial" << endl;
        testMaterial();
        cAtomicData cAtomicData;
        cAtomicData.prepare();
        double sigma_Aluminum  = cAtomicData.getTotalCrossSection(13, 0.0279);
        cout << sigma_Aluminum << endl;

        cMaterial water;
        water.setName("water particle");
        water.setDensity(1);
        water.addElement(1, 2);
        water.addElement(8, 1);

        // for energy level of 50KeV or 0.05MeV
        double totalAttCoeff = water.getAttCoeff(0.05);
        cout << "Total attenuation coefficient: " << totalAttCoeff << " cm e-3 .b"<< endl;
        cout << "Tests_task2_successful_run" << endl;


    } catch (exception &exc) {
        cout << "ERROR:_" << exc.what() << endl;
    } catch (...) {
        cout << "ERROR:_Unexpected_error" << endl;
    }
    return 0;
}