#include <iostream>
#include "task1/Tests.h"
#include "task2/Tests.h"
#include "task3/Tests.h"
#include "task4/Tests.h"
#include "task5/Tests.h"

//void testMedImage();
void testRandom();
using namespace std;

int main() {
    try {
//        testTask1();
//        testTask2();
//        testParticle();
//        testMinDestToOrigin();
//        testH2O();
//        testTask2();
//        testParticle();
//        testMedImage();
        testSimulation();
//        testRandom();

    } catch (exception &exc) {
        cout << "ERROR:_" << exc.what() << endl;
    } catch (...) {
        cout << "ERROR:_Unexpected_error" << endl;
    }
    return 0;
}

