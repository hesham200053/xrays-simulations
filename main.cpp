#include <iostream>
#include "task1/Tests.h"
#include "task2/Tests.h"
#include "task3/Tests.h"

using namespace std;

int main() {
    try {
        // task 1
        // testTask1();
//        testTask2();
        testParticle();
//        testMinDestToOrigin();
//        testH2O();
        // task 2
        // testTask2();
        // testParticle();

    } catch (exception &exc) {
        cout << "ERROR:_" << exc.what() << endl;
    } catch (...) {
        cout << "ERROR:_Unexpected_error" << endl;
    }
    return 0;
}