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
        cout << "Tests_task2_successful_run" << endl;


    } catch (exception &exc) {
        cout << "ERROR:_" << exc.what() << endl;
    } catch (...) {
        cout << "ERROR:_Unexpected_error" << endl;
    }
    return 0;
}