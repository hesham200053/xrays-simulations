//
// Created by Hesham Othman on 18.11.19.
//



#include "random.h"

void testRandom();
void testRandom() {
    cRandKiss kiss((uint32_t)time(NULL));
    cRandNormal normal(kiss);
    normal.prepare();
    cout << normal.rand() << endl;
    unsigned data[100] = {0};
    double randoms[1000000] = {0};
    for (int i = 0; i < 1000000 ; ++i) {
        double r = normal.rand();
        randoms[i] = r;
    }

}
