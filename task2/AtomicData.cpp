//
// Created by Hesham Othman on 07.10.19.
//

using namespace std;
#include "AtomicData.h"
#include <iostream>
#include "fstream"
#include <cmath>


// static initializations .. is this right?
vector<tPoint2d> cAtomicData::tcs[100];
float cAtomicData::A[100];
bool cAtomicData::prepared = false;

cAtomicData::cAtomicData() = default;
void cAtomicData::prepare() {
    if (prepared == true) return;

    ifstream inp;
    bool valueFound = false;
    // open file
    inp.open("/Users/Hesham/CLionProjects/xrays-simulations/task2/totalCrossSection.dat", ios::binary);
    if(!inp.is_open()) throw runtime_error(" totalCrossSection.dat couldn't be opened");

    // 1- title
    char title[12];
    inp.read(title, 11);
    title[11] = 0;
    cout << "title: " << title << endl;

    // 2- length of the name Nn
    uint32_t Ns;
    inp.read((char*) &Ns, sizeof(Ns));
    // cout << "Ns: " << Ns << endl;

    // 3- totalCrossSection name
    char *charOfName = new char[Ns+1];
    inp.read(charOfName, Ns);
    charOfName[Ns] = 0;
    // cout << "Total cross section: " << charOfName << endl;

    for (int i = 0 ; i < 100 ; i++) {
        // atomic number Z
        uint32_t Z;
        inp.read((char*) &Z, sizeof(Z));
        // cout << "Z: " << Z << endl;

        // std atomic weight Ar
        float Ar;
        inp.read((char*) &Ar, sizeof(Ar));
        // cout << "Ar: " << Ar << endl;
        A[i] = Ar;

        // number of entries
        uint32_t Nt;
        inp.read((char*) &Nt, sizeof(Nt));
        // cout << "Nt: " << Nt << endl;

        for (int j = 0; j < Nt; j++) {
            // energy
            float en;
            inp.read((char*) &en, sizeof(en));

            // cross section
            float cs;
            inp.read((char*) &cs, sizeof(cs));
            // for debugging reasons
            tPoint2d t2d = {en, cs};

            tcs[i].insert(tcs[i].end(), t2d);
            // cout << j << ": " << en << ", " << cs << endl;
        }
    }
    prepared = true;
}

double cAtomicData::getStdAtomicWeight(unsigned Z) {
    if (Z > 100 || Z < 1 ) throw runtime_error("atomic_number_out_range");
    return A[Z - 1];
}

double cAtomicData::getTotalCrossSection(unsigned Z, double energy) {
    if (Z > 100 || Z < 1 ) throw runtime_error("atomic_number_out_range");
    // removed to save memeory.
    // vector<tPoint2d> values = tcs[Z - 1];
    // energy divided by 1000 as we get it in KeV but vector entries are in MeV
    float sigma_k, sigma_k_plus_1, e_k, e_k_plus_1;
    for( int i = 0; i != tcs[Z - 1].size(); i++) {
        if ((energy/1000) > tcs[Z - 1][i].x) continue;
        e_k = tcs[Z - 1][i-1].x;
        sigma_k = tcs[Z - 1][i-1].y;

        e_k_plus_1 = tcs[Z - 1][i].x;
        sigma_k_plus_1 = tcs[Z - 1][i].y;
        break;
    }
    return sigma_k * ::exp(log(sigma_k_plus_1 / sigma_k) * (log((energy/1000)/e_k)/log(e_k_plus_1/e_k)));
}
