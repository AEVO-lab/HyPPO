#ifndef LONGCOUNTER_H
#define LONGCOUNTER_H

#include <vector>
#include <string>
#include "div/util.h"

using namespace std;

class LongCounter
{
public:
    LongCounter(int nbCells, int base);

    void Increment();

    int GetFromRight(int i);
    int GetFromLeft(int i);

    string ToString();

private:
    int nbCells;
    vector<int> counters;
    int base;
};

#endif // LONGCOUNTER_H
