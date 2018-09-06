#include "longcounter.h"

LongCounter::LongCounter(int nbCells, int base)
{
    this->nbCells = nbCells;
    this->base = base;

    for (int i = 0; i < this->nbCells; i++)
    {
        this->counters.push_back(0);
    }
}


int LongCounter::GetFromRight(int i)
{
    return counters[i];
}

int LongCounter::GetFromLeft(int i)
{
    return counters[ counters.size() - 1 - i ];
}


void LongCounter::Increment()
{
    bool done = false;
    int cindex = 0;

    while (!done)
    {
        if (counters.size() <= cindex)
            done = true;

        counters[cindex] += 1;
        if (counters[cindex] >= base)
        {
            counters[cindex] = 0;
            cindex++;
        }
        else
        {
            done = true;
        }
    }
}


string LongCounter::ToString()
{
    string ret = "";
    for (int i = 0; i < counters.size(); i++)
    {
        if (ret != "")
            ret += " ";
        ret += Util::ToString(this->GetFromLeft(i));
    }
    return ret;
}
