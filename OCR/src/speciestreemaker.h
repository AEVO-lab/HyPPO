#ifndef SPECIESTREEMAKER_H
#define SPECIESTREEMAKER_H

#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include "cluster.h"
#include "div/longcounter.h"

using namespace std;

class SpeciesTreeMaker
{
public:
    SpeciesTreeMaker();


    string BuildSpeciesTree(vector<Cluster*> &clusters);


    vector< pair< set<string>, set<string> > > GetAllSplits(set<string> &species);


    //return val = <cost, newick>
    pair<int, string> EvaluateSpeciesSet(set<string> &species, vector<Cluster*> &clusters);

    //return val = <cost, newick>
    pair<int, string> EvaluateSplit(set<string> &s1, set<string> &s2, vector<Cluster*> &clusters, vector<set<string> > &intersections);


    string GetSpeciesSetKey(set<string> species);





    int verbose;
    bool applyCladesFilter;

private:
    unordered_map<string, pair<int, string> > speciesSetCache;


    void BuildPossibleClades(vector<Cluster*> &clusters, set<string> &allSpecies);
    bool IsPossibleClade(set<string> clade);

    set<string> GetComplement(set<string> myset, set<string> universe);

    set<string> GetMostCommonSubclade(vector<Cluster*> &clusters, set<string> &species);

    set<string> GetKeySpecies(string key);

    set< set<string> > possibleClades;

};

#endif // SPECIESTREEMAKER_H
