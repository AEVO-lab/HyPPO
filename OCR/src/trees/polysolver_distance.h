#ifndef POLYSOLVER_DISTANCE_H
#define POLYSOLVER_DISTANCE_H

#include "trees/node.h"

#include <map>
#include <vector>
#include <unordered_map>

#include <string>

#include "genespeciestreeutil.h"
#include "newicklex.h"

using namespace std;

class PSDCachedResolution;


/**
  PSD stands for PolySolverDistance.

  **/
class PSDCacheManager
{
public:

    unordered_map<string, PSDCachedResolution*> cache;


    static PSDCacheManager* Instance()
    {
        static PSDCacheManager instance;

        return &instance;
    }

    void CacheResolution(Node* orig, Node* resolved, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);

    string GetKey(Node* orig);

    PSDCachedResolution* GetCachedResolution(Node* orig, Node* speciesTree);
    PSDCachedResolution* GetCachedResolution(string key);

    void CleanUp();

};


class PSDCachedResolution
{
public:
    Node* speciesTree;
    unordered_map<Node*,Node*> lcaMapping;
    Node* resolution;

    ~PSDCachedResolution()
    {
        delete resolution;
    }
};

/**
  Information specific to each node in the temp trees we're correcting
  **/
class PSDNodeInfo : public NodeInfo
{
public:
    virtual void CopyFrom(NodeInfo* n);
    virtual NodeInfo* GetClone();

    string psdName;
    double totalBranchLength;
    bool isLoss;
    int dlScore;
};

/**
  Represents a cell of the dynamic programming table we're building.
  Each cell holds its trees, which are then given to subsequent cells that will need to join these trees.
  Either the trees were formed by k speciations, by taking k + 1 trees and joining two under a dup, or
  taking k - 1 trees adding a loss.  In the first case, we store the trees in spec_subtrees,
  in the two latter cases (which are mutually exclusive), we store them in duploss_subtrees.
  For now, one of them is filled.  In the future, both might be filled if the two choices lead
  to equivalent solutions.
  **/
class PSDCell
{
public:
    vector<Node*> duploss_subtrees;
    vector<Node*> spec_subtrees;
    int CVal;
    int MVal;

    ~PSDCell()
    {
        for (int i = 0; i < duploss_subtrees.size(); i++)
        {
            delete duploss_subtrees[i];
        }
        duploss_subtrees.clear();

        for (int i = 0; i < spec_subtrees.size(); i++)
        {
            delete spec_subtrees[i];
        }
        spec_subtrees.clear();

    }

    vector<Node*> GetDefaultSubtrees()
    {
        if (spec_subtrees.size() > 0)
            return spec_subtrees;
        return duploss_subtrees;
    }
};


/**
  The class that solves a polytomy minimizing dup/loss and joining trees in an NJ fashion.
  **/
class PolySolverDistance
{

private:
    map<Node*, map<Node*, double> > GetSubtreesDistances(vector<Node*> subtrees, map<Node*, map<Node*, double> > &distances);

    vector<Node*> MakeSpecSubtrees(Node* s, int k, set<Node*> speciesLeaves, map<Node*, vector<Node*> > &speciesGenes, map<Node*, map<int, PSDCell*> > &cells, unordered_map<Node*, Node*> &lcaMapping,  map<Node*, map<Node*, double> > &distances);
    vector<Node*> MakeOneDuplication(Node* s, int k, set<Node*> speciesLeaves, map<Node*, vector<Node*> > &speciesGenes, map<Node*, map<int, PSDCell*> > &cells, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances);
    vector<Node*> MakeOneLoss(Node* s, int k, map<Node*, vector<Node*> > &speciesGenes, map<Node*, map<int, PSDCell*> > &cells, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances);


    void PrintScoresTable(Node* speciesTree, int maxK, map<Node*, map<int, PSDCell*> > &cells);

    //for each polytomy, we need to work with the 'linked species tree', which has no useless leaves.  To see this, try the algorithm on G = ((a,b), c,d), S = ((a,b), (c,d))
    set<Node*> ComputeLinkSpeciesLeaves(vector<Node*> polytomyChildren, unordered_map<Node*, Node*> &lcaMapping, Node* speciesTree);
public:

    //Node* SolvePolytomies(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping, map<Node*, map<Node*, double> > &distances);

    PolySolverDistance()
    {
        verbose = 0;
        useCache = true;
    }

    Node* SolvePolytomies(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances);

    /**
      Makes a binary tree out of the given leaves that minimizes the DL score.  The trees are joined using the distances.
      If leaves.size() = n, distances should contain n maps of n - 1 elements.
      Distances between n1, n2 are accessed through distances[n1][n2] OR distances[n2][n1], so both sides need to be filled.
      **/
    Node* SolvePolytomy(vector<Node*> leaves, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances);

    int verbose;
    bool useCache;


};

#endif // POLYSOLVER_DISTANCE_H
