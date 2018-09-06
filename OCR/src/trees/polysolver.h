#ifndef RECONCILERV2_H
#define RECONCILERV2_H


#include <unordered_map>
#include <string>
#include <sstream>
#include <unordered_set>
#include <map>
#include "node.h"
#include "treeiterator.h"

#include <queue>

using namespace std;

class PolySolver;


/**
  This class essentially represents a row in the dynamic programming table,
  defined by nb(s), break1 and break2 (see paper for details).
  It will also hold the nb of dups/losses required for the species it corresponds to.
  **/
class PolyInfo
{
public:
    int nb;
    int break1;
    int break2;
    int minval;
    int dups;
    int losses;

    PolyInfo()
    {
        nb = 0; break1 = 0; break2 = 0; minval = 0;
    }

    int GetMinCost(int k);
    string ToString();
};

/**
  For a given set of genes, this class is a wrapper around a species tree,
  and contains the functionality necessary to traverse the nodes of the species tree
  that are relevant to the set of genes.
  **/
class LinkedSpeciesTreeUtil
{
private:
    Node* polytomyRoot;
    Node* speciesTreeRoot;
    PolySolver* rec;
    unordered_set<Node*> speciesLeaves;
    unordered_set<Node*> traversedNodes;


    void FindLeaves();
    void AddLeaf(Node* n);
    bool AddInternalNode(Node* n);

public:
    LinkedSpeciesTreeUtil(Node* polytomyRoot, PolySolver* rec);
    //unordered_set<Node*>* GetSpeciesLeaves();
    PostOrderTreeIterator* GetPostOrderIterator();
    PreOrderTreeIterator* GetPreOrderIterator();
    void CloseIterator(TreeIterator* it);
    bool IsLeaf(Node* n);
    Node* GetSpeciesTreeRoot();
};


/**
  This class is a singleton, and has only one main method : SolvePolytomies(),
  which finds a resolution for a non-binary gene tree.
  **/
class PolySolver
{
private:
    unordered_map<Node*, Node*> genesMapping;
    unordered_map<Node*, Node*> genesResolutions;


    /**
      Solves a single polytomy and returns the result.  The returned value needs deletion.
      **/
    Node* SolvePolytomy(Node* polyroot);
    /**
      Build a resolution according to the values computed in polyinfos, for each node of the linked species tree
      **/
    Node* BuildResolution(Node* polyroot, LinkedSpeciesTreeUtil &sp, unordered_map<Node*, PolyInfo*> &polyinfos);

    /**
      Computes the number of dups and losses we'll need to infer
      for each node of the species tree.
      It basically corresponds to backtracking the costs table in order to learn these.

      s is the current species (row)
      k is the number of subtrees rooted at s we want after duplicating/adding losses with s
      This recurses into children of s...the initial call should be done with the root of S with k=1
      (we want one binary subtree rooted at the root of s in the end)
      **/
    void CalcDupLoss(Node* s, int k, LinkedSpeciesTreeUtil &sp, unordered_map<Node*, PolyInfo*> &polyinfos);

    PolySolver(){}
public:



    static PolySolver* Instance()
    {
        static PolySolver instance;

        return &instance;
    }

    ~PolySolver();

    /**
      Finds a binary resolution of geneTree in which the dup + loss cost is minimized.
      geneLeavesSpeciesMapping is a map from a gene tree leaf to a species tree node
      Caller has to delete returned pointer.
      WATCH OUT : the species tree must have been created with maintainTreeInfo = true
      **/
    Node* SolvePolytomies(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping);

    /**
      Returns the node of the species tree that g is mapped to.
      Should not be used by user - serves to LinkedSpeciesTreeUtil
      **/
    Node* GetGeneMapping(Node* g);
    //int GetNbLosses();
    //int GetNbDuplications();



    /**
      Tries to determine the gene/species mapping automatically, if your labels have a standard format.
      Say that speciesTree has a leaf labeled "SPECIESX", and geneTree has a leaf labeled "SPECIESX_MYGENE_NAME".
      Then this function splits SPECIESX_MYGENE_NAME by the given separator (default "_"),
      and checks if the first split item matches a species tree label.
      In this case, SPECIESX_MYGENE_NAME maps to SPECIESX because the first substring before the first "_" is SPECIESX.
      **/
    unordered_map<Node*, Node*> GetGeneSpeciesMappingByPrefix(Node* geneTree, Node* speciesTree, string separator = "_");


    /**
      If a node n in geneTree has a branch value < minThreshold, this function contracts n with its parent.
      **/
    void RestrictGeneTreeByBranchSupport(Node* geneTree, double minThreshold);
};


/**
  Function used by RestrictGeneTreeByBranchSupport.  Do not call directly !
  **/
bool PolySolver__RestrictFunction(Node* n, void* args);


#endif // RECONCILERV2_H
