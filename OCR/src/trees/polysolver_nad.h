#ifndef NADCORRECTOR_H
#define NADCORRECTOR_H

#include "trees/node.h"
#include "trees/genespeciestreeutil.h"

#include "trees/newicklex.h"

#include <map>
#include <unordered_map>

#include <iostream>

#include <utility>

using namespace std;

class SADNADNode;
/**
  A SADNADGraph stands for [S][AD][NAD]Node.
  It turns each subtree of a polytomy into a node, connecting them by the S/AD/NAD
  relationship between them.

  Note : this is a naive, slow implementation.
  **/
class SADNADGraph
{
public:
    SADNADGraph();
    ~SADNADGraph();

    void AddNode(Node* n);


    //set<Node*> GetADComponentOf(Node* n);

    /**
      Each node of the graph is part of an AD connected component.
      We store the components by assigning them a canonical node
      (a node that always represents the AD components), which gets
      updated as the graph evolves.
      This method returns the representant of the AD-component containing n
      **/
    Node* GetADComponentRepresentantOf(Node* n);

    /**
      Tells if there is an AD path between n1 and n2
      **/
    bool HaveSameADComponent(Node* n1, Node* n2);


    /**
      Create the SADNADNode objects based on the nodes in the vector (roots of disjoint subtrees) and link them by the appropriate edge type.
      This must be called before the graph can be used.
      **/
    void BuildGraph(vector<Node*> nodes, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);

    /**
      Applies a join on n1, n2 and updates the graph.  This removes n1, n2 from the graph and creates a new node p
      that is the parent of n1 and n2.
      **/
    void MergeNodes(Node* n1, Node* n2, Node* newNode, unordered_map<Node*, Node*> &lcaMapping);


    /**
      Gene tree subtrees are mapped to a SADNADNode object.  This returns the map.
      **/
    const map<Node*, SADNADNode*> GetNodes(){return nodes;}

    int GetNbADComponents();

    /**
      Get the number of AD-components after we add the specified additionalSEdges
      **/
    int GetNbADComponents(set<pair<Node*, Node*> > &additionalSEdges);

    vector<pair<Node *, Node *> > GetUsefulSpeciationEdges();

    void PrintGraph();

    /**
      true iff n1 and n2 are related by speciation
      **/
    bool HasSEdge(Node* n1, Node* n2);

private:
    map<Node*, SADNADNode*> nodes;
    void AddEdge(Node* n1, Node* n2, string edgeType);
    string GetEdgeType(Node* n1, Node* n2, unordered_map<Node*, Node*> &lcaMapping);

    //map<Node*, set<Node*> > ADComponents;   //the key is a representant of the component

};



/**
  SADNADNodes are used internally by SADNADGraph mainly to hold the relationship between nodes.
  **/
class SADNADNode
{
public:
    set<Node*> S_Neighbors;
    set<Node*> AD_Neighbors;
    set<Node*> NAD_Neighbors;

    Node* ADComponent;

    string GetEdgeTypeWith(Node* n)
    {
        if (S_Neighbors.find(n) != S_Neighbors.end())
            return "S";
        if (AD_Neighbors.find(n) != AD_Neighbors.end())
            return "AD";
        if (NAD_Neighbors.find(n) != NAD_Neighbors.end())
            return "NAD";
        throw "No edge defined !";
    }

    void RemoveNeighbor(Node* n)
    {
        S_Neighbors.erase(n);
        AD_Neighbors.erase(n);
        NAD_Neighbors.erase(n);
    }
};


/**
  Object returned by a correction.  Contains
  correction : the actual corrected tree
  firstPolySize and secondPolySize : the number of children of the 2 polytomies created
  nadCladeGenes : the original leaves under the node we corrected
  **/
class PolySolverCorrectionInfo
{
public:
    Node* correction;
    int firstPolySize;
    int secondPolySize;
    vector<string> nadCladeGenes;
};


class PolySolverNAD
{
public:
    PolySolverNAD();





    /**
      Returns a copy of geneTree in which the highest NAD is corrected, or NULL if no correction is applied.
      The returned value must be deleted by caller, unless NULL.
      **/
    PolySolverCorrectionInfo CorrectHighestNAD(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping);


    /**
      Makes a copy of geneTree, then applies correction to node (assuming it is a NAD - otherwise problems will arise).
      **/
    PolySolverCorrectionInfo CorrectNodeByMultifurcation(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping, Node* node);


    /**
      Generates a random polytomy with k children, and a random species tree and compares the heuristic with the optimal solution.
      **/
    void PerformRandomTest(int k, int verbose = 0);


    /**
      Tries all possibilities to find the best clique set that connects all AD components.  Returns the minimum # of AD components we can make.
      Takes time.
      **/
    int FindBestSolution(SADNADGraph &graph, vector<pair<Node*, Node*> > availableSpecs, set<pair<Node*, Node*> > chosenSpecs, int verbose = 0);

    /**
      For debugging purposes
      **/
    int bestSolSoFar;
    set<pair<Node*, Node*> > bestChoiceSoFar;
     int last_nb_ad_components;  //use with care after calling SolvePolytomy
     pair<Node*, Node*> GetRandomPolytomy(int k, int verbose);
     Node* SolvePolytomy(vector<Node*> leaves, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);
private:

     Node* SolvePolytomies(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping);


     Node* SolvePolytomy(Node* polytomyRoot, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);

    pair<Node*, unordered_map<Node*, Node*> > PolytomizeNAD(Node* nadNode, Node* speciesTree, unordered_map<Node*, Node*> lcaMapping);

    void SpeciateCleverly(SADNADGraph &graph, set<Node*> &gLeftGuys, set<Node*> &gRightGuys, Node* curSpecies, unordered_map<Node*, Node*> &lcaMapping);

    vector<vector<Node*> > GetSortedLocalADComponents(SADNADGraph &graph, set<Node*> guys, set<Node*> otherGuys);





};

#endif // NADCORRECTOR_H
