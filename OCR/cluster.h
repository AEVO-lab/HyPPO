#ifndef CLUSTER_H
#define CLUSTER_H


#include <unordered_map>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>

#include "div/util.h"
#include "trees/node.h"
#include "trees/genespeciestreeutil.h"

using namespace std;

class GeneTree
{

public:

    GeneTree();

    GeneTree(Node* geneTreeRoot, Node* speciesTree, unordered_map<Node*, Node*> &lca_mapping);

    bool IsDupUnderLCAMapping(Node* node);

    Node* GetLCAMapping(Node* geneTreeNode);


    Node* speciesTree;
    Node* geneTreeRoot;
    unordered_map<Node*, Node*> lca_mapping;
};


class Cluster
{

public:

	static string speciesSeparator;
	static int speciesIndex;

    Cluster(vector<string> &geneNames, Node* speciesTree);

    ~Cluster();

    string GetGeneSpecies(string g);

    set<string> GetSpeciesSet();

    int GetNbGenes();

    string GetGeneName(int i);

    string GetSpeciesName(string geneLabel);

    Node* GetLCASpecies();

    GeneTree GetGeneTree();

    void BuildGeneTree();

    void SetGeneTree(Node* gtree, Node* speciesTree);


    void UpdateLCAMapping();


    string ToString();

    bool HasSpecies(string sp);

    set<string> GetSpeciesIntersection(set<string> species);

    Cluster* GetCompressedCluster(map<string, set<string> > &blocks, string new_gene_name);


    static set<string> GetAllSpecies(vector<Cluster*> clusters);
    static map<string, set<string> > GetBlocks(vector<Cluster *> &clusters);

		
private:

    void Init(vector<string> &geneNames, Node* speciesTree);

    vector<string> geneNames;
    map<string, string> geneToSpecies;
    set<string> speciesSet;
    Node* speciesTree;
    Node* lcaSpecies;
    GeneTree geneTree;

    map< string, set<string> > blocks;


};




#endif // CLUSTER_H
