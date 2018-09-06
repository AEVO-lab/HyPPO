#ifndef SUPERGENETREEMAKER_H
#define SUPERGENETREEMAKER_H

#include <string>
#include <vector>
#include <unordered_map>

#include "trees/node.h"
#include "div/util.h"
#include "trees/genespeciestreeutil.h"


using namespace std;

//NOTE: trees must have internal nodes labeled uniquely to use this class
class TreeLabelIntersectionInfo
{
public:
    TreeLabelIntersectionInfo();

    void ComputeAllIntersections(vector<Node*> trees);

    void ComputeIntersections(Node *tree1, Node *tree2);

    void AddIntersection(string lbl1, string lbl2);

    bool Intersect(string lbl1, string lbl2);

    bool Intersect(vector<Node*> trees);

    bool IsPartitionIntersecting(vector<Node*> treesLeft, vector<Node*> treesRight);

private:
    unordered_set<string> intersecting_nodes;   //key is inter1-inter2 (separated by this)
    string key_separator;



};

class SuperGeneTreeMaker
{
public:
    SuperGeneTreeMaker();

    //returns a supertree + DL cost
    pair<Node*, int> GetSuperGeneTreeMinDL(vector<Node *> &trees, vector<unordered_map<Node *, Node *> > &lca_mappings,
                                       Node* speciesTree, bool mustPreserveDupSpec, bool isFirstCall = true);

private:


    void ApplyNextConfig(vector<int> &counters);

    TreeLabelIntersectionInfo intersectionInfo;

    unordered_map<string, pair<Node*, int> > recursionCache;


};

#endif // SUPERGENETREEMAKER_H
