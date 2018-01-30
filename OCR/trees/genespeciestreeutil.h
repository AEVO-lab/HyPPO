#ifndef GENESPECIESTREEUTIL_H
#define GENESPECIESTREEUTIL_H

#include "trees/node.h"
#include "trees/newicklex.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <iostream>

using namespace std;

class GeneSpeciesTreeUtil
{
private:
    GeneSpeciesTreeUtil();

public:
    static GeneSpeciesTreeUtil* Instance()
    {
        static GeneSpeciesTreeUtil instance;

        return &instance;
    }

    int LASTNBDUPS;
    int LASTNBLOSSES;

    unordered_map<Node*, Node*> GetLCAMapping(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> &geneLeavesSpeciesMapping);

    unordered_map<Node*, Node*> GetLCAMapping(Node *geneTree, Node *speciesTree, string geneLabelSeparator, int speciesIndex);

    unordered_set<Node*> GetGeneTreeSpecies(Node *geneTree, unordered_map<Node*, Node*> &lcaMapping);

    vector<Node*> GetGenesSpecies(vector<Node*> genes, unordered_map<Node*, Node*> &lcaMapping);

    unordered_map<Node*, Node*> GetGeneSpeciesMappingByLabel(Node* geneTree, Node* speciesTree, string separator = "_", int speciesIndex = 0);

    Node* CopyTreeWithNodeMapping(Node* tree, unordered_map<Node*, Node*> &yourMapping, unordered_map<Node*,Node*> &mappingToFill);

    bool HaveCommonSpecies(Node* tree1, Node* tree2, unordered_map<Node*, Node*> &mapping);

    Node* GetSingleNodeLCAMapping(Node* n, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);

    void PrintMapping(Node* g, unordered_map<Node*, Node*> &mapping);

    vector<Node*> GetNADNodes(Node* g, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);


    void RelabelGenes(Node* geneTree, string search = ";;", string replace = "__");

    void RelabelGenesByIndex(Node* geneTree, string separator, int indexToKeep);

    int GetDLScore(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> lcaMapping);


    int GetNbLossesOnBranch(Node* speciesDown, Node* speciesUp, bool isDupTop);


    void PruneSpeciesTreeFromLCAMapping(Node* speciesTree, Node* geneTree, unordered_map<Node *, Node *> lca_mapping);
    void PruneSpeciesTree(Node* speciesTree, set<Node*> speciesToKeep);

    vector<Node*> GetGeneTreeHighestSpeciations(Node* geneTree, Node* speciesTree,
                                                 unordered_map<Node*, Node*> lca_mapping);

    /**
     * @brief IsNodeDup Checks if parsimony would infer a duplication at geneTreeNode.  Works in the non-binary case.
     * @param geneTreeNode
     * @param lca_mapping
     * @return
     */
    bool IsNodeDup(Node* geneTreeNode, unordered_map<Node*, Node*> &lca_mapping);

    void LabelInternalNodesWithLCAMapping(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> lca_mapping);
    void LabelInternalNodesUniquely(Node* tree);
    void LabelInternalNodesUniquely(vector<Node*> trees);

    string GetPrunedSpeciesTreeNewick(string gcontent, string scontent);
};

#endif // GENESPECIESTREEUTIL_H
