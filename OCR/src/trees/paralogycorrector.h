#ifndef PARALOGYCORRECTOR_H
#define PARALOGYCORRECTOR_H

#include "node.h"




#include <map>
#include <string>
#include <vector>
#include <set>


#include "newicklex.h"


using namespace std;


class ParalogyCorrector
{

private:



    /**
      This is used to correct subtrees recursively - does the same as correct gene tree, except that the lcaMapping, forbiddenClades and
      map of genes by name are given so we don't recompute them each time.
      **/
    Node* CorrectGeneTree(Node* geneTree, Node* speciesTree, map<Node*, Node*> &geneSpeciesMapping, vector<pair<string, string> > &orthologs, map<Node*, Node*> &lcaMapping,
                          set<Node*> &forbiddenClades, map<string, Node*> &genesByName);


    /**
      Solves the subproblem of CorrectGeneTree : given a set of subtrees, reconnect them in order to maximize orthology relationships.
    **/
    Node* GetMaxOrthologyTree(set<Node*> &subtrees, map<Node*, Node*> &lcaMapping, Node* speciesTree);

    /**
      Makes orthlogs symmetric by adding everything that's missing
    **/
    void EnsureOrthologsSymnmetry(vector<pair<string, string> > &orthologs);

    /**
        Does a post order traversal from current and fills highest with the highest preservable nodes (the first encoutered that are not forbidden)
    **/
    void FindHighestPreservableNodes(Node* current, set<Node*> &highest, set<Node*> &forbiddenNodes);


    //Used by an old version of the algorithm - kept for historical purposes
    //void DoDFS(set<Node*> &unvisited, set<Node*> &curCC, map<Node*, set<Node*> > edges, Node* lastVisited);

public:
    ParalogyCorrector();

    /**

      Frees geneTree from bad paralogs, given orthologs, the list of gene names that are required to be orthologs.
      The returned Node* is the root of the corrected tree, and needs to be deleted by caller.
      This method takes care of reconciling the geneTree and identifying bad duplications.

      Node* geneTree : the root of the gene tree to correct
      Node* speciesTree : the root of the species tree
      map<Node*, Node*> &geneSpeciesMapping : a map where the key is a gene leaf, the value is a species tree leaf.  It's the mapping between genes and their species
      vector<pair<string, string> > &orthologs : the list of required orthologs, given by gene name pairs which are the labels of gene tree leaves.  Doesn't need to be symmetric.

      **/
    Node* CorrectGeneTree(Node* geneTree, Node* speciesTree, map<Node*, Node*> &geneSpeciesMapping, vector<pair<string, string> > &orthologs);


};

#endif // PARALOGYCORRECTOR_H
