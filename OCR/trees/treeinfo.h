#ifndef TREEINFO_H
#define TREEINFO_H


#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include "node.h"
#include "div/util.h"

#define STATE_SPECIATION 0
#define STATE_DUPLICATION 1
#define STATE_LOSS 2

class Node;

using namespace std;


/**
  A TreeInfo is associated with a tree, and maintains information on the whole tree.
  If a tree was created with maintainTreeInfo = true, any node in the tree can returns the tree info using
  Node::GetTreeInfo().
  The main usage is to find LCA quickly (if the height < 64),
  getting the max degree of the tree, finding a node by label quickly.
  **/
class TreeInfo
{
private:

    unordered_map<string, Node*> labelNodes;
    unordered_map<uint64, Node*> pathBitNodes;
    void GiveNodeSomeLabel(Node* n, int labelCounter);
    Node* root;
    int degree;

public:
    TreeInfo(Node* root);

    /**
      Returns, as expected, the root of the tree
      **/
    Node* GetRoot();

    /**
      Returns the highest number of children a node has in the tree.
      **/
    int GetDegree();

    /**
      Method called automatically by Node when something gets inserted.
      **/
    void OnNodeInserted(Node* n, int pos);


    /**
      Method called automatically by Node when something gets deleted.
      **/
    void OnNodeDeleted();



    /**
      Method called automaically by Node when a label is changed.
      **/
    void OnLabelChanged(Node* n, string oldLabel, string newLabel);

    /**
      Parses a tree to update internal variables of the tree info.
      Useful if tree info and tree got desynchronized.
      **/
    void ParseTree(Node* node = NULL, bool nameEmptyLabels = false, bool computeDepth = false, bool computePathBits = false);

    /**
      Returns the lowest common ancestor of n1, n2 in constant time.
      **/
    Node* GetLCA(Node* n1, Node* n2);

    /**
      Returns the lowest common ancestor of all passed nodes, in O(|nodes|) time.
      **/
    Node* GetLCA(vector<Node*> &nodes);

    /**
      Path bits are used to find LCA in constant time.  Each node has a unique path bit string.
      This is mainly used by GetLCA
      **/
    Node* GetNodeByPathBits(uint64 pathbits);

    /**
      Find a node by its label, in constant time
      WARNING: if multiple nodes share the same label, undefined behavior
      **/
    Node* GetNodeByLabel(string label);

};

#endif // TREEINFO_H
