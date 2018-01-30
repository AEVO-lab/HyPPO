#ifndef NEWICKLEX_H
#define NEWICKLEX_H

#include <string>
#include "node.h"
#include "div/util.h"

#include <iostream>
#include <set>

using namespace std;

class Node;

class NewickLex
{
public:

    /**
      Takes a Newick string and returns the root of a new tree.\n
      Does not try to validate anything, and assumes format correctness.\n
      The accepted Newick format is basic : it doesn't consider special parameters
      such as branch length or anything of the sort - the whole string corresponding to a node
      in the Newick string
      becomes the label of a node in the tree.
      User has to delete returned value. \n
      Set maintainTreeInfo = true if you want each node of the tree to hold a treeInfo,
      which mainly serves to accelerate LCA finding and searching a node by label. \n
      See Node constructor and TreeInfo class for more information.
    **/
    static Node* ParseNewickString(string& str, bool maintainTreeInfo = false);

    /**
      Converts a tree to a Newick string, naming the nodes using Node::GetLabel().
      If addBranchLengthToLabel is true, the branch length will be appended to the outputted label
      for each node (with a "-" between the label and the branch length)
    **/
    static string ToNewickString(Node* root, bool addBranchLengthToLabel = false, bool addInternalNodesLabel = true);



    static string GetCaterpillarNewick(set<string> labels);

private:
    static int ReadNodeChildren(string &str, int revstartpos, Node* curNode);

    static void WriteNodeChildren(string &str, Node* curNode, bool addBranchLengthToLabel, bool addInternalNodesLabel);

    static void ParseLabel(Node* node, string label);
};



#endif // NEWICKLEX_H
