#ifndef TREEITERATOR_H
#define TREEITERATOR_H

#include "node.h"
#include <unordered_set>


using namespace std;

class Node;


/**
  Base class for a tree iterator.
  **/
class TreeIterator
{
public:
    TreeIterator(Node* root);

    /**
      Set leavesOnly to true if you want the iterator to only go through leaves
      **/
    void SetLeavesOnly(bool leavesOnly);

    /**
      Set the leaves of interest - tree will only iterate through leaves in the set or their ancestors
      **/
    void SetLeaves(unordered_set<Node*>* leaves);

    /**
      Returns the next node in the iteration list
      **/
    virtual Node* next();

    /**
      Removes the current node from the tree, and returns the next one in the iteration.
      By deleting a node n, the parent of n becomes the parent of n's children.
      **/
    virtual Node* DeleteCurrent();
protected:
    Node* curNode;
    Node* root;
    bool leavesOnly;
    unordered_set<Node*>* leaves;

    bool IsLeaf(Node* n);

};


/**
  As implied by the class name, iterates over the nodes in a post order traversal order.
  Set Node::GetPostOrderIterator() for usage example
  **/
class PostOrderTreeIterator : public TreeIterator
{
public:
    PostOrderTreeIterator(Node* root);

    Node* next();
    Node* DeleteCurrent();

};


/**
  As implied by the class name, iterates over the nodes in a preorder traversal order.
  DeleteCurrent is not implemented.
  **/
class PreOrderTreeIterator : public TreeIterator
{
public:
    PreOrderTreeIterator(Node* root);

    Node* next();


};






#endif // TREEITERATOR_H
