#ifndef NODE_H
#define NODE_H

#include "div/define.h"

#include <vector>
#include <string>
#include "treeiterator.h"
#include "treeinfo.h"
#include <set>
#include <unordered_map>

class TreeIterator;
class TreeInfo;

using namespace std;




/**
  Inheritable class used for Node::NodeInfo
  **/
class NodeInfo
{
public:
    virtual void CopyFrom(NodeInfo* n) = 0;
    virtual NodeInfo* GetClone() = 0;
};


/**
  A tree node.
  The root has a NULL parent.
  The leaves have no children.
  Only the root should be created/destroyed by user.  The creation/destruction of descendents is handled by
  the tree (descendents are deleted when the root is deleted). \n
  Nodes can have a TreeInfo or not.  Recommended is NOT.  Refer to the Node constructor, and to the TreeInfo class for information.
  **/
class Node
{
protected:
    vector<Node*> children;
    string label;
    //string mappingLabel;    //non generic but so important in reconciliation
    Node* parent;

    TreeInfo* treeInfo;

    int depth;
    uint64 pathBits;
    int state;

    unordered_map<string, string> customFields;


    double branchLength;

    Node(TreeInfo* treeInfo = NULL);

    virtual Node* CreateNode(TreeInfo* treeInfo);


    Node* SetAsRootInCopy(Node* ignore);
public:
    /**
      Set maintainTreeInfo = true if you want every node of the tree to
      have a TreeInfo (and hence GetTreeInfo() won't return NULL).\n
      Only one TreeInfo object is maintained and shared by every Node.\n
      Setting to true is not recommended if the tree gets a lot of nodes removed/inserted.
      **/
    Node(bool maintainTreeInfo = false);

    virtual ~Node();


    /**
      Add a child to the current node, and returns the newly created node.
      **/
    Node* AddChild();

    /**
      Inserts a child to the current node at the specified position, and returns the newly created node.
      **/
    Node* InsertChild(int pos);

    int GetNbChildren();
    Node* GetChild(int pos);

    /**
      Get/set Node label, which mainly used to export the tree to Newick.
      **/
    void SetLabel(string lbl);
    string GetLabel();

    //void SetMappingLabel(string lbl);
    //string GetMappingLabel();

    /**
      Returns NULL iif the node is the root
      **/
    Node* GetParent();

    /**
      Get previous child of node's parent, or NULL if node is the first child.\n
      Watch out : this method may need to go through the whole vector of parent's children
      **/
    Node* GetLeftSibling();

    /**
      Get next child of node's parent, or NULL if node is the last child.\n
      Watch out : this method may need to go through the whole vector of parent's children
      **/
    Node* GetRightSibling();


    /**
      Get an iterator that traverses the node and its descendants is post-order.
      Set leavesOnly = true to only iterate through leaves.\n
      WATCH OUT : iterator needs to be closed by calling Node::CloseIterator (or deleting the iterator). \n
      Typical usage :
      @code
      TreeIterator* it = my_node->GetPostOrderIterator();
      while (Node* n = it->next())
      {
        //do fascinating stuff
      }
      my_node->CloseIterator(it);
      @endcode
      **/
    TreeIterator* GetPostOrderIterator(bool leavesOnly = false);




    /**
      Get an iterator that traverses the node and its descendants is pre-order.
      Set leavesOnly = true to only iterate through leaves.
      WATCH OUT : iterator needs to be closed by calling Node::CloseIterator (or deleting the iterator). \n
      Typical usage :
      @code
      TreeIterator* it = my_node->GetPreOrderIterator();
      while (Node* n = it->next())
      {
        //do fascinating stuff
      }
      my_node->CloseIterator(it);
      @endcode
      **/
    TreeIterator* GetPreOrderIterator(bool leavesOnly = false);


    /**
      Closes an iterator previously created by GetPostOrderIterator or GetPreOrderIterator.\n
      Must be called to clear iterator from memory !
      **/
    void CloseIterator(TreeIterator* it);


    /**
      Returns TreeInfo associated with the tree the node belongs to, or NULL if maintainTreeInfo = false
      */
    TreeInfo* GetTreeInfo();


    /**
      Set the distance from the root of the node.  Handled by treeInfo, if there is one.
      */
    void SetDepth(int depth);

    /**
      Get the distance from the root of the node.  Handled by treeInfo, if there is one.  Otherwise, default is zero.
      */
    int GetDepth();

    /**
      Path bits are used to find LCA between nodes in constant time. Handled by treeInfo, if there is one.
      **/
    void SetPathBits(uint64 pathBits);

    /**
      Path bits are used to find LCA between nodes in constant time. Handled by treeInfo, if there is one.
      **/
    uint64 GetPathBits();


    /**
      Returns true iif node is the root (i.e. if parent is NULL)
      */
    bool IsRoot();

    /**
      Returns true iif node is a leaf (i.e. children.size() == 0)
      */
    bool IsLeaf();

    /**
      Set the parent pointer of the node.  Use with caution.
      */
    void SetParent(Node* parent);

    /**
      Add a tree as a new child of the node.
      The added node is the root of the subtree.  It will have its parent
      reaffacted, whether it previously had a parent or not.\n

      WATCH OUT : TreeInfo does not work with this.  AddSubTree can only be used with tree without treeInfo
      **/
    void AddSubTree(Node* node);


    void DeleteTreeInfo();

    /**
      Remove a node that belongs to the children of the node.  This child DOES NOT get deleted, and has its parent set to NULL.
      Since the caller has access to the node, he is expected to delete it.
      */
    void RemoveChild(Node* node);

    /**
      Clear the children vector of the node.  The child Node* objects ARE NOT DELETED.
      */
    void RemoveChildren(bool setParentToNullOnChildren = true);

    /**
      Returns one of STATE_LOSS, STATE_SPECIATION or STATE_DUPLICATION.
      @todo This is non-generic and doesn't belong here.
      */
    int GetState();

    /**
      Sets state to one of STATE_LOSS, STATE_SPECIATION or STATE_DUPLICATION.
      @todo This is non-generic and doesn't belong here.
      */
    void SetState(int state);

    /**
      Special object handled by user.  Gets deleted with the node.
      **/
    NodeInfo* nodeInfo;

    /**
      Copy the properties and nodeInfo of passed node.  Also copies children
      and descendants of the node (not the ancestors), UNLESS the descendants
      is marked as 'not-to-copy' in ignoreNodes.
      **/
    void CopyFrom(Node *n, set<Node*> ignoreNodes = set<Node*>());


    /**
      Takes a function which is called for every node in the tree.
      The function should return true if the received node is to be kept, or false if
      it need to be removed from the tree.
      Removed nodes are deleted from memory.
      arg is a custom argument sent on each call to fncDelete
      **/
    void Restrict(bool (*fncDelete)(Node*, void*), void *arg);


    /**
      This is the not-so-clever non-constant-but-log-time implementation
      **/
    Node* FindLCAWith(Node* n);

    /**
      Calls FindLCAWith as many times as needed.
      **/
    static Node* FindLCA(vector<Node*> nodes);


    /**
      Returns true iif ancestor is on the path between current node and the root (inclusively)
      **/
    bool HasAncestor(Node* ancestor);


    /**
      Creates a node p that becomes the parent of this and sibling, and p becomes the child of the
      previous parent of this and sibling.
      Returns the newly created parent
      NOTE : Does nothing if sibling isn't a sibling.
      NOTE : Do not use with maintainTreeInfo = true
      **/
    Node* InsertParentWith(Node* sibling);

    void SetBranchLength(double length);
    double GetBranchLength();


    Node* GetNodeWithLabel(string label, bool ignoreCase = false);

    set<Node*> GetLeafSet();
    vector<Node*> GetLeafVector();

    set<string> GetLeafLabels();

    Node* GetLeafByLabel(string label);

    /**
     * @brief Returns a copy of the children vector (manipulating this vector will not affect the tree)
     * @return
     */
    vector<Node*> GetChildrenVector();
    /**
      DELETES (clears from memory) the descendants of this node that have only one child
      **/
    void DeleteSingleChildDescendants();


    void BinarizeRandomly();

    /**
     * @brief Returns a list of all the nodes in the tree, ordered by their visiting order in a post-order traversal.
     * The main use is to avoid confusion when the user wants to manipulate the tree while iterating over it.
     * @return
     */
    vector<Node*> GetPostOrderedNodes();


    /**
     * @brief SetAsRootInCopy Copies the current tree and reroots it at current node.
     * This method name is ugly enough, but at least it does exactly what it says.
     * @return The copied + rerooted tree.  Caller MUST delete returned value.
     *
     */
    Node* SetAsRootInCopy();


    Node* SetRootOnParentEdgeInCopy();


    /**
     * @brief GraftOnParentEdge
     * @param nodeToGraft
     */
    Node* GraftOnParentEdge(Node* nodeToGraft);




    void SetCustomField(string name, string val);

    string GetCustomField(string name);


    int GetNbLeaves();

    static void RestrictToLeafset(Node* root, set<Node*> leavesToKeep);
};


bool Node__RestrictToLeafsetFunction(Node* n, void *arg);

#endif // NODE_H
