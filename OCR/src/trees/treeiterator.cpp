#include "treeiterator.h"

TreeIterator::TreeIterator(Node* root)
{
    this->curNode = NULL;
    this->root = root;
    this->leavesOnly = false;
    this->leaves = NULL;
}

void TreeIterator::SetLeavesOnly(bool leavesOnly)
{
    this->leavesOnly = leavesOnly;
}

void TreeIterator::SetLeaves(unordered_set<Node*>* leaves)
{
    this->leaves = leaves;
}


Node* TreeIterator::next()
{
    return NULL; //NOT IMPLEMENTED
}

Node* TreeIterator::DeleteCurrent()
{
    return NULL; //NOT IMPLEMENTED
}


bool TreeIterator::IsLeaf(Node* n)
{
    if (n->IsLeaf())
        return true;
    if (this->leaves && this->leaves->count(n) > 0)
        return true;

    return false;
}






//////////////////////////////////////////////////
//Postorder
//////////////////////////////////////////////////
PostOrderTreeIterator::PostOrderTreeIterator(Node* root) : TreeIterator(root)
{

}

Node* PostOrderTreeIterator::next()
{
    if (curNode == NULL)
    {
        //goto leftmost leaf
        curNode = root;
        while(!this->IsLeaf(curNode))
            curNode = curNode->GetChild(0);

    }
    else if (curNode == root)
    {
        return NULL;
    }
    else
    {
        Node* r = curNode->GetRightSibling();
        if (r == NULL)
        {
            curNode = curNode->GetParent();
        }
        else
        {
            curNode = r;
            while(!this->IsLeaf(curNode))
                curNode = curNode->GetChild(0);
        }

    }

    if (leavesOnly && !this->IsLeaf(curNode))
    {
        return this->next();
    }
    else
    {
        return curNode;
    }
}



Node* PostOrderTreeIterator::DeleteCurrent()
{
    //TODO : THIS PART IS VERY FRAGILE, MEMORY-WISE
    //       Do anything not in this order and a bug/leak is likely to appear

    if (curNode->IsRoot())
    {
        throw "No you can't delete the root like this - just delete it yourself";
        return NULL;
    }

    Node* deadNode = curNode;

    //sets curNode
    this->next();

    Node* pops = deadNode->GetParent();

    pops->RemoveChild(deadNode);

    for (int i = 0; i < deadNode->GetNbChildren(); i++)
    {
        pops->AddSubTree(deadNode->GetChild(i));
    }

    deadNode->RemoveChildren(false);

    delete deadNode;

    return curNode;
}


//////////////////////////////////////////////////
//Preorder
//////////////////////////////////////////////////
PreOrderTreeIterator::PreOrderTreeIterator(Node* root) : TreeIterator(root)
{

}

Node* PreOrderTreeIterator::next()
{
    if (curNode == NULL)
    {
        curNode = root;
    }
    else if (curNode == root)
    {
        if (root->GetNbChildren() == 0)
            return NULL;
        curNode = root->GetChild(0);
    }
    else
    {
        if (!this->IsLeaf(curNode))
        {
            curNode = curNode->GetChild(0);
        }
        else
        {
            if (curNode->GetRightSibling())
            {
                curNode = curNode->GetRightSibling();
            }
            else
            {
                //find first parent with right sibling
                bool ok = true;

                while (ok)
                {
                    if (curNode == root)
                    {
                        curNode = NULL;
                        ok = false;
                    }
                    else
                    {
                        curNode = curNode->GetParent();
                        if (curNode->GetRightSibling())
                        {
                            curNode = curNode->GetRightSibling();
                            ok = false;
                        }
                    }
                }
            }
        }
    }

    if (leavesOnly && !this->IsLeaf(curNode))
    {
        return this->next();
    }
    else
    {
        return curNode;
    }
}





