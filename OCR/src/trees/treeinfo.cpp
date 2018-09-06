#include "treeinfo.h"

TreeInfo::TreeInfo(Node* root)
{
    this->root = root;
    this->degree = 1;
}


void TreeInfo::OnLabelChanged(Node* n, string oldLabel, string newLabel)
{
    if (this->labelNodes.find(oldLabel) != labelNodes.end())
    {
        labelNodes.erase(oldLabel);
    }

    if (n->GetLabel() != "")
    {
        this->labelNodes[n->GetLabel()] = n;
        this->labelNodes[Util::ToUpper(n->GetLabel())] = n;
    }


}


void TreeInfo::OnNodeInserted(Node* n, int pos)
{
    if (!n->IsRoot())
    {
        if (n->GetParent()->GetNbChildren() > degree)
            degree = n->GetParent()->GetNbChildren();
    }

    /*if (this->nameEmptyLabels && n->GetLabel() == "")
    {
        GiveNodeSomeLabel(n);
    }

    if (n->GetLabel() != "")
    {
        this->labelNodes[n->GetLabel()] = n;
    }


    if (this->retainDepth)
    {
        if (n->IsRoot())
            n->SetDepth(0);
        else
            n->SetDepth(n->GetParent()->GetDepth() + 1);
    }


    //TODO : We assume the tree is binary
    if (this->retainBinaryLCAPath)
    {
        uint64 pathbits = 0;
        uint64 buf = 1;

        //The path gets the parent's path, with the depth bit changed to 1
        //and the depth - 1 bit set to 0 if it is the leftmost node
        if (n->IsRoot())
        {
            buf = (buf << (uint64)63);
            pathbits |= buf;
        }
        else
        {
            pathbits = n->GetParent()->GetPathBits();
            pathbits |= (buf << (uint64)(63 - n->GetDepth()));

            if (pos == 0)
                pathbits &= ~(buf << (uint64)(63 - n->GetDepth() + 1));
        }

        n->SetPathBits(pathbits);


    }*/

}



void TreeInfo::GiveNodeSomeLabel(Node* n, int labelCounter)
{
    stringstream ss;
    ss<<"N"<<labelCounter;
    n->SetLabel(ss.str());

}

Node* TreeInfo::GetRoot()
{
    return root;
}


int TreeInfo::GetDegree()
{
    return degree;
}


void TreeInfo::OnNodeDeleted()
{
    ParseTree(this->GetRoot(), false, true, true);
}


void TreeInfo::ParseTree(Node* node, bool nameEmptyLabels, bool computeDepth, bool computePathBits)
{
    if (!node)
        node = this->root;

    computeDepth = computeDepth || computePathBits;

    if (computePathBits)
        pathBitNodes.clear();

    TreeIterator* it = node->GetPreOrderIterator();

    int cpt = 1;
    while (Node* n = it->next())
    {
        if (nameEmptyLabels && n->GetLabel() == "")
        {
            GiveNodeSomeLabel(n, cpt);
            cpt++;
        }

        if (computeDepth)
        {
            if (n->IsRoot())
                n->SetDepth(0);
            else
                n->SetDepth(n->GetParent()->GetDepth() + 1);
        }

        if (computePathBits && degree == 2)
        {
            uint64 pathbits = 0;
            uint64 buf = 1;

            //The path gets the parent's path, with the depth bit changed to 1
            //and the depth - 1 bit set to 0 if it is the leftmost node
            if (n->IsRoot())
            {
                buf = (buf << (uint64)63);
                pathbits |= buf;
            }
            else
            {
                pathbits = n->GetParent()->GetPathBits();
                pathbits |= (buf << (uint64)(63 - n->GetDepth()));

                if (n->GetRightSibling())
                    pathbits &= ~(buf << (uint64)(63 - n->GetDepth() + 1));
            }

            n->SetPathBits(pathbits);
            pathBitNodes[pathbits] = n;
        }


    }

    node->CloseIterator(it);
}



Node* TreeInfo::GetNodeByPathBits(uint64 pathbits)
{
    if (pathBitNodes.find(pathbits) != pathBitNodes.end())
        return pathBitNodes[pathbits];
    else
        return NULL;
}

Node* TreeInfo::GetNodeByLabel(string label)
{
    if (labelNodes.find(label) != labelNodes.end())
        return labelNodes[label];
    else
    {
        if (labelNodes.find(Util::ToUpper(label)) != labelNodes.end())
            return labelNodes[Util::ToUpper(label)];
        else
        {
            if (labelNodes.find(Util::ToLower(label)) != labelNodes.end())
                return labelNodes[Util::ToLower(label)];
            else
                return NULL;
        }
    }
}




Node* TreeInfo::GetLCA(Node* n1, Node* n2)
{
    if (n1->GetPathBits() == 0)
        this->ParseTree(NULL, false, true, true);

    uint64 path1 = n1->GetPathBits();
    uint64 cmp = path1 ^ n2->GetPathBits();

    //same path => same node
    if (cmp == 0)
        return n1;

    uint64 res = 0;

    bool agreeing = true;
    uint64 cpt = 63;

    while (agreeing)
    {
        uint64 buf = 1;
        buf = buf << cpt;

        if ((cmp & buf) == 0)
        {
            res |= (path1 & buf);

            cpt--;
        }
        else
        {
            res |= buf;
            agreeing = false;
        }
    }

    //check if disagreeing bit was a path bit or not
    cpt = 64 - cpt - 1;

    //if not, one is ancestor of the other.  Return the highest node
    if (cpt > n1->GetDepth() || cpt > n2->GetDepth())
    {
        if (n1->GetDepth() < n2->GetDepth())
            return n1;
        else
            return n2;
    }
    else
    {
        return pathBitNodes[res];
    }

}

Node* TreeInfo::GetLCA(vector<Node*> &nodes)
{
    if (nodes.size() == 0)
        return NULL;
    if (nodes.size() == 1)
        return nodes[0];

    Node* cur = GetLCA(nodes[0], nodes[1]);
    for (int i = 2; i < nodes.size(); i++)
    {
        cur = GetLCA(cur, nodes[i]);
    }

    return cur;
}

