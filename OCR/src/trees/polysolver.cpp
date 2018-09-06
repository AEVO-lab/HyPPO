#include "polysolver.h"


#include <iostream>

/*ReconcilerV2::ReconcilerV2(Node* geneTree, Node* speciesTree)
{
    this->geneTree = geneTree;
    this->speciesTree = speciesTree;
    this->resolution = NULL;
}*/

PolySolver::~PolySolver()
{
    /*if (resolution)
        delete resolution;*/
    //TODO : anything else ?  Only possibility is the genes resolutions map, but they're in the resolution
}

Node* PolySolver::GetGeneMapping(Node* g)
{
    if (this->genesMapping.find(g) != genesMapping.end())
        return genesMapping[g];

    return NULL;
}





Node *PolySolver::SolvePolytomies(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping)
{
    //nbLosses = 0;
    //nbDups = 0;

    TreeIterator* it = geneTree->GetPostOrderIterator();

    Node* resolution = NULL;

    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            //commented part below tried to map a gene to the one species with the same label
            //This is now the problem of the caller with geneLeavesToSpeciesLeaves
            //Node* m = speciesTree->GetTreeInfo()->GetNodeByLabel(n->GetLabel());

            //TODO : CANNOT MAP LEAF TO SPECIES TREE - WE HAVE A PROBLEM
            //if (!m)
            //    return NULL;

            if (geneLeavesSpeciesMapping.find(n) == geneLeavesSpeciesMapping.end())
            {
                throw "A gene is not mapped to the species tree";
            }

            this->genesMapping[n] = geneLeavesSpeciesMapping[n];
            this->genesResolutions[n] = geneLeavesSpeciesMapping[n];
        }
        else
        {
            vector<Node*> v;

            for (int i = 0; i < n->GetNbChildren(); i++)
            {
                Node* mapped = this->genesMapping[n->GetChild(i)];
                v.push_back(mapped);
            }

            Node* lca = speciesTree->GetTreeInfo()->GetLCA(v);
            this->genesMapping[n] = lca;


            Node* solved = SolvePolytomy(n);
            this->genesResolutions[n] = solved;

            if (n->IsRoot())
                resolution = solved;
        }
    }
    geneTree->CloseIterator(it);


    //cleanup everything for the next call
    genesMapping.clear();
    genesResolutions.clear();

    return resolution;
}



Node* PolySolver::SolvePolytomy(Node* polyroot)
{

    //compute nb(s) for each node of S involded with G
    unordered_map<Node*, PolyInfo*> polyinfos;

    for (int i = 0; i < polyroot->GetNbChildren(); i++)
    {
        Node* s = GetGeneMapping(polyroot->GetChild(i));

        if (polyinfos.count(s) == 0)
        {
            polyinfos[s] = new PolyInfo();
        }

        PolyInfo* pi = polyinfos[s];
        pi->nb++;
        polyinfos[s] = pi;
    }


    ////////////////////////////////////////////////////////////////////////////////
    //compute cup function values after finding the polytomy linked species tree

    LinkedSpeciesTreeUtil sp(polyroot, this);
    TreeIterator* it = sp.GetPostOrderIterator();
    while (Node* s = it->next())
    {
        //found a node for which nb(s) = 0
        if (polyinfos.count(s) == 0)
        {
            polyinfos[s] = new PolyInfo();
        }


        if (sp.IsLeaf(s))
        {
            PolyInfo* pi = polyinfos[s];
            if (pi->nb == 0)
            {
                pi->break1 = 1; pi->break2 = 1; pi->minval = 1;
            }
            else
            {
                pi->break1 = pi->nb;  pi->break2 = pi->nb;  pi->minval = 0;
            }
            polyinfos[s] = pi;
        }
        else
        {
            PolyInfo* curpoly = polyinfos[s];
            PolyInfo* lpoly = polyinfos[s->GetChild(0)];
            PolyInfo* rpoly = polyinfos[s->GetChild(1)];
            int l1 = lpoly->break1;
            int l2 = lpoly->break2;
            int r1 = rpoly->break1;
            int r2 = rpoly->break2;
            ///////////////////////////////////////////////////////
            // Cup function specific stuff - see the paper for an explanation
            if (l2 < r1)
            {
                curpoly->break1 = l2;
                curpoly->break2 = r1;
                curpoly->minval = r1 - l2;
            }
            else if (l1 < r1 && r1 <= l2 && l2 <= r2)
            {
                curpoly->break1 = r1;
                curpoly->break2 = l2;
            }
            else if (l1 < r1 && l2 > r2)
            {
                curpoly->break1 = r1;
                curpoly->break2 = r2;
            }
            else if (r1 <= l1 && l1 <= r2 && r1 <=  l2 && l2 <= r2)
            {
                curpoly->break1 = l1;
                curpoly->break2 = l2;
            }
            else if (r1 <= l1 && l1 <= r2 && l2 > r2)
            {
                curpoly->break1 = l1;
                curpoly->break2 = r2;
            }
            else if (l1 > r2)
            {
                curpoly->break1 = r2;
                curpoly->break2 = l1;
                curpoly->minval = l1 - r2;
            }

            curpoly->minval += (lpoly->minval + rpoly->minval);
            curpoly->break1 += curpoly->nb;
            curpoly->break2 += curpoly->nb;

            ///////////////////////////////////////////////////////


        }
    }
    sp.CloseIterator(it);
    ///////////////////////////////////////////////////////////////////////////////////////


    CalcDupLoss(sp.GetSpeciesTreeRoot(), 1, sp, polyinfos);

    PolyInfo* p = polyinfos[sp.GetSpeciesTreeRoot()];
    //cout<<"Root break1="<<p->break1<<" minval="<<p->minval<<endl;

    Node* ultimateRoot = BuildResolution(polyroot, sp, polyinfos);

    //CLEAR ALL POLYINFO
    unordered_map<Node*, PolyInfo*>::iterator itdel = polyinfos.begin();
    while (itdel != polyinfos.end())
    {
        delete (*itdel).second;
        itdel = polyinfos.erase(itdel);
    }


    return ultimateRoot;
}



Node* PolySolver::BuildResolution(Node* polyroot, LinkedSpeciesTreeUtil &sp, unordered_map<Node*, PolyInfo*> &polyinfos)
{

    //specie forest maintains all subtrees mapped to a given s \in S
    unordered_map<Node*, vector<Node*>*> specieForest;


    //INSERT nb(s) subtrees of each s into the specie forest
    for (int i = 0; i < polyroot->GetNbChildren(); i++)
    {
        Node* child = polyroot->GetChild(i);
        Node* mapped = GetGeneMapping(child);

        if (specieForest.count(mapped) == 0)
        {
            specieForest[mapped] = new vector<Node*>();
        }
        vector<Node*>* curvec = specieForest[mapped];

        Node* n = NULL;
        if (child->IsLeaf())
        {
            n = new Node(false);
            //n->CopyFrom(mapped);
            n->CopyFrom(child);
        }
        else
        {
            n = this->genesResolutions[child];
        }


        curvec->push_back(n);
    }


    //DEBUG
    /*unordered_map<Node*, vector<Node*>*>::iterator itx = specieForest.begin();

    while (itx != specieForest.end())
    {
        vector<Node*>* list = itx->second;
        for (int i = 0; i < list->size(); i++)
        {
            Util::DebugOut(list->at(i)->GetLabel());
        }
        itx++;
    }
    Util::DebugOut("");
    Util::DebugOut("");*/
    //END DEBUG



    TreeIterator* it = sp.GetPostOrderIterator();
    while (Node* s = it->next())
    {

        PolyInfo* curpoly = polyinfos[s];


        if (specieForest.count(s) == 0)
        {
            specieForest[s] = new vector<Node*>();
        }
        vector<Node*>* curvec = specieForest[s];

        //these are the losses
        /*for (int i = 0; i < curpoly->losses; i++)
        {
            Node* n = new Node(false);
            n->SetLabel(s->GetLabel());
            n->SetState(STATE_LOSS);
            nbLosses++;
            curvec->push_back(n);
        }*/

        //these are the duplications
        //up to this point, curvec is the forest of the best k-speciation resolution at s (WTF ???)
        //what's left to do is check if some of these need to be joined by duplications
        //This algo gives a resolution in normal form (see Chang & Eulenstein, 2006)
        //i.e. all dups in n form a caterpillar
        for (int i = 0; i < curpoly->dups; i++)
        {
            Node* n1 = curvec->at(0);
            Node* n2 = curvec->at(1);

            Node* joiner = new Node(false);
            //joiner->SetLabel(n1->GetLabel());
            joiner->SetState(STATE_DUPLICATION);
            //nbDups++;
            joiner->AddSubTree(n1);
            joiner->AddSubTree(n2);

            //the 2 first subtrees are replaced by the newly created root
            curvec->erase(curvec->begin());
            curvec->erase(curvec->begin());

            curvec->insert(curvec->begin(), joiner);
        }

        specieForest[s] = curvec;

        //If we're on the right side of a subtree, the the left side has been taken care of...
        //We just create as many (parent) specs as we can - it there's a disbalance
        //it means there were losses, that didn't get added
        //(because the "these are the losses" part was commented out)
        Node* sibl = s->GetLeftSibling();

        if (sibl && s != sp.GetSpeciesTreeRoot())
        {
            Node* parent = s->GetParent();
            if (specieForest.count(parent) == 0)
            {
                specieForest[parent] = new vector<Node*>();
            }

            vector<Node*>* parentvec = specieForest[parent];
            vector<Node*>* siblvec = specieForest[sibl];

            vector<Node*>* largestVec = curvec;
            vector<Node*>* smallestVec = siblvec;
            if (siblvec->size() > curvec->size())
            {
                largestVec = siblvec;
                smallestVec = curvec;
            }

            int smallestSize = smallestVec->size();
            //this part adds all speciatable stuff we can do
            for (int i = 0; i < smallestSize; i++)
            {
                Node* joiner = new Node(false);
                joiner->SetLabel(parent->GetLabel());

                joiner->AddSubTree(smallestVec->at(0));
                joiner->AddSubTree(largestVec->at(0));

                smallestVec->erase(smallestVec->begin());
                largestVec->erase(largestVec->begin());

                parentvec->push_back(joiner);
            }

            int largestSize = largestVec->size();
            //and this part are stuff we couldn't speciate
            for (int i = 0; i < largestSize; i++)
            {
                parentvec->push_back(largestVec->at(0));
                largestVec->erase(largestVec->begin());
            }

            specieForest[parent] = parentvec;
        }

    }

    sp.CloseIterator(it);

    Node* ultimateRoot = specieForest[sp.GetSpeciesTreeRoot()]->at(0);

    //TODO : DELETE ALL TEMP VECTORS
    /*unordered_map<Node*, vector<Node*>*>::iterator itdel = specieForest.begin();
    while (itdel != specieForest.end())
    {
        delete (*itdel).second;
        itdel = specieForest.erase(itdel);
    }*/

    return ultimateRoot;

}




void PolySolver::CalcDupLoss(Node* s, int k, LinkedSpeciesTreeUtil &sp, unordered_map<Node*, PolyInfo*> &polyinfos)
{
    PolyInfo* curpoly = polyinfos[s];


    if (sp.IsLeaf(s))
    {
        if (k > curpoly->nb)
        {
            curpoly->dups = 0;
            curpoly->losses = k - curpoly->nb;
        }
        else
        {
            curpoly->dups = curpoly->nb - k;
            curpoly->losses = 0;
        }
    }
    else
    {
        int nbs = curpoly->nb;
        PolyInfo* lpoly = polyinfos[s->GetChild(0)];
        PolyInfo* rpoly = polyinfos[s->GetChild(1)];

        //CASE 1 : We can attain k nodes at s by speciating k - nbs children of s, then adding the nbs
        //that were already there
        //Of course k - nbs > 0, otherwise this would imply speciating 0 or less children of s
        //Taking this path implies that we are in he plateau of minimums
        //THIS IS WHERE CHOICES CAN BE MADE : WE MIGHT WANT TO TRAVEL TO A NEIGHBOR VALUE THAT IS ALSO MINIMUM
        if (k - nbs > 0 && curpoly->GetMinCost(k) ==
                            lpoly->GetMinCost(k - nbs) +
                            rpoly->GetMinCost(k - nbs))
        {
            curpoly->dups = 0;
            curpoly->losses = 0;
            CalcDupLoss(s->GetChild(0), k - nbs, sp, polyinfos);
            CalcDupLoss(s->GetChild(1), k - nbs, sp, polyinfos);
        }
        //CASE 2 : we have too many s nodes if we want k -> we'll reduce it with duplications
        //break one is the first value of the plateau
        else if (k < curpoly->break1)
        {
            curpoly->dups = curpoly->break1 - k;
            curpoly->losses = 0;
            CalcDupLoss(s->GetChild(0), curpoly->break1 - nbs, sp, polyinfos);
            CalcDupLoss(s->GetChild(1), curpoly->break1 - nbs, sp, polyinfos);
        }
        //CASE 3 : we need more s nodes -> add losses
        //break two is the last value of the plateau
        else if (k > curpoly->break2)
        {
            curpoly->dups = 0;
            curpoly->losses = k - curpoly->break2;
            CalcDupLoss(s->GetChild(0), curpoly->break2 - nbs, sp, polyinfos);
            CalcDupLoss(s->GetChild(1), curpoly->break2 - nbs, sp, polyinfos);
        }
    }

}


/*int ReconcilerV2::GetNbLosses()
{
    return nbLosses;
}

int ReconcilerV2::GetNbDuplications()
{
    return nbDups;
}*/




//////////////////////////////////////////////////////////////
// LINKEDSPECIESTREEUTIL
//////////////////////////////////////////////////////////////
LinkedSpeciesTreeUtil::LinkedSpeciesTreeUtil(Node* polytomyRoot, PolySolver* rec)
{
    this->polytomyRoot = polytomyRoot;
    this->rec = rec;
    this->speciesTreeRoot = rec->GetGeneMapping(polytomyRoot);
    FindLeaves();
}



Node* LinkedSpeciesTreeUtil::GetSpeciesTreeRoot()
{
    return speciesTreeRoot;
}

void LinkedSpeciesTreeUtil::FindLeaves()
{
    for (int i = 0; i < polytomyRoot->GetNbChildren(); i++)
    {
        Node* g = polytomyRoot->GetChild(i);

        Node* s = rec->GetGeneMapping(g);

        AddLeaf(s);

        bool go = true;
        while (s != this->speciesTreeRoot && go)
        {
            Node* sibl = s->GetRightSibling();
            if (!sibl)
                sibl = s->GetLeftSibling();

            AddLeaf(sibl);

            s = s->GetParent();

            if (!AddInternalNode(s))
            {
                //s was already added...we can stop
                go = false;
            }

        }
    }
}



void LinkedSpeciesTreeUtil::AddLeaf(Node* n)
{
    if (this->traversedNodes.count(n) == 0)
    {
        this->speciesLeaves.insert(n);
    }
}

bool LinkedSpeciesTreeUtil::AddInternalNode(Node* n)
{
    if (this->speciesLeaves.count(n) > 0)
    {
        this->speciesLeaves.erase(n);
    }

    if (this->traversedNodes.count(n) > 0)
    {
        return false;
    }
    else
    {
        traversedNodes.insert(n);
        return true;
    }

}

bool LinkedSpeciesTreeUtil::IsLeaf(Node* n)
{
    return (this->speciesLeaves.count(n) > 0);
}

PostOrderTreeIterator* LinkedSpeciesTreeUtil::GetPostOrderIterator()
{
    PostOrderTreeIterator* p = new PostOrderTreeIterator(this->speciesTreeRoot);
    p->SetLeaves(&speciesLeaves);
    return p;
}

PreOrderTreeIterator* LinkedSpeciesTreeUtil::GetPreOrderIterator()
{
    PreOrderTreeIterator* p = new PreOrderTreeIterator(this->speciesTreeRoot);
    p->SetLeaves(&speciesLeaves);
    return p;
}



void LinkedSpeciesTreeUtil::CloseIterator(TreeIterator* it)
{
    delete it;
}

//////////////////////////////////////////////////////////////
// POLYINFO
//////////////////////////////////////////////////////////////
int PolyInfo::GetMinCost(int k)
{
    if (k < this->break1)
        return this->minval + this->break1 - k;
    else if (k > this->break2)
        return this->minval + k - this->break2;
    else
        return this->minval;
}


string PolyInfo::ToString()
{
    stringstream ss;
    ss<<"NB="<<this->nb<<", B1="<<this->break1<<", B2="<<this->break2<<
        ", MIN="<<this->minval<<", DUPS="<<this->dups<<", LOS="<<this->losses;
    return ss.str();
}





unordered_map<Node*, Node*> PolySolver::GetGeneSpeciesMappingByPrefix(Node* geneTree, Node* speciesTree, string separator)
{
    unordered_map<Node*, Node*> mapping;

    TreeIterator* it = geneTree->GetPostOrderIterator(true);
    while (Node* g = it->next())
    {
        vector<string> sz = Util::Split(g->GetLabel(), separator);

        Node* s = speciesTree->GetTreeInfo()->GetNodeByLabel(sz[0]);
        if (s)
        {
            mapping[g] = s;
        }
    }
    geneTree->CloseIterator(it);

    return mapping;
}




void PolySolver::RestrictGeneTreeByBranchSupport(Node* geneTree, double minThreshold)
{
    void* args = (void*)&minThreshold;
    geneTree->Restrict(&PolySolver__RestrictFunction, args);
}


bool PolySolver__RestrictFunction(Node* n, void* args)
{
    if (n->IsRoot() || n->IsLeaf())
        return true;

    double threshold = *(double*)args;
    double blen = n->GetBranchLength();
    return blen >= threshold;
}
