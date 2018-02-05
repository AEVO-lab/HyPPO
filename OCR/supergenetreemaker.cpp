#include "supergenetreemaker.h"








SuperGeneTreeMaker::SuperGeneTreeMaker()
{

}



int NCALLS = 0;
pair<Node*, int> SuperGeneTreeMaker::GetSuperGeneTreeMinDL(vector<Node *> &trees, vector<unordered_map<Node *, Node *> > &lca_mappings,
                                                Node* speciesTree, bool mustPreserveDupSpec, bool isFirstCall)
{
    if (isFirstCall)
    {

        GeneSpeciesTreeUtil::Instance()->LabelInternalNodesUniquely(trees);
        this->intersectionInfo = TreeLabelIntersectionInfo();
        this->intersectionInfo.ComputeAllIntersections(trees);

        recursionCache.clear();

        //for (int t = 0; t < trees.size(); t++)
        {
            //cout<<"DEBUG INFO!"<<endl<<endl;
            //cout<<NewickLex::ToNewickString(trees[t])<<endl;
            //GeneSpeciesTreeUtil::Instance()->PrintMapping(trees[t], lca_mappings[t]);
            //cout<<endl<<"END DEBUG INFO"<<endl<<endl;
        }

    }


    //check if cached already
    string strtrees = "";
    for (int i = 0; i < trees.size(); i++)
    {
        strtrees += trees[i]->GetLabel() + "-";
    }
    if (recursionCache.find(strtrees) != recursionCache.end())
    {
        pair<Node*, int> sol = recursionCache[strtrees];
        //return a copy: we don't want the cached one to get deleted
        //TODO: eliminate the need for all this copying stuff
        Node* nret = new Node(false);
        nret->CopyFrom(sol.first);
        return make_pair(nret, sol.second);
    }


    /*NCALLS++;

    if (NCALLS % 100 == 0)
    {
        cout<<"NCALLS="<<NCALLS<<" "<<strtrees<<endl;
    }*/

    if (trees.size() == 0)
        throw "No tree was given";

    //I used to put a stop to the recursion in the case of a single tree
    //but not anymore.  The cost of the tree needs to be computed.
    //The recursion does that.
    /*if (trees.size() == 1)
    {
        Node* copy = new Node(false);
        copy->CopyFrom(trees[0]);
        return make_pair(copy, 0);
    }*/

    bool dostop = true;
    string leaflabel = "";
    //the real recursion stop is when there is a single gene
    for (int i = 0; i < trees.size(); i++)
    {
        if (!trees[i]->IsLeaf())
        {
            dostop = false;
            break;  //evil break out of for loop
        }
        else
        {
            if (i == 0)
                leaflabel = trees[i]->GetLabel();
            else if (trees[i]->GetLabel() != leaflabel)
            {
                dostop = false;
                break;
            }

        }
    }

    if (dostop)
    {
        Node* newnode = new Node(false);
        newnode->CopyFrom(trees[0]);    //if we end up here all trees are the same, so take any
        return make_pair(newnode, 0);
    }





    int nbTrees = trees.size();

    //each tree has 4 ways of being split...each way has its index
    //0 = send whole tree left, 1 = send left side left
    //2 = send whole tree right, 3 = send right side left
    //The deal is, if we fix config 0 for the last tree,
    //all 4^{nbtrees - 1}config combos for other trees are non-redundant
    //Then, same if we fix config 1.
    //Afterwards, at config 2 or 3 for the last tree,
    //everything is symmetric to something we've seen before.

    vector<int> counters(nbTrees, 0);

    ApplyNextConfig(counters);   //skip the all-zeros config

    bool done = false;


    Node* currentBestSol = NULL; //we only keep this one
    int currentMinDL = 9999999;  //hopfully enough

    while (!done)
    {
        //evaluate current config corresponding to counters
        vector<Node*> treesLeft;
        vector<Node*> treesRight;
        vector<Node*> speciesLeft;  //should be unordered_set, but meh...
        vector<Node*> speciesRight;


        bool isConfigFine = true;   //some configs are impossible if a tree is a leaf

        //this loop sends what goes left to the left, and what goes right right
        //we also check that no leaf gets split
        for (int c = 0; c < nbTrees; c++)
        {
            Node* tree = trees[c];
            unordered_map<Node*, Node*> lca_mapping = lca_mappings[c];

            //we can't split a leaf
            if (tree->IsLeaf() && (counters[c] == 1 || counters[c] == 3))
            {
                isConfigFine = false;
                break;  //evil break out of for loop
            }

            switch (counters[c])
            {
                case (0):   //all left
                {
                    treesLeft.push_back(tree);
                    speciesLeft.push_back(lca_mapping[tree]);

                    if (!lca_mapping[tree])
                    {
                        cout<<"LCA PROBLEM IN CASE 0"<<endl;
                        throw "LCA mapping problem here";
                    }
                }
                break;
                case(1):    //left goes left
                {
                    Node* t1 = tree->GetChild(0);
                    Node* t2 = tree->GetChild(1);

                    if (!lca_mapping[t2] || !lca_mapping[t1])
                    {
                        cout<<"LCA PROBLEM IN CASE 1"<<endl;
                        cout<<NewickLex::ToNewickString(speciesTree)<<endl;
                        GeneSpeciesTreeUtil::Instance()->PrintMapping(tree, lca_mapping);
                        throw "LCA mapping problem here";
                    }

                    treesLeft.push_back(t1);
                    speciesLeft.push_back(lca_mapping[t1]);
                    treesRight.push_back(t2);
                    speciesRight.push_back(lca_mapping[t2]);

                }
                break;
                case(2):    //all right
                {
                    treesRight.push_back(tree);
                    speciesRight.push_back(lca_mapping[tree]);

                    if (!lca_mapping[tree])
                    {
                        cout<<"LCA PROBLEM IN CASE 2"<<endl;
                        throw "LCA mapping problem here";
                    }
                }
                break;
                case(3):    //right goes left
                {
                    Node* t1 = tree->GetChild(0);
                    Node* t2 = tree->GetChild(1);

                    if (!lca_mapping[t2] || !lca_mapping[t1])
                    {
                        cout<<"LCA PROBLEM IN CASE 3"<<endl;
                        throw "LCA mapping problem here";
                    }

                    treesLeft.push_back(t2);
                    speciesLeft.push_back(lca_mapping[t2]);

                    treesRight.push_back(t1);
                    speciesRight.push_back(lca_mapping[t1]);

                }
                break;
            }
        }

        if (isConfigFine && treesLeft.size() > 0 && treesRight.size() > 0 &&
                !this->intersectionInfo.IsPartitionIntersecting(treesLeft, treesRight))
            /*!this->intersectionInfo.Intersect(treesLeft) && //DOES NOT WORK !  WE MUST CHECK LEFT VS RIGHT
            !this->intersectionInfo.Intersect(treesRight))*/
        {
            //eval losses on created branches + dup, and recurse
            Node* species_lca_left = speciesTree->FindLCA(speciesLeft);
            Node* species_lca_right = speciesTree->FindLCA(speciesRight);
            Node* species_lca = species_lca_left->FindLCAWith(species_lca_right);

            bool isdup = false;
            if (species_lca == species_lca_left || species_lca == species_lca_right)
                isdup = true;

            //verify dup/spec preservation if we need to.
            if (mustPreserveDupSpec)
            {
                //so all trees that got split (counter 1 or 3) must agree with isdup
                for (int tt = 0; tt < trees.size(); tt++)
                {
                    if (counters[tt] == 1 || counters[tt] == 3)
                    {
                        bool tt_isdup = GeneSpeciesTreeUtil::Instance()->IsNodeDup(trees[tt], lca_mappings[tt]);

                        if (tt_isdup != isdup)
                            isConfigFine = false;
                    }
                }

            }

            if (isConfigFine)   //might become false because of check just above
            {
                int losses_left = GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(species_lca_left, species_lca, isdup);

                int losses_right = GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(species_lca_right, species_lca, isdup);

                int localDLCost = losses_left + losses_right + (isdup ? 1 : 0);

                //we have to make sure that trees/mapping indices match.
                //That or make a new tree structure, which I won't
                vector< unordered_map<Node*, Node*> > lca_mappings_left;
                int lcacpt = 0;
                for (int tl = 0; tl < treesLeft.size(); tl++)
                {
                    Node* stmp = NULL;
                    do
                    {
                        stmp = lca_mappings[lcacpt][treesLeft[tl]];
                        lcacpt++;
                    }
                    while (!stmp);
                    lca_mappings_left.push_back(lca_mappings[lcacpt - 1]);
                }

                vector< unordered_map<Node*, Node*> > lca_mappings_right;
                lcacpt = 0;
                for (int tr = 0; tr < treesRight.size(); tr++)
                {
                    Node* stmp = NULL;
                    do
                    {
                        stmp = lca_mappings[lcacpt][treesRight[tr]];
                        lcacpt++;
                    }
                    while (!stmp);
                    lca_mappings_right.push_back(lca_mappings[lcacpt - 1]);
                }




                pair<Node*, int> res_left = GetSuperGeneTreeMinDL(treesLeft, lca_mappings_left, speciesTree, mustPreserveDupSpec, false);
                pair<Node*, int> res_right = GetSuperGeneTreeMinDL(treesRight, lca_mappings_right, speciesTree, mustPreserveDupSpec, false);

                int ttlCost = res_left.second + res_right.second + localDLCost;


                //update min if need be
                if (ttlCost < currentMinDL)
                {
                    if (currentBestSol)
                        delete currentBestSol;
                    currentBestSol = new Node(false);
                    currentBestSol->AddSubTree(res_left.first);
                    currentBestSol->AddSubTree(res_right.first);
                    currentMinDL = ttlCost;
                }
                else
                {
                    delete res_left.first;
                    delete res_right.first;
                }
            }

        }


        ApplyNextConfig(counters);

        if (counters[nbTrees - 1] > 1)   //see long comment above
            done = true;

    }


    //make a copy for the cache - the original might get deleted
    Node* cachedBestSol = new Node(false);
    cachedBestSol->CopyFrom(currentBestSol);
    recursionCache[strtrees] = make_pair(cachedBestSol, currentMinDL);


    if (isFirstCall)
    {
        //cleanup cache
        for (unordered_map<string, pair<Node*, int> >::iterator it = recursionCache.begin(); it != recursionCache.end(); it++)
        {
            Node* cached = (*it).second.first;
            delete cached;
        }
    }

    return make_pair(currentBestSol, currentMinDL);

}



void SuperGeneTreeMaker::ApplyNextConfig(vector<int> &counters)
{
    bool done = false;
    int cindex = 0;

    while (!done)
    {
        if (counters.size() <= cindex)
            done = true;

        counters[cindex] += 1;
        if (counters[cindex] > 3)
        {
            counters[cindex] = 0;
            cindex++;
        }
        else
        {
            done = true;
        }
    }
}











TreeLabelIntersectionInfo::TreeLabelIntersectionInfo()
{
    key_separator = "%*!?!*%"; //hopefully no gene label has this as a substring
}

void TreeLabelIntersectionInfo::ComputeAllIntersections(vector<Node*> trees)
{
    for (int i = 0; i < trees.size(); i++)
    {
        for (int j = i + 1; j < trees.size(); j++)
        {
            ComputeIntersections(trees[i], trees[j]);
        }
    }
}

void TreeLabelIntersectionInfo::ComputeIntersections(Node* tree1, Node* tree2)
{
    TreeIterator* it = tree1->GetPostOrderIterator();

    while (Node* n1 = it->next())
    {
        TreeIterator* it2 = tree2->GetPostOrderIterator();
        while (Node* n2 = it2->next())
        {
            bool intersect = false;
            if (n1->IsLeaf() && n2->IsLeaf())
            {
                if (n1->GetLabel() == n2->GetLabel())
                    intersect = true;
            }
            else if (!n1->IsLeaf())
            {
                Node* n1_1 = n1->GetChild(0);
                Node* n1_2 = n1->GetChild(1);
                if ( Intersect(n1_1->GetLabel(), n2->GetLabel())  ||  Intersect(n1_2->GetLabel(), n2->GetLabel()))
                {
                    intersect = true;
                }
            }
            else   //n1 is leaf, n2 is not
            {
                Node* n2_1 = n2->GetChild(0);
                Node* n2_2 = n2->GetChild(1);
                if ( Intersect(n2_1->GetLabel(), n1->GetLabel())  ||  Intersect(n2_2->GetLabel(), n1->GetLabel()))
                {
                    intersect = true;
                }
            }

            if (intersect)
            {
                AddIntersection(n1->GetLabel(), n2->GetLabel());
            }
        }

        tree2->CloseIterator(it2);
    }

    tree1->CloseIterator(it);
}




void TreeLabelIntersectionInfo::AddIntersection(string lbl1, string lbl2)
{
    intersecting_nodes.insert(lbl1 + key_separator + lbl2);
    intersecting_nodes.insert(lbl2 + key_separator + lbl1);
}

bool TreeLabelIntersectionInfo::Intersect(string lbl1, string lbl2)
{
    string key = lbl1 + key_separator + lbl2;
    return (intersecting_nodes.find(key) != intersecting_nodes.end());
}


bool TreeLabelIntersectionInfo::Intersect(vector<Node*> trees)
{
    for (int i = 0; i < trees.size(); i++)
    {
        for (int j = i + 1; j < trees.size(); j++)
        {
            if (Intersect(trees[i]->GetLabel(), trees[j]->GetLabel()))
            {
                return true;
            }
        }
    }

    return false;
}


bool TreeLabelIntersectionInfo::IsPartitionIntersecting(vector<Node*> treesLeft, vector<Node*> treesRight)
{
    for (int i = 0; i < treesLeft.size(); i++)
    {
        for (int j = 0; j < treesRight.size(); j++)
        {
            if (Intersect(treesLeft[i]->GetLabel(), treesRight[j]->GetLabel()))
                return true;
        }
    }

    return false;
}
