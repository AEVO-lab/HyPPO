#include "polysolver_distance.h"






void PSDCacheManager::CacheResolution(Node* orig, Node* resolved, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping)
{
    string key = GetKey(orig);

    if (!GetCachedResolution(key))
    {
        PSDCachedResolution* p = new PSDCachedResolution();
        p->speciesTree = speciesTree;

        p->resolution = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(orig, lcaMapping, p->lcaMapping);
        cache[key] = p;
    }
}

string PSDCacheManager::GetKey(Node* orig)
{
    vector<string> leafLabels;

    TreeIterator* it = orig->GetPostOrderIterator(true);
    while (Node* n = it->next())
    {
        leafLabels.push_back(n->GetLabel());
    }
    orig->CloseIterator(it);

    std::sort(leafLabels.begin(), leafLabels.end());

    string key = "";

    for (int i = 0; i < leafLabels.size(); i++)
    {
        if (key != "")
            key += ",";
        key += leafLabels[i];
    }

    return key;
}

PSDCachedResolution* PSDCacheManager::GetCachedResolution(Node* orig, Node* speciesTree)
{
    string key = GetKey(orig);
    return GetCachedResolution(key);
}

PSDCachedResolution* PSDCacheManager::GetCachedResolution(string key)
{
    if (cache.find(key) == cache.end())
        return NULL;
    return cache[key];
}

void PSDCacheManager::CleanUp()
{
    for (unordered_map<string, PSDCachedResolution*>::iterator it = cache.begin(); it != cache.end(); it++)
    {
        PSDCachedResolution* p = (*it).second;
        delete p;
    }
    cache.clear();
}


void PSDNodeInfo::CopyFrom(NodeInfo* n)
{
    this->psdName = ((PSDNodeInfo*)n)->psdName;
    this->totalBranchLength = ((PSDNodeInfo*)n)->totalBranchLength;
    this->isLoss = ((PSDNodeInfo*)n)->isLoss;
}

NodeInfo* PSDNodeInfo::GetClone()
{
    PSDNodeInfo* p = new PSDNodeInfo();
    p->CopyFrom(this);
    return p;
}



Node* PolySolverDistance::SolvePolytomies(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances)
{
    unordered_map<Node*, Node*> mylcaMapping;

    //we work on a copy of the gene tree not to destroy it...could this be prevented in order to increase speed ?
    Node* curResolution = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(geneTree, lcaMapping, mylcaMapping);

    //This first pass here finds every highest subtree that was previously resolved and stored in the cache.
    //It's a very naive O(n^2) pass but it's certainly worth it !
    if (useCache)
    {
        TreeIterator* it = curResolution->GetPreOrderIterator();
        vector<pair<Node*, PSDCachedResolution*> > replacements;
        while (Node* n = it->next())
        {
            if (!n->IsLeaf()  && !n->IsRoot())
            {
                bool hasResolvedAncestor = false;
                //verify that n's ancestor is not going to get replaced
                for (int i = 0; i < replacements.size() && !hasResolvedAncestor; i++)
                {
                    Node* p = replacements[i].first;
                    if (n->HasAncestor(p))
                        hasResolvedAncestor = true;
                }

                if (!hasResolvedAncestor)
                {
                    PSDCachedResolution* cached = PSDCacheManager::Instance()->GetCachedResolution(n, speciesTree);
                    if (cached)
                    {
                        replacements.push_back(make_pair(n, cached));
                    }
                }
            }
        }
        curResolution->CloseIterator(it);


        for (int i = 0; i < replacements.size(); i++)
        {

            Node* n = replacements[i].first;
            PSDCachedResolution* p = replacements[i].second;



            Node* res = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(p->resolution, p->lcaMapping, mylcaMapping);

            if (verbose > 1)
            {
                cout<<"Used cache on "<<PSDCacheManager::Instance()->GetKey(n)<<" subtree"<<endl;
                cout<<"Replaced by "<<NewickLex::ToNewickString(p->resolution)<<endl;
                /*vector<Node*> v = res->GetPostOrderedNodes();

                for (int k = 0; k < v.size(); k++)
                {
                    cout<<"lca["<<v[k]->GetLabel()<<"] = "<<mylcaMapping[v[k]]->GetLabel()<<endl;
                }*/
            }

            //note : n cannot be the root
            Node* parent = n->GetParent();
            parent->RemoveChild(n);
            parent->AddSubTree(res);
            delete n;

        }
    }




    //STARTRESOLVING
    vector<Node*> geneNodes = curResolution->GetPostOrderedNodes();
    for (int i = 0; i < geneNodes.size(); i++)
    {
        Node* n = geneNodes[i];
        if (n->GetNbChildren() > 2)
        {
            vector<Node*> children = n->GetChildrenVector();
            Node* res = this->SolvePolytomy(children, speciesTree, mylcaMapping, distances);



            //replace n by res
            if (n->IsRoot())
            {
                delete curResolution;
                curResolution = res;
            }
            else
            {
                if (useCache)
                    PSDCacheManager::Instance()->CacheResolution(n, res, speciesTree, mylcaMapping);
                Node* parent = n->GetParent();
                parent->RemoveChild(n);
                parent->AddSubTree(res);
                delete n;
            }
        }
    }


    return curResolution;


}


Node* PolySolverDistance::SolvePolytomy(vector<Node*> leaves, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances)
{
    //first, cut off unnecessary ancestor in the species tree
    vector<Node*> relevantSpecies = GeneSpeciesTreeUtil::Instance()->GetGenesSpecies(leaves, lcaMapping);
    Node* highestSpecies = speciesTree->FindLCA(relevantSpecies);
    speciesTree = highestSpecies;

    set<Node*> speciesLeaves = this->ComputeLinkSpeciesLeaves(leaves, lcaMapping, speciesTree);

    if (!speciesTree->IsRoot() && verbose > 0)
    {
        cout<<"RELOCATED SPECIES ROOT"<<endl;
    }

    //maxNbS is the count of genes in the species S that has the most genes
    int maxNbS = 0;

    int nbSpeciesDone = 0;

    //first key : species, second : k
    //M_{s, k} is cells[s][k]
    map<Node*, map<int, PSDCell*> > cells;

    //This one contains, for each species s, the list of genes mapped to s
    map<Node*, vector<Node*> > speciesGenes;
    vector<Node*>::iterator lfit = leaves.begin();
    while (lfit != leaves.end())
    {
        Node* leaf = (*lfit);
        Node* sp = lcaMapping[leaf];
        if (speciesGenes.find(sp) == speciesGenes.end())
        {
            speciesGenes[sp] = vector<Node*>();
        }

        vector<Node*> v = speciesGenes[sp];
        v.push_back(leaf);
        speciesGenes[sp] = v;

        if (v.size() > maxNbS)
            maxNbS = v.size();

        lfit++;
    }


    //and here we go...we construct the table here
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* s = it->next())
    {
        //special case : s is uninteresting because it's below a linked species tree leaf
        //TODO : this takes useless time
        bool interesting = true;
        for (set<Node*>::iterator it = speciesLeaves.begin(); it != speciesLeaves.end(); it++)
        {
            Node* sleaf = (*it);
            if (s != sleaf && s->HasAncestor(sleaf))
            {
                interesting = false;
            }
        }
        if (!interesting)
        {
            continue;
        }

        int nbS = speciesGenes[s].size();

        int leftMinK = 0;
        int rightMinK = 0;
        int minC = 99999;

        //Compute all C_{s,k}
        for (int k = 1; k <= maxNbS * 2; k++)
        {
            PSDCell* cell = new PSDCell();

            //s isn't a leaf (of speciesLeaves) case => value of cell is based on 2 upper rows
            if (speciesLeaves.find(s) == speciesLeaves.end() && cells.find(s->GetChild(0)) != cells.end() && cells.find(s->GetChild(1)) != cells.end())
            {
                if (k - nbS > 0)
                {
                    PSDCell* upcell1 = cells[s->GetChild(0)][k - nbS];
                    PSDCell* upcell2 = cells[s->GetChild(1)][k - nbS];
                    cell->CVal = upcell1->MVal + upcell2->MVal;

                    //The C values form a plateau of minimums.  This remembers at which k the plateau starts and stops
                    if (cell->CVal < minC)
                    {
                        minC = cell->CVal;
                        leftMinK = k;
                    }
                    if (leftMinK != 0 && cell->CVal == minC)
                    {
                        rightMinK = k;
                    }
                }
                else
                {
                    cell->CVal = 99999;
                }
            }
            //s is a leaf case
            else
            {
                if (k == nbS)
                {
                    cell->CVal = 0;
                    leftMinK = k;
                    rightMinK = k;
                    minC = 0;
                }
                else
                    cell->CVal = 99999;
            }

            cells[s][k] = cell;

        }

        //NOTE that at this point, leftMinK and rightMinK and minC might still be their default value.
        //This is the case when s is a leaf, and nbS is zero.  Then, every cell is the product of a loss addition.
        //We call it the ALL-IS-LOSS-SPECIAL-CASE


        //Then, compute all M_{s, k}
        //All this magic is in the papeeeeeer !  M is a plateau function will the same mins as C

        //What happens is that we start creating the subtrees for the cells where C_{s,k} is at the plateau.
        //These cells can only correspond to multiple speciations, and we can construct them unambiguously.
        //Then, we build the subtrees from leftMinK to 1 (these cells NEED the subtrees of C_{s,leftMinK})
        if (leftMinK > 0 && rightMinK > 0)
        {
            for (int k = leftMinK; k <= rightMinK; k++)
            {
                PSDCell* cell = cells[s][k];
                cell->MVal = minC;
                vector<Node*> specSubtrees = this->MakeSpecSubtrees(s, k, speciesLeaves, speciesGenes, cells, lcaMapping, distances);
                cell->spec_subtrees = specSubtrees;

                //debug
                /*for (int j = 0; j < cell->spec_subtrees.size(); j++)
                {
                    Node* nn = cell->spec_subtrees[j];
                    cout<<"s="<<s->GetLabel()<<" "<<nn->GetLabel()<<endl;
                }*/
            }
        }

        for (int k = leftMinK - 1; k >= 1; k--)
        {
            //BUILD, ACCORDING TO cells[s, k + 1]
            PSDCell* cell = cells[s][k];
            cell->MVal = minC + (leftMinK - k);
            vector<Node*> dupSubtrees = this->MakeOneDuplication(s, k, speciesLeaves, speciesGenes, cells, lcaMapping, distances);
            cell->duploss_subtrees = dupSubtrees;


            //special case when making the dup/loss has same cost as making the specs
            if (cell->MVal == cell->CVal)
            {
                vector<Node*> specSubtrees = this->MakeSpecSubtrees(s, k, speciesLeaves, speciesGenes, cells, lcaMapping, distances);
                cell->spec_subtrees = specSubtrees;
            }


            /*vector<Node*> trees = cell->GetDefaultSubtrees();
            cout<<"CELL["<<s->GetLabel()<<","<<k<<"] : CVal="<<cell->CVal<<" MVal="<<cell->MVal<<" DEFTREES="<<trees.size()<<
                  " STREES="<<cell->spec_subtrees.size()<<" DTREES="<<cell->duploss_subtrees.size()<<endl;*/
        }


        for (int k = rightMinK + 1; k <= maxNbS * 2; k++)
        {
            PSDCell* cell = cells[s][k];

            //ALL-IS-LOSS-SPECIAL-CASE : on the first pass, we make as if the smallest value was C[s, 0] = 0
            if (k == 1)
            {
                minC = 0;
                rightMinK = 0;  //should already be zero, but here for clarification
            }


            cell->MVal = minC + (k - rightMinK);

            vector<Node*> lossSubtrees = this->MakeOneLoss(s, k, speciesGenes, cells, lcaMapping, distances);
            cell->duploss_subtrees = lossSubtrees;


            //special case when making the dup/loss has same cost as making the specs
            if (cell->MVal == cell->CVal)
            {
                vector<Node*> specSubtrees = this->MakeSpecSubtrees(s, k, speciesLeaves, speciesGenes, cells, lcaMapping, distances);
                cell->spec_subtrees = specSubtrees;
            }


            /*vector<Node*> trees = cell->GetDefaultSubtrees();
            cout<<"CELL["<<s->GetLabel()<<","<<k<<"] : CVal="<<cell->CVal<<" MVal="<<cell->MVal<<" DEFTREES="<<trees.size()<<
                  " STREES="<<cell->spec_subtrees.size()<<" DTREES="<<cell->duploss_subtrees.size()<<endl;*/

        }

        //cout<<"LK="<<leftMinK<<" RK="<<rightMinK<<endl;




        nbSpeciesDone++;

        if (nbSpeciesDone % 10 == 0 && verbose > 0)
        {
            cout<<nbSpeciesDone<<" species done"<<endl;
        }

    }
    speciesTree->CloseIterator(it);

    //PrintScoresTable(speciesTree, leaves.size(), cells);

    //copy the solution at M[root of S, 0], and delete every cell and temp subtree
    vector<Node*> rootTrees = cells[speciesTree][1]->GetDefaultSubtrees();

    if (verbose > 0)
        cout<<"score="<<cells[speciesTree][1]->MVal<<endl;


    Node* retTree = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(rootTrees[0], lcaMapping, lcaMapping);

    if (!retTree->nodeInfo)
    {
        retTree->nodeInfo = new PSDNodeInfo();
    }
    ((PSDNodeInfo*)retTree->nodeInfo)->dlScore = cells[speciesTree][1]->MVal;

    //clean up everything
    for (map<Node*, map<int, PSDCell*> >::iterator cellit = cells.begin(); cellit != cells.end(); cellit++)
    {
        map<int, PSDCell*> cellmap = (*cellit).second;

        for (map<int, PSDCell*>::iterator cellmapit = cellmap.begin(); cellmapit != cellmap.end(); cellmapit++)
        {
            PSDCell* cell = (*cellmapit).second;
            delete cell;
        }
        cellmap.clear();
    }
    cells.clear();

    return retTree;
}


vector<Node*> PolySolverDistance::MakeSpecSubtrees(Node* s, int k, set<Node*> speciesLeaves, map<Node*, vector<Node*> > &speciesGenes, map<Node*, map<int, PSDCell*> > &cells, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances)
{
    vector<Node*> builtSubtrees;


    if (verbose > 1)
    {
        vector<Node*> s_leaves = s->GetLeafVector();
        string str_leaves = "";
        for (int i = 0; i < s_leaves.size(); i++)
        {
            if (str_leaves != "")
                str_leaves += ",";
            str_leaves += s_leaves[i]->GetLabel();
        }
        cout<<"------------------------------------------"<<endl;
        cout<<"M[s, "<<k<<"] : making k speciations, s = ("<<str_leaves<<")"<<endl;

        if (speciesLeaves.find(s) != speciesLeaves.end())
            cout<<"(nothing to do, s is a leaf)"<<endl;
    }


    //here we copy every gene subtree mapped to s that were given to us
    vector<Node*> v = speciesGenes[s];
    for (int i = 0; i < v.size(); i++)
    {
        Node* leaf = new Node(false);
        leaf->CopyFrom(v[i]);
        builtSubtrees.push_back(leaf);
        lcaMapping[leaf] = lcaMapping[v[i]];
    }


    //and here, we make spec subtrees by joining subtrees mapped to children of s
    //The case where cells[s->GetChild(0)] is not there is when s is a leaf of the
    //considered species tree
    if (speciesLeaves.find(s) == speciesLeaves.end() && cells.find(s->GetChild(0)) != cells.end() && cells.find(s->GetChild(1)) != cells.end())
    {



        //how does that work ?
        //We need to make k spec subtrees with the root mapped to s
        //and we already have nbS of them in speciesGenes - thus we
        //join k - nbS left/right children of s
        int nbS = speciesGenes[s].size();
        PSDCell* c1 = cells[s->GetChild(0)][k - nbS];
        PSDCell* c2 = cells[s->GetChild(1)][k - nbS];

        vector<Node*> c1trees = c1->GetDefaultSubtrees();
        vector<Node*> c2trees = c2->GetDefaultSubtrees();

        //copy these trees in a set...we'll gradually join the trees and remove them from this set,
        //so in the end, it's faster to just make a set copy
        set<Node*> c1trees_set = set<Node*>(c1trees.begin(), c1trees.end());
        set<Node*> c2trees_set = set<Node*>(c2trees.begin(), c2trees.end());

        vector<Node*> alltrees = c1trees;
        for (int i = 0; i < c2trees.size(); i++)
        {
            alltrees.push_back(c2trees[i]);
        }
        map<Node*, map<Node*, double> > treesDistances = this->GetSubtreesDistances(alltrees, distances);


        // VERBOSE OUTPUT -------------------------------------------
        map<Node*, string> tempTreeNames;
        if (verbose > 1)
        {
            cout<<"Candidate gene trees : "<<endl;

            int cnt = 0;
            for (map<Node*, map<Node*, double> >::iterator vit = treesDistances.begin(); vit != treesDistances.end(); vit++)
            {
                string tname = "G" + Util::ToString(cnt);
                cout<<tname<<" = "<<NewickLex::ToNewickString((*vit).first)<<endl;
                tempTreeNames[(*vit).first] = tname;

                cnt++;
            }

            cout<<"Distances : "<<endl;
            for (int i = 0; i < c1trees.size(); i++)
            {
                for (int j = 0; j < c2trees.size(); j++)
                {
                    cout<<tempTreeNames[c1trees[i]]<<"  "<<tempTreeNames[c2trees[j]]<<"  "<<treesDistances[c1trees[i]][c2trees[j]]<<endl;
                }
            }
        }
        // /VERBOSE OUTPUT -------------------------------------------


        int cptJoin = 1;




        while (c1trees_set.size() > 0 && c2trees_set.size() > 0)
        //for (int i = 1; i <= k - nbS; i++)
        {
            //find the two trees with min distance
            //TODO : this is clearly suboptimal - we should sort all pairs only once
            double minDist = 9999999;
            Node* t1min = NULL;
            Node* t2min = NULL;

            //Here we use c1it (that is, c-one-i-t, and not the vulgar CLIT as some may think)
            for (set<Node*>::iterator c1it = c1trees_set.begin(); c1it != c1trees_set.end(); c1it++)
            {
                Node* t1 = (*c1it);
                for (set<Node*>::iterator c2it = c2trees_set.begin(); c2it != c2trees_set.end(); c2it++)
                {
                    Node* t2 = (*c2it);




                    if (treesDistances[t1][t2] < minDist)
                    {
                        minDist = treesDistances[t1][t2];
                        t1min = t1;
                        t2min = t2;
                    }
                }
            }

            c1trees_set.erase(t1min);
            c2trees_set.erase(t2min);


            Node* t1 = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(t1min, lcaMapping, lcaMapping);
            Node* t2 = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(t2min, lcaMapping, lcaMapping);

            Node* newspec = new Node(false);
            newspec->AddSubTree(t1);
            newspec->AddSubTree(t2);
            builtSubtrees.push_back(newspec);
            lcaMapping[newspec] = s;

            if (verbose > 1)
            {
                cout<<"JOIN "<<cptJoin<<" : "<<tempTreeNames[t1min]<<"   "<<tempTreeNames[t2min]<<endl;
                cptJoin++;
            }
        }


        //it's possible there aren't as many left trees as right trees (eg one was added as a loss - so it's not really there)
        //In this case, we copy leftover trees in this cell
        set<Node*>::iterator c1trees_it = c1trees_set.begin();
        while (c1trees_it != c1trees_set.end())
        {
            Node* toCopy = (*c1trees_it);
            Node* t1 = new Node(false);
            t1->CopyFrom(toCopy);
            builtSubtrees.push_back(t1);
            lcaMapping[t1] = lcaMapping[toCopy];
            //c1trees_it = c1trees_set.erase(c1trees_it);
            c1trees_set.erase(c1trees_it++);
        }

        set<Node*>::iterator c2trees_it = c2trees_set.begin();
        while (c2trees_it != c2trees_set.end())
        {
            Node* toCopy = (*c2trees_it);
            Node* t2 = new Node(false);
            t2->CopyFrom(toCopy);
            builtSubtrees.push_back(t2);
            lcaMapping[t2] = lcaMapping[toCopy];
            //c2trees_it = c2trees_set.erase(c2trees_it);
            c2trees_set.erase(c2trees_it++);
        }

    }

    if (verbose > 1)
        cout<<"------------------------------------------"<<endl;

    return builtSubtrees;
}



vector<Node*> PolySolverDistance::MakeOneDuplication(Node* s, int k, set<Node*> speciesLeaves, map<Node*, vector<Node*> > &speciesGenes, map<Node*, map<int, PSDCell*> > &cells, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances)
{

    if (verbose > 1)
    {
        vector<Node*> s_leaves = s->GetLeafVector();
        string str_leaves = "";
        for (int i = 0; i < s_leaves.size(); i++)
        {
            if (str_leaves != "")
                str_leaves += ",";
            str_leaves += s_leaves[i]->GetLabel();
        }
        cout<<"------------------------------------------"<<endl;
        cout<<"M[s, "<<k<<"] : one dup required, s = ("<<str_leaves<<")"<<endl;
    }

    vector<Node*> builtSubtrees;

    PSDCell* cellright = cells[s][k + 1];
    vector<Node*> righttrees = cellright->GetDefaultSubtrees();


    map<Node*, map<Node*, double> > treesDistances = this->GetSubtreesDistances(righttrees, distances);




    //find the 2 trees with min distance
    double minDist = 9999999;
    Node* t1min = NULL;
    Node* t2min = NULL;
    for (int i = 0; i < righttrees.size(); i++)
    {
        for (int j = i + 1; j < righttrees.size(); j++)
        {
            double sumt1 = 0;
            double sumt2 = 0;

            for (int k = 0; k < righttrees.size(); k++)
            {
                sumt1 += treesDistances[righttrees[i]][righttrees[k]];
                sumt2 += treesDistances[righttrees[j]][righttrees[k]];
            }

            double dd = (double)(righttrees.size() - 2) * treesDistances[righttrees[i]][righttrees[j]] - sumt1 - sumt2;

            //if (treesDistances[righttrees[i]][righttrees[j]] < minDist)
            if (dd < minDist)
            {
                //minDist = treesDistances[righttrees[i]][righttrees[j]];
                minDist = dd;
                t1min = righttrees[i];
                t2min = righttrees[j];
            }
        }
    }


    // VERBOSE OUTPUT -------------------------------------------
    if (verbose > 1)
    {

        cout<<"Candidate gene trees : "<<endl;

        map<Node*, string> tempTreeNames;

        int cnt = 0;
        for (map<Node*, map<Node*, double> >::iterator vit = treesDistances.begin(); vit != treesDistances.end(); vit++)
        {
            string tname = "G" + Util::ToString(cnt);
            if (speciesLeaves.find(s) != speciesLeaves.end())
            {
                tname = NewickLex::ToNewickString((*vit).first);
                cout<<tname<<endl;
            }
            else
            {

                cout<<tname<<" = "<<NewickLex::ToNewickString((*vit).first)<<endl;
            }
            tempTreeNames[(*vit).first] = tname;

            cnt++;
        }

        cout<<"Distances : "<<endl;
        for (int i = 0; i < righttrees.size(); i++)
        {
            for (int j = i + 1; j < righttrees.size(); j++)
            {
                cout<<tempTreeNames[righttrees[i]]<<"  "<<tempTreeNames[righttrees[j]]<<"  "<<treesDistances[righttrees[i]][righttrees[j]]<<endl;
            }
        }

        cout<<"Chosen join : "<<tempTreeNames[t1min]<<"  "<<tempTreeNames[t2min]<<endl;
    }
    // /VERBOSE OUTPUT -------------------------------------------


    if (t1min && t2min)
    {
        Node* t1 = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(t1min, lcaMapping, lcaMapping);
        Node* t2 = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(t2min, lcaMapping, lcaMapping);

        Node* newdup = new Node(false);
        newdup->AddSubTree(t1);
        newdup->AddSubTree(t2);
        builtSubtrees.push_back(newdup);
        lcaMapping[newdup] = s;
    }

    for (int i = 0; i < righttrees.size(); i++)
    {
        Node* t = righttrees[i];
        if (t != t1min && t != t2min)
        {
            Node* copied = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(righttrees[i], lcaMapping, lcaMapping);
            builtSubtrees.push_back(copied);
        }
    }


    if (verbose > 1)
        cout<<"------------------------------------------";

    return builtSubtrees;
}


vector<Node*> PolySolverDistance::MakeOneLoss(Node* s, int k, map<Node*, vector<Node*> > &speciesGenes, map<Node*, map<int, PSDCell*> > &cells, unordered_map<Node*, Node*> &lcaMapping, map<Node*, map<Node*, double> > &distances)
{
    //TODO : this just blindly adds a loss node in s
    //TODO : NO ACTUALLY IT DOESN'T CREATE A NODE FOR THE LOSS - We just take the same trees as on the left cells
    vector<Node*> builtSubtrees;

    if (k > 1)
    {
        PSDCell* cellleft = cells[s][k - 1];
        vector<Node*> lefttrees = cellleft->GetDefaultSubtrees();

        for (int i = 0; i < lefttrees.size(); i++)
        {
            Node* copied = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(lefttrees[i], lcaMapping, lcaMapping);
            builtSubtrees.push_back(copied);
        }
    }


    /*Node* newloss = new Node(false);
    PSDNodeInfo* info = new PSDNodeInfo();
    info->isLoss = true;
    newloss->nodeInfo = info;
    builtSubtrees.push_back(newloss);
    lcaMapping[newloss] = s;*/

    return builtSubtrees;
}


map<Node*, map<Node*, double> > PolySolverDistance::GetSubtreesDistances(vector<Node*> subtrees, map<Node*, map<Node*, double> > &distances)
{
    //bool allHomo = true;
    map<string, Node*> allMyLeaves;
    for (int i = 0; i < subtrees.size(); i++)
    {
        vector<Node*> tleaves = subtrees[i]->GetLeafVector();
        for (int j = 0; j < tleaves.size(); j++)
        {
            allMyLeaves[tleaves[j]->GetLabel()] = tleaves[j];

            /*string lbl = tleaves[j]->GetLabel();
            if (lbl.find("sapien") == string::npos)
            {
                allHomo = false;
            }*/
        }
    }

    map<Node*, map<Node*, double> > mydistances;
    //copy the distances, we'll be playing with them a lot locally.
    //TODO : this here could be improved
    for (map<Node*, map<Node*, double> >::iterator dit = distances.begin(); dit != distances.end(); dit++)
    {
        Node* n = (*dit).first;

        if (allMyLeaves.find(n->GetLabel()) != allMyLeaves.end())
        {
            Node* myN1 = allMyLeaves[n->GetLabel()];
            map<Node*, double> dz = (*dit).second;
            map<Node*, double> dzcopy;

            for (map<Node*, double>::iterator dzit = dz.begin(); dzit != dz.end(); dzit++)
            {
                Node* n2 = (*dzit).first;

                if (allMyLeaves.find(n2->GetLabel()) != allMyLeaves.end())
                {
                    Node* myN2 = allMyLeaves[n2->GetLabel()];
                    dzcopy[myN2] = (*dzit).second;
                }
            }

            mydistances[myN1] = dzcopy;
        }
    }


    for (int i = 0; i < subtrees.size(); i++)
    {
        Node* t = subtrees[i];

        TreeIterator* it = t->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            //if we get here, we haven't computed the distance from n to the rest...
            //but recursively, we have previously computed the distances from n's children to the rest
            if (!n->IsLeaf() && mydistances.find(n) == mydistances.end())
            {
                Node* c1 = n->GetChild(0);
                Node* c2 = n->GetChild(1);

                double sumC1 = 0, sumC2 = 0;
                for (map<Node*, map<Node*, double> >::iterator sumIt = mydistances.begin(); sumIt != mydistances.end(); sumIt++)
                {
                    sumC1 += mydistances[c1][(*sumIt).first];
                    sumC2 += mydistances[c2][(*sumIt).first];
                }

                //cout<<"dist_nc1 = 1/2 * "<<mydistances[c1][c2]<<"+1/(2*("<<mydistances.size()<<" - 2))*("<<sumC1<<"-"<<sumC2<<")"<<endl;
                double dist_nc1 = 1.0/2.0 * mydistances[c1][c2] + 1.0/(2.0*(mydistances.size() - 2.0)) * (sumC1 - sumC2);
                double dist_nc2 = mydistances[c1][c2] - dist_nc1;

                //and here, we compute the dist from x to n for all current x.
                //but we can't add the distances into mydistances since we'll be iterating over it.
                map<Node*, double> distancesToN;
                for (map<Node*, map<Node*, double> >::iterator toNIt = mydistances.begin(); toNIt != mydistances.end(); toNIt++)
                {
                    Node* x = (*toNIt).first;
                    if (x != c1 && x != c2)
                    {
                        distancesToN[x] = 1.0/2.0 * (mydistances[x][c1] + mydistances[x][c2] - mydistances[c1][c2]);
                    }
                }

                mydistances[n] = map<Node*, double>();

                for (map<Node*, double>::iterator disttoNIt = distancesToN.begin(); disttoNIt != distancesToN.end(); disttoNIt++)
                {
                    Node* x = (*disttoNIt).first;
                    mydistances[x][n] = (*disttoNIt).second;
                    mydistances[n][x] = (*disttoNIt).second;
                }

                mydistances[c1][n] = dist_nc1;
                mydistances[n][c1] = dist_nc1;
                mydistances[c2][n] = dist_nc2;
                mydistances[n][c2] = dist_nc2;

                mydistances.erase(c1);
                mydistances.erase(c2);

            }
        }
        t->CloseIterator(it);
    }


    //debug
    /*if (allHomo)
    {
    for (map<Node*, map<Node*, double> >::iterator ttit = mydistances.begin(); ttit != mydistances.end(); ttit++)
    {
        Node* n1 = (*ttit).first;
        map<Node*, double> dd = (*ttit).second;

        for (map<Node*, double>::iterator ddit = dd.begin(); ddit != dd.end(); ddit++)
        {
            Node* n2 = (*ddit).first;
            double thed = (*ddit).second;

            cout<<n1->GetLabel()<<" "<<n2->GetLabel()<<" "<<thed<<endl;
        }
    }
    }*/

    allMyLeaves.clear();
    return mydistances;
}



void PolySolverDistance::PrintScoresTable(Node* speciesTree, int maxK, map<Node*, map<int, PSDCell*> > &cells)
{
    string str = "";
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    int cnts = 0;
    while (Node* s = it->next())
    {
        cnts++;
        set<Node*> leaves = s->GetLeafSet();
        string sLabel = "";
        for (set<Node*>::iterator lit = leaves.begin(); lit != leaves.end(); lit++)
        {

            sLabel += (*lit)->GetLabel().substr(0, 1) + " ";
        }
        string line = "";

        line += sLabel + " : ";

        if (cells.find(s) != cells.end())
        {
            map<int, PSDCell*> sCells = cells[s];

            for (int k = 1; k <= maxK; k++)
            {
                if (sCells.find(k) != sCells.end())
                {
                    PSDCell* cell = sCells[k];
                    //line += Util::ToString(k) + "/" + Util::ToString(cell->CVal) + "/" + Util::ToString(cell->MVal) + "   |   ";
                    line += Util::ToString(cell->MVal) + "   |   ";
                }
            }
        }

        if (str != "")
            str += "\n";

        str += line;

    }
    speciesTree->CloseIterator(it);

    cout<<"CNTS="<<cnts<<endl;
    cout<<str<<endl;
}



set<Node*> PolySolverDistance::ComputeLinkSpeciesLeaves(vector<Node*> polytomyChildren, unordered_map<Node*, Node*> &lcaMapping, Node* speciesTree)
{
    set<Node*> leafset;
    for (int p = 0; p < polytomyChildren.size(); p++)
    {
        Node* g = polytomyChildren[p];
        Node* s = lcaMapping[g];

        leafset.insert(s);

        //yes, in case we miss something, we add every sibling of an ancestor of s
        //they need to be there for the table computation
        //If they're useless, the next while loop removes them

        while (s != speciesTree)
        {
            Node* sibl = s->GetLeftSibling();
            if (!sibl)
                sibl = s->GetRightSibling();
            leafset.insert(sibl);
            s = s->GetParent();
        }

    }


    //remove leaves that are actually strict ancestors of other chosen leaves
    set<Node*>::iterator it1 = leafset.begin();
    while (it1 != leafset.end())
    {

        Node* s = (*it1);

        if (this->verbose > 1)
        {
            cout<<"Considering s="<<s->GetLabel()<<endl;
        }

        bool removed = false;
        for (set<Node*>::iterator it2 = leafset.begin(); it2 != leafset.end() && !removed; it2++)
        {
            Node* s2 = (*it2);

            if (s != s2 && s2->HasAncestor(s))
            {
                it1 = leafset.erase(it1);
                removed = true;
            }
        }

        if (!removed)
        {
            it1++;
        }


    }


    return leafset;

}
