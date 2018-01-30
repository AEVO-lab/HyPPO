#include "paralogycorrector.h"


/**
  This is privately used to delete degree two internal nodes or useless leaves
  **/
bool ParalogyCorrector__RestrictFunction(Node* n, void *arg)
{
    if (n->IsLeaf())
    {
        set<Node*>* leaves = (set<Node*>*)arg;
        return (leaves->find(n) != leaves->end());
    }
    else
    {
        return (n->IsRoot() || n->GetNbChildren() > 1);
    }
}


ParalogyCorrector::ParalogyCorrector()
{
}


void ParalogyCorrector::EnsureOrthologsSymnmetry(vector<pair<string, string> > &orthologs)
{
    set<string> prev;   //prev for previously encountered

    for (int i = 0; i < orthologs.size(); i++)
    {
        string rep = orthologs[i].first + ":;:" + orthologs[i].second;
        string comp = orthologs[i].second + ":;:" + orthologs[i].first;

        if (prev.find(comp) == prev.end())
        {
            prev.insert(rep);
        }
        else
        {
            prev.erase(comp);
        }

    }

    //Everything that's left in prev has no symmetric counterpart - so we add it
    for (set<string>::iterator it = prev.begin(); it != prev.end(); it++)
    {
        vector<string> sz = Util::Split(*it, ":;:");
        orthologs.push_back(make_pair(sz[1], sz[0]));
    }
}


Node* ParalogyCorrector::CorrectGeneTree(Node* geneTree, Node* speciesTree, map<Node*, Node*> &geneSpeciesMapping, vector<pair<string, string> > &orthologs)
{
    map<Node*, Node*> lcaMapping;
    set<Node*> forb;
    map<string, Node*> genesByName;
    return CorrectGeneTree(geneTree, speciesTree, geneSpeciesMapping, orthologs, lcaMapping, forb, genesByName);
}





Node* ParalogyCorrector::CorrectGeneTree(Node* geneTree, Node* speciesTree, map<Node*, Node*> &geneSpeciesMapping, vector<pair<string, string> > &orthologs, map<Node*, Node*> &lcaMapping, set<Node*> &forbiddenClades, map<string, Node*> &genesByName)
{

    //base case : leaves need no correction
    if (geneTree->IsLeaf())
    {
        Node* n = new Node(false);
        n->CopyFrom(geneTree);
        return n;
    }

    //This only needs to be called on first top call
    //Here we do the lca mapping and identify unpreservable subtrees (called forbiddenClades)
    if (geneTree->IsRoot())
    {
        EnsureOrthologsSymnmetry(orthologs);


        //first pass to get LCA mapping
        //Suboptimal : takes something like O(n^2) as the lca is not constant
        TreeIterator* it = geneTree->GetPostOrderIterator();
        while (Node* g = it->next())
        {
            if (g->IsLeaf())
            {
                lcaMapping[g] = geneSpeciesMapping[g];
                genesByName[g->GetLabel()] = g;
            }
            else
            {
                Node* left = lcaMapping[g->GetChild(0)];
                lcaMapping[g] = left->FindLCAWith(lcaMapping[g->GetChild(1)]);
            }
        }
        geneTree->CloseIterator(it);


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //STEP 1 : identify forbidden subtrees

        //For each false paralogs, find h_{a,b}, the highest node on path(a, LCA(a,b)) that has an lcaMapping
        //on path(m(a), LCA(m(a), m(b)))
        //In other words, h_{a,b} and its descendents are on the side of LCA(m(a), m(b))
        //To do this, we make gptr go up gradually on path(a, LCA(a, b)) in G,
        //and stop when m(gptr) crosses LCA(m(a), m(b))
        for (int i = 0; i < orthologs.size(); i++)
        {
            Node* g1 = genesByName[orthologs[i].first];
            Node* g2 = genesByName[orthologs[i].second];
            Node* glca = NULL;

            if (g1 && g2)
                glca = g1->FindLCAWith(g2);

            //if lca is already a speciation, don't bother with this constraint
            if (glca && (lcaMapping[glca] == lcaMapping[glca->GetChild(0)] || lcaMapping[glca] == lcaMapping[glca->GetChild(1)]))
            {
                Node* gptr = g1->GetParent();

                Node* s1 = geneSpeciesMapping[g1];
                Node* s2 = geneSpeciesMapping[g2];
                Node* slca = s1->FindLCAWith(s2);
                Node* sptr = s1;

                bool foundHighest = false;

                while (gptr != glca)
                {
                    if (!foundHighest)
                        {
                        while (sptr != lcaMapping[gptr])
                        {
                            sptr = sptr->GetParent();

                            if (sptr == slca)
                                foundHighest = true;
                        }
                    }

                    if (foundHighest)
                    {
                        forbiddenClades.insert(gptr);
                    }

                    gptr = gptr->GetParent();
                }
            }
        }

        //special case : no forbidden guys => everything is allright, return NULL
        if (forbiddenClades.size() == 0)
        {
            return NULL;
        }

    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //STEP 2 : identify highest preservable nodes

    set<Node*> highestPreservables;
    set<Node*> highestCorrectedSubtrees;    //result from correcting highestPreservable subtrees
    map<Node*, Node*> highestCorrectedMapping;  //maps highest preservables to their corrected version


    //avoid the root as a highest preservable
    for (int i = 0; i < geneTree->GetNbChildren(); i++)
    {
        FindHighestPreservableNodes(geneTree->GetChild(i), highestPreservables, forbiddenClades);
    }


    //Correct each highest preservable subtree
    for (set<Node*>::iterator itpres = highestPreservables.begin(); itpres != highestPreservables.end(); itpres++)
    {
        //Note : subpres is a corrected copy of the highestPreservable
        Node* subpres = CorrectGeneTree((*itpres), speciesTree, geneSpeciesMapping, orthologs, lcaMapping, forbiddenClades, genesByName);
        highestCorrectedMapping[(*itpres)] = subpres;
        lcaMapping[subpres] = lcaMapping[(*itpres)];

        highestCorrectedSubtrees.insert(subpres);
    }




    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STEP 3 : PURE MAX ORTHOLOGY : connects every corrected highest preservable maximizing orthology
    // s copies the subtrees, so we can then clean out the corrected trees from memory
    Node* s = GetMaxOrthologyTree(highestCorrectedSubtrees, lcaMapping, speciesTree);
    lcaMapping[s] = lcaMapping[geneTree];

    set<Node*>::iterator itHCS = highestCorrectedSubtrees.begin();
    while (itHCS != highestCorrectedSubtrees.end())
    {
        delete (*itHCS);
        itHCS++;
    }

    highestPreservables.clear();
    highestCorrectedSubtrees.clear();
    highestCorrectedMapping.clear();

    return s;

}





//Used by an old version of the algorithm - kept for historical purposes
/*void ParalogyCorrector::DoDFS(set<Node*> &unvisited, set<Node*> &curCC, map<Node*, set<Node*> > edges, Node* lastVisited)
{
    set<Node*> neigh = edges[lastVisited];

    set<Node*>::iterator it = neigh.begin();
    while (it != neigh.end())
    {
        Node* n = (*it);

        if (unvisited.find(n) != unvisited.end())
        {
            curCC.insert(n);
            unvisited.erase(n);
            DoDFS(unvisited, curCC, edges, n);
        }

        it++;
    }
}*/





Node* ParalogyCorrector::GetMaxOrthologyTree(set<Node*> &subtrees, map<Node*, Node*> &lcaMapping, Node* speciesTree)
{
    //As in the paper, we graft our subtrees where they belong on the species tree.  Then we clean out unnecessary stuff.

    map<Node*, Node*> speciesToSkeletonMapping;

    //TODO : we should only copy speciesTree strating from LCA(subtrees) - we only add more nodes to delete later on
    Node* skeleton = new Node(false);
    Node* retTree = skeleton;   //retTree always points to the root of the tree to return, which might be other than skeleton
    skeleton->CopyFrom(speciesTree);

    set<Node*> originalLeaves = skeleton->GetLeafSet();

    //TODO : is this a stupid way of mapping nodes of skeleton to nodes of S ?  It does work at least
    TreeIterator* itSkel = skeleton->GetPostOrderIterator();
    TreeIterator* itSp = speciesTree->GetPostOrderIterator();

    while (Node* skel = itSkel->next())
    {
        Node* sp = itSp->next();

        lcaMapping[skel] = sp;
        speciesToSkeletonMapping[sp] = skel;
    }

    skeleton->CloseIterator(itSkel);
    speciesTree->CloseIterator(itSp);



    set<Node*>::iterator subIt = subtrees.begin();
    while (subIt != subtrees.end())
    {
        Node* toInsert = new Node(false);
        toInsert->CopyFrom((*subIt));
        lcaMapping[toInsert] = lcaMapping[(*subIt)];

        //paper proves we need to put subtree right above r, as defined below
        Node* toInsertMapping = lcaMapping[toInsert];
        Node* r = speciesToSkeletonMapping[lcaMapping[toInsert]];

        if (r->IsRoot())
        {
            retTree = new Node(false);
            retTree->AddSubTree(skeleton);
            retTree->AddSubTree(toInsert);
        }
        else
        {
            Node* x = new Node(false);
            Node* p = r->GetParent();
            p->RemoveChild(r);

            x->AddSubTree(r);
            x->AddSubTree(toInsert);

            p->AddSubTree(x);

        }

        subIt++;
    }



    //we restrict to the complement of the leafset, otherwise it causes problems - believe me
    //It's easier to restrict than to remove, because briefly, some internal nodes become leaves
    set<Node*> newLeaves = retTree->GetLeafSet();
    for (set<Node*>::iterator oit = originalLeaves.begin(); oit != originalLeaves.end(); oit++)
    {
        newLeaves.erase(*oit);
    }

    retTree->Restrict(&ParalogyCorrector__RestrictFunction, (void*)&newLeaves);

    //special case when the root wasn't deleted by has only one child.  REMOVE IT !
    if (retTree->GetNbChildren() == 1)
    {
        Node* good = new Node(false);
        good->CopyFrom(retTree->GetChild(0));
        delete retTree;
        retTree = good;
    }



    return retTree;
}




void ParalogyCorrector::FindHighestPreservableNodes(Node* current, set<Node*> &highest, set<Node*> &forbiddenNodes)
{
    if (forbiddenNodes.find(current) != forbiddenNodes.end())
    {
        for (int i = 0; i < current->GetNbChildren(); i++)
        {
            FindHighestPreservableNodes(current->GetChild(i), highest, forbiddenNodes);
        }
    }
    else
    {
        highest.insert(current);
    }
}



