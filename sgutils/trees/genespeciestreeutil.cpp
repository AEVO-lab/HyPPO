#include "genespeciestreeutil.h"

//This is used for pruning species trees.  The Node::Restrict function takes a fct pointer as argument.
bool SpeciesTreePruning__RestrictFunction(Node* n, void *arg)
{
    set<Node*>* leavesToKeep = (set<Node*>*)arg;
    if (n->IsLeaf())
    {
        return (leavesToKeep->find(n) != leavesToKeep->end());
    }
    else
    {
        return (n->IsRoot() || n->GetNbChildren() > 1);
    }
}






GeneSpeciesTreeUtil::GeneSpeciesTreeUtil()
{
}



unordered_map<Node*, Node*> GeneSpeciesTreeUtil::GetLCAMapping(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> &geneLeavesSpeciesMapping)
{
    unordered_map<Node*, Node*> lcaMapping;
    TreeIterator* it = geneTree->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        if (g->IsLeaf())
        {
            lcaMapping[g] = geneLeavesSpeciesMapping[g];
        }
        else
        {
            vector<Node*> v;
            for (int i = 0; i < g->GetNbChildren(); i++)
            {
                v.push_back(lcaMapping[g->GetChild(i)]);
            }
            lcaMapping[g] = speciesTree->FindLCA(v);
        }
    }
    geneTree->CloseIterator(it);

    return lcaMapping;
}



unordered_map<Node*, Node*> GeneSpeciesTreeUtil::GetLCAMapping(Node *geneTree, Node *speciesTree, string geneLabelSeparator, int speciesIndex)
{

    unordered_map<Node*, Node*> m = this->GetGeneSpeciesMappingByLabel(geneTree, speciesTree, geneLabelSeparator, speciesIndex);

    return this->GetLCAMapping(geneTree, speciesTree, m);
}


unordered_set<Node*> GeneSpeciesTreeUtil::GetGeneTreeSpecies(Node *geneTree, unordered_map<Node*, Node*> &lcaMapping)
{
    unordered_set<Node*> species;
    TreeIterator* it = geneTree->GetPostOrderIterator(true);
    while (Node* g = it->next())
    {

        species.insert(lcaMapping[g]);
    }
    geneTree->CloseIterator(it);

    return species;
}


vector<Node*> GeneSpeciesTreeUtil::GetGenesSpecies(vector<Node*> genes, unordered_map<Node*, Node*> &lcaMapping)
{
    vector<Node*> species;
    for (int i = 0; i < genes.size(); i++)
    {
        species.push_back(lcaMapping[genes[i]]);
    }

    return species;
}

//-gl 1 -s "/u/lafonman/Projects/PolytomySolverDistance_1.2.4/PolytomySolverDistance_1.2.4/example_files/nonstar/Compara.73.species_tree" -g "/u/lafonman/Projects/PolytomySolverDistance_1.2.4/PolytomySolverDistance_1.2.4/example_files/nonstar/famille_1.start_tree" -d "/u/lafonman/Projects/PolytomySolverDistance_1.2.4/PolytomySolverDistance_1.2.4/example_files/nonstar/famille_1.dist" -r none -n -v

unordered_map<Node*, Node*> GeneSpeciesTreeUtil::GetGeneSpeciesMappingByLabel(Node* geneTree, Node* speciesTree, string separator, int speciesIndex)
{
    unordered_map<Node*, Node*> mapping;

    TreeIterator* it = geneTree->GetPostOrderIterator(true);
    while (Node* g = it->next())
    {
        string lbl = g->GetLabel();
        vector<string> sz = Util::Split(g->GetLabel(), separator);

        if (sz.size() == 0 || (speciesIndex > 0 && sz.size() < speciesIndex - 1))   //there was a bug when speciesIndex = 0
        {
            cout<<"Gene label "<<lbl<<" malformed"<<endl<<flush;
            throw "Gene label " + lbl + " malformed.";
        }
        Node* s = speciesTree->GetLeafByLabel(sz[speciesIndex]);

        if (s)
        {
            mapping[g] = s;
        }
        else
        {
            string msg = "Could not find species for gene " + lbl +
                    "  S=" + NewickLex::ToNewickString(speciesTree);
            cout<<msg<<endl;
            throw msg;
        }
    }
    geneTree->CloseIterator(it);

    return mapping;
}



Node* GeneSpeciesTreeUtil::CopyTreeWithNodeMapping(Node* tree, unordered_map<Node*, Node*> &yourMapping, unordered_map<Node*,Node*> &mappingToFill)
{
    Node* newTree = new Node(false);

    newTree->CopyFrom(tree);

    TreeIterator* it = tree->GetPostOrderIterator();
    TreeIterator* itNew = newTree->GetPostOrderIterator();

    while (Node* n = it->next())
    {
        Node* n2 = itNew->next();

        //n and n2 advance together, so they should point to corresponding nodes
        if (yourMapping.find(n) != yourMapping.end())
            mappingToFill[n2] = yourMapping[n];

    }
    tree->CloseIterator(it);
    newTree->CloseIterator(itNew);

    //debugging
    //GeneSpeciesTreeUtil::Instance()->PrintMapping(newTree, mappingToFill);

    return newTree;

}


bool GeneSpeciesTreeUtil::HaveCommonSpecies(Node* tree1, Node* tree2, unordered_map<Node*, Node*> &mapping)
{
    unordered_set<Node*> left = GetGeneTreeSpecies(tree1, mapping);
    unordered_set<Node*> right = GetGeneTreeSpecies(tree2, mapping);


    for (unordered_set<Node*>::iterator itLeft = left.begin(); itLeft != left.end(); itLeft++)
    {
        for (unordered_set<Node*>::iterator itRight = right.begin(); itRight != right.end(); itRight++)
        {

            if ((*itLeft) == (*itRight))
            {
                return true;
            }
        }
    }

    return false;
}


Node* GeneSpeciesTreeUtil::GetSingleNodeLCAMapping(Node* n, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping)
{
    vector<Node*> v;
    for (int i = 0; i < n->GetNbChildren(); i++)
    {
        v.push_back(lcaMapping[n->GetChild(i)]);
    }

    return speciesTree->FindLCA(v);
}


void GeneSpeciesTreeUtil::PrintMapping(Node* tree, unordered_map<Node*, Node*> &mapping)
{
    TreeIterator* it = tree->GetPostOrderIterator();

    while (Node* n = it->next())
    {
        if (mapping[n])
            cout<<n->GetLabel()<<" -> "<<mapping[n]->GetLabel()<<endl;
        else
            cout<<n->GetLabel()<<" -> "<<"NULL"<<endl;
    }
    tree->CloseIterator(it);

}




vector<Node*> GeneSpeciesTreeUtil::GetNADNodes(Node* g, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping)
{
    vector<Node*> res;
    //TODO : this is dirty and sloooow
    TreeIterator* it = g->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            //first check if it's a duplication, lca mapping classic rule
            if (lcaMapping[n->GetChild(0)] == lcaMapping[n] || lcaMapping[n->GetChild(1)] == lcaMapping[n])
            {
                if (!GeneSpeciesTreeUtil::Instance()->HaveCommonSpecies(n->GetChild(0), n->GetChild(1), lcaMapping))
                {
                    res.push_back(n);
                }
            }
        }
    }
    g->CloseIterator(it);

    return res;
}



int GeneSpeciesTreeUtil::GetDLScore(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> lcaMapping)
{
    int nbDups = 0;
    int nbLosses = 0;

    TreeIterator* it = geneTree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            Node* ns = lcaMapping[n];
            //n->SetLabel(lcaMapping[n]->GetLabel());
            bool isDup = false;
            for (int i = 0; i < n->GetNbChildren(); i++)
            {
                if (lcaMapping[n] == lcaMapping[n->GetChild(i)])
                {
                    isDup = true;
                }
            }

            if (isDup)
                nbDups++;

            //count everything that's missing between node and child edge

            for (int i = 0; i < n->GetNbChildren(); i++)
            {
                Node* c = n->GetChild(i);

                Node* cs = lcaMapping[c];


                int pathCount = 0;
                while (cs != ns)
                {
                    cs = cs->GetParent();
                    pathCount++;
                }
                nbLosses += (isDup ? pathCount : pathCount - 1);
            }
        }


    }
    geneTree->CloseIterator(it);

    LASTNBDUPS = nbDups;
    LASTNBLOSSES = nbLosses;

    return nbDups + nbLosses;
}


void GeneSpeciesTreeUtil::RelabelGenes(Node* geneTree, string search, string replace)
{
    TreeIterator* it = geneTree->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        /*vector<string> s = Util::Split(g->GetLabel(), separator);
        if (s.size() > geneNameIndex)
        {
            g->SetLabel(s[geneNameIndex]);
        }*/
        g->SetLabel(Util::ReplaceAll(g->GetLabel(), search, replace));
    }
    geneTree->CloseIterator(it);
}



void GeneSpeciesTreeUtil::RelabelGenesByIndex(Node* geneTree, string separator, int indexToKeep)
{
    TreeIterator* it = geneTree->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        vector<string> s = Util::Split(g->GetLabel(), separator);
        if (s.size() > indexToKeep)
        {
            g->SetLabel(s[indexToKeep]);
        }

    }
    geneTree->CloseIterator(it);
}


bool GeneSpeciesTreeUtil::IsNodeDup(Node* geneTreeNode, unordered_map<Node*, Node*> &lca_mapping)
{
    if (geneTreeNode->GetNbChildren() <= 1)
        return false;
    else if (geneTreeNode->GetNbChildren() == 2)
    {
        return (lca_mapping[geneTreeNode] == lca_mapping[ geneTreeNode->GetChild(0) ]
                ||
                lca_mapping[geneTreeNode] == lca_mapping[ geneTreeNode->GetChild(1) ]);
    }
    else
    {
        //check that all species are separated
        for (int i = 0; i < geneTreeNode->GetNbChildren(); i++)
        {
            for (int j  = i + 1; j < geneTreeNode->GetNbChildren(); j++)
            {
                Node* sp1 = lca_mapping[geneTreeNode->GetChild(i)];
                Node* sp2 = lca_mapping[geneTreeNode->GetChild(j)];

                Node* splca = sp1->FindLCAWith(sp2);

                if (splca == sp1 || splca == sp2)
                    return false;

            }
        }

        return true;
    }
}



//NOTE: requires tree info
int GeneSpeciesTreeUtil::GetNbLossesOnBranch(Node* speciesDown, Node* speciesUp, bool isDupTop)
{
    int depthDiff = speciesDown->GetDepth() - speciesUp->GetDepth();

    if (isDupTop)
        return depthDiff;
    else
    {
        if (depthDiff <= 0)
            throw "Speciation is inconsistent";
        return depthDiff - 1;
    }
}



void GeneSpeciesTreeUtil::PruneSpeciesTreeFromLCAMapping(Node* speciesTree, Node* geneTree, unordered_map<Node*, Node*> lca_mapping)
{
    set<Node*> speciesToKeep;

    TreeIterator* it = geneTree->GetPostOrderIterator(true);
    while (Node* leaf = it->next())
    {
        speciesToKeep.insert(lca_mapping[leaf]);
    }
    geneTree->CloseIterator(it);

    PruneSpeciesTree(speciesTree, speciesToKeep);
}


void GeneSpeciesTreeUtil::PruneSpeciesTree(Node* speciesTree, set<Node*> speciesToKeep)
{
    speciesTree->Restrict(&SpeciesTreePruning__RestrictFunction, (void*)&speciesToKeep);
}





vector<Node*> GeneSpeciesTreeUtil::GetGeneTreeHighestSpeciations(Node* geneTree, Node* speciesTree,
                                             unordered_map<Node*, Node*> lca_mapping)
{
    vector<Node*> specs;

    //root of gene tree is a spec => return it
    if (!GeneSpeciesTreeUtil::Instance()->IsNodeDup(geneTree, lca_mapping))
    {
        specs.push_back(geneTree);
    }
    //root is a dup => check for specs below
    else
    {
        for (int i = 0; i < geneTree->GetNbChildren(); i++)
        {
            vector<Node*> ch_specs = GetGeneTreeHighestSpeciations(geneTree->GetChild(i), speciesTree, lca_mapping);
            specs.insert(specs.end(), ch_specs.begin(), ch_specs.end());
        }
    }

    return specs;
}

void GeneSpeciesTreeUtil::LabelInternalNodesUniquely(Node* tree)
{
    vector<Node*> v;
    v.push_back(tree);
    LabelInternalNodesUniquely(v);
}

void GeneSpeciesTreeUtil::LabelInternalNodesUniquely(vector<Node*> trees)
{
    int cpt = 0;

    for (int t = 0; t < trees.size(); t++)
    {
        TreeIterator* it = trees[t]->GetPostOrderIterator();

        while (Node* n = it->next())
        {
            if (!n->IsLeaf())
            {
                string lbl = n->GetLabel();
                if (lbl != "")
                    lbl += "-";
                lbl += Util::ToString(cpt);

                n->SetLabel(lbl);

                cpt++;
            }
        }

        trees[t]->CloseIterator(it);
    }
}


void GeneSpeciesTreeUtil::LabelInternalNodesWithLCAMapping(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> lca_mapping)
{
    TreeIterator* it = geneTree->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        if (!g->IsLeaf())
            g->SetLabel(lca_mapping[g]->GetLabel());
    }

    geneTree->CloseIterator(it);
}


string GeneSpeciesTreeUtil::GetPrunedSpeciesTreeNewick(string gcontent, string scontent)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);



    speciesTree->DeleteTreeInfo();
    GeneSpeciesTreeUtil::Instance()->PruneSpeciesTreeFromLCAMapping(speciesTree, geneTree, lcamap);

    //special case: root has only one child
    if (speciesTree->GetNbChildren() == 1)
    {
        Node* prevroot = speciesTree;
        speciesTree = prevroot->GetChild(0);
        prevroot->RemoveChild(speciesTree);
        delete prevroot;
    }

    //this heavy gymnastic sucks, but there's no quick alternative
    string snewick = NewickLex::ToNewickString(speciesTree);
    delete speciesTree;
    delete geneTree;

    return snewick;
}
