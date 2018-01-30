#include "cluster.h"



//string g_cluster_species_separator = "__";
//int g_cluster_species_index = 0;
string Cluster::speciesSeparator = "__";
int Cluster::speciesIndex = 0;

GeneTree::GeneTree()
{
    this->geneTreeRoot = NULL;
    this->speciesTree = NULL;
}

GeneTree::GeneTree(Node* geneTreeRoot, Node* speciesTree, unordered_map<Node*, Node*> &lca_mapping)
{
    this->geneTreeRoot = geneTreeRoot;
    this->speciesTree = speciesTree;
    this->lca_mapping = lca_mapping;
}

bool GeneTree::IsDupUnderLCAMapping(Node* node)
{
    for (int i = 0; i < node->GetNbChildren(); i++)
    {
        if (GetLCAMapping(node) == GetLCAMapping(node->GetChild(i)))
            return true;
    }
    return false;
}

Node* GeneTree::GetLCAMapping(Node* geneTreeNode)
{
    return lca_mapping[geneTreeNode];
}





Cluster::Cluster(vector<string> &geneNames, Node* speciesTree)
{
    Init(geneNames, speciesTree);
}

Cluster::~Cluster()
{
    //TODO: maybe delete gene tree stuff
}


void Cluster::Init(vector<string> &geneNames, Node* speciesTree)
{
    this->geneNames = geneNames;
    this->speciesTree = speciesTree;
    this->lcaSpecies = NULL;


    for (int g = 0; g < this->GetNbGenes(); g++)
    {
        string sp = GetSpeciesName(this->GetGeneName(g));
        this->geneToSpecies[this->GetGeneName(g)] = sp;

        if (this->speciesSet.find(sp) == this->speciesSet.end())
        {
            this->speciesSet.insert(sp);
        }
    }
}

string Cluster::GetGeneSpecies(string g)
{
    return this->geneToSpecies[g];
}

set<string> Cluster::GetSpeciesSet()
{
    return this->speciesSet;
}

int Cluster::GetNbGenes()
{
    return geneNames.size();
}

string Cluster::GetSpeciesName(string geneLabel)
{
    vector<string> pz = Util::Split(geneLabel, Cluster::speciesSeparator, false);
    return pz[Cluster::speciesIndex];
}

string Cluster::GetGeneName(int i)
{
    return geneNames[i];
}


Node* Cluster::GetLCASpecies()
{

    if (!lcaSpecies)    //make a cache
    {
        //TODO: handle errors
        vector<Node*> presentSpecies;
        for (int i = 0; i < geneNames.size(); i++)
        {

            string spname = GetSpeciesName(geneNames[i]);
            Node* snode = speciesTree->GetLeafByLabel(spname);

            if (!snode)
            {
                cout<<"COULD NOT FIND SPECIES "<<spname<<"...this is bad"<<endl;
                throw "COULD NOT FIND SPECIES " + spname;
            }

            presentSpecies.push_back( snode );
        }

        lcaSpecies = speciesTree->FindLCA(presentSpecies);

    }

    return lcaSpecies;
}

GeneTree Cluster::GetGeneTree()
{
    if (!geneTree.geneTreeRoot)
    {
        BuildGeneTree();
    }

    return geneTree;
}

void Cluster::BuildGeneTree()
{




    //We just copy the species tree, and keep leaves of species that have at least one gene
    //Then we add genes accordingly (all genes from the same species are just clustered
    //together at the bottom - we are kind of assuming one gene per species here).




    //TODO: handle errors
    unordered_map<string, vector<string> > genesPerSpecies;
    for (int i = 0; i < geneNames.size(); i++)
    {
        string spname = GetSpeciesName(geneNames[i]);

        if (genesPerSpecies.find( spname ) == genesPerSpecies.end())
        {
            genesPerSpecies[spname] = vector<string>();
        }
        //three lines below because C++ is C++, and makes copies of vectors whenever it gets a chance
        vector<string> tmp = genesPerSpecies[ spname ];
        tmp.push_back(geneNames[i]);
        genesPerSpecies[ spname ] = tmp;

    }

    Node* splca = GetLCASpecies();

    //NOTE: scopy will become the gene tree root after relabeling
    Node* scopy = new Node(false);
    scopy->CopyFrom(splca);

    set<Node*> leavesToKeep;
    TreeIterator* it = scopy->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            n->SetCustomField("event", "speciation");
        }
        else
        {
            //does this species have 0 genes?
            if (genesPerSpecies.find(n->GetLabel()) == genesPerSpecies.end())
            {
                //markedForDeletion.push_back(n);   //will delete later, not nice to delete while iterating	//ML: idk why it's commented out
            }
            else
            {
                vector<string> genes = genesPerSpecies[n->GetLabel()];
                if (genes.size() == 1)
                {
                    //geneSpeciesMapping[n] = spnode;
                    n->SetLabel(genes[0]);
                    leavesToKeep.insert(n);
                }
                else
                {

                    n->SetLabel("Dup");
                    n->SetCustomField("event", "duplication");


                    for (int i = 0; i < genes.size(); i++)
                    {
                        Node* ch = n->AddChild();
                        ch->SetLabel(genes[i]);
                        leavesToKeep.insert(ch);
                        //geneSpeciesMapping[ch] = spnode;
                    }
                }
            }
        }
    }
    scopy->CloseIterator(it);

    Node::RestrictToLeafset(scopy, leavesToKeep);

    //cleanup useless leaves
    /*for (int i = 0; i < markedForDeletion.size(); i++)
    {
        Node* n = markedForDeletion[i];
        n->GetParent()->RemoveChild(n);
        delete n;
    }
    scopy->DeleteSingleChildDescendants();*/

    unordered_map<Node*, Node*> lca_mapping;

    lca_mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(scopy, speciesTree, Cluster::speciesSeparator, Cluster::speciesIndex);

    //NOTE: recall that scopy is the root of the gene tree now
    GeneTree g(scopy, speciesTree, lca_mapping);

    this->geneTree = g;
}

void Cluster::SetGeneTree(Node* gtree, Node* speciesTree)
{
    unordered_map<Node*, Node*> lca_mapping;
    lca_mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(gtree, speciesTree, Cluster::speciesSeparator, Cluster::speciesIndex);
    GeneTree g(gtree, speciesTree, lca_mapping);
    this->geneTree = g;
}


void Cluster::UpdateLCAMapping()
{
    unordered_map<Node*, Node*> lca_mapping;
    lca_mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(this->geneTree.geneTreeRoot, this->geneTree.speciesTree,
                                                                 Cluster::speciesSeparator, Cluster::speciesIndex);
    this->geneTree.lca_mapping = lca_mapping;

    lcaSpecies = NULL;
}


string Cluster::ToString()
{
    string str = "";
    for (int i = 0; i < this->GetNbGenes(); i++)
    {
        str += this->GetGeneName(i) + " ";
    }
    return str;
}


set<string> Cluster::GetSpeciesIntersection(set<string> species)
{
    return Util::GetSetIntersection(this->speciesSet, species);
}



Cluster *Cluster::GetCompressedCluster(map<string, set<string> > &blocks, string new_gene_name)
{
    vector<string> gene_names;

    for (map<string, set<string> >::iterator it = blocks.begin(); it != blocks.end(); it++)
    {
        string name = (*it).first;
        set<string> block = (*it).second;

        if (this->GetSpeciesIntersection(block).size() > 0)
        {
			string newname = name + Cluster::speciesSeparator + new_gene_name;
			if (Cluster::speciesIndex == 1)
			{
				newname = new_gene_name + Cluster::speciesSeparator + name;
			}
            gene_names.push_back(newname);
        }
    }

    return new Cluster(gene_names, NULL);
}





set<string> Cluster::GetAllSpecies(vector<Cluster*> clusters)
{
    set<string> species;
    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];
        for (int g = 0; g < cl->GetNbGenes(); g++)
        {
            string sp = cl->GetGeneSpecies(cl->GetGeneName(g));

            if (species.find(sp) == species.end())
                species.insert(sp);
        }
    }
    return species;
}


bool Cluster::HasSpecies(string sp)
{
    for (set<string>::iterator it = this->speciesSet.begin(); it != this->speciesSet.end(); it++)
    {
        if (sp == (*it))
            return true;
    }
    return false;
}


map< string, set<string> > Cluster::GetBlocks(vector<Cluster*> &clusters)
{
    map< string, set<string> > blocks;

    set<string> l_speciesSet = Cluster::GetAllSpecies(clusters);
	
	
    while (l_speciesSet.size() > 0)
    {
        string sp = (*l_speciesSet.begin());

        set<string> curinter;
        set<string> forbidden;

        for (int c = 0; c < clusters.size(); c++)
        {
            Cluster* cl = clusters[c];

            if (cl->HasSpecies(sp))
            {
                if (curinter.size() == 0)
                {
                    curinter = cl->GetSpeciesSet();
                }
                else
                {
                    curinter = cl->GetSpeciesIntersection(curinter);
                }
            }
            else
            {
                set<string> sx = cl->GetSpeciesSet();
                for (set<string>::iterator it = sx.begin(); it != sx.end(); it++)
                {
                    if (forbidden.find(*it) == forbidden.end())
                        forbidden.insert(*it);
                }
            }

        }

        for (set<string>::iterator it = forbidden.begin(); it != forbidden.end(); it++)
        {
            curinter.erase(*it);
        }


        string name = "";
        name = NewickLex::GetCaterpillarNewick(curinter);
        /*for (set<string>::iterator it = curinter.begin(); it != curinter.end(); it++)
        {
            if (name != "") name += ";;";
            name += (*it);
        }
        name = "B" + name;*/

        blocks[name] = curinter;

        for (set<string>::iterator it = curinter.begin(); it != curinter.end(); it++)
        {
            l_speciesSet.erase(*it);
        }
    }

    return blocks;

}



