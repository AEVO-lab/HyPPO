//////////////////////////////
// THIS FILE IS FOR ARCHIVE PURPOSE ONLY

#ifdef OH_HI_PLEASE_IGNORE_ME

#include "trees/node.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"

#include <iostream>
#include <string>
#include <map>
#include <time.h>
#include <algorithm>
#include "div/tinydir.h"


using namespace std;


string g_gene_label_separator = "__";
int g_gene_label_species_index = 1;


string GetSpeciesName(string geneLabel)
{
    vector<string> pz = Util::Split(geneLabel, g_gene_label_separator, false);
    return pz[g_gene_label_species_index];
}

double GetWeight(unordered_map<string, double> weights, string gene1, string gene2)
{
    string key = gene1 + ";;" + gene2;

    if (weights.find(key) == weights.end())
        return 0.0;

    return weights[key];
}


class GeneTree
{

public:

    GeneTree(Node* geneTreeRoot, Node* speciesTree, unordered_map<Node*, Node*> &lca_mapping)
    {
        this->geneTreeRoot = geneTreeRoot;
        this->speciesTree = speciesTree;
        this->lca_mapping = lca_mapping;
    }

    Node* GetLCAMapping(Node* geneTreeNode)
    {
        return lca_mapping[geneTreeNode];
    }

    //TODO: MUST GRAFT AT SPECIFIC NODE SO THAT LCAMAPPING REMAINS OK
    //SEE main() FOR USAGE
    Node* GraftDuplicationGeneTree(Node* node_to_graft_above, GeneTree &geneTree)
    {
        Node* newnode = NULL;
        if (!node_to_graft_above->IsRoot())
        {
            newnode = node_to_graft_above->GraftOnParentEdge(geneTree.geneTreeRoot);
        }
        else
        {
            newnode = new Node(false);
            newnode->AddSubTree(this->geneTreeRoot);
            newnode->AddSubTree(geneTree.geneTreeRoot);
            this->geneTreeRoot = newnode;
        }

        newnode->SetCustomField("event", "duplication");
        lca_mapping[newnode] = geneTree.GetLCAMapping(geneTree.geneTreeRoot);
        lca_mapping.insert(geneTree.lca_mapping.begin(), geneTree.lca_mapping.end());

        return newnode;
    }

    Node* speciesTree;
    Node* geneTreeRoot;
    unordered_map<Node*, Node*> lca_mapping;
};


class OrthologyCluster
{

public:
    OrthologyCluster(vector<string> &geneNames, Node* speciesTree)
    {
        this->geneNames = geneNames;
        this->speciesTree = speciesTree;
    }

    ~OrthologyCluster()
    {

    }

    int GetNbGenes()
    {
        return geneNames.size();
    }

    string GetGeneName(int i)
    {
        return geneNames[i];
    }


    Node* GetLCASpecies()
    {
        if (!lcaSpecies)    //make a cache
        {
            //TODO: handle errors
            vector<Node*> presentSpecies;
            for (int i = 0; i < geneNames.size(); i++)
            {

                string spname = GetSpeciesName(geneNames[i]);
                Node* snode = speciesTree->GetTreeInfo()->GetNodeByLabel(spname);
                presentSpecies.push_back( snode );
            }
            lcaSpecies = speciesTree->FindLCA(presentSpecies);
        }

        return lcaSpecies;
    }

    GeneTree BuildGeneTree()
    {
        //We just copy the species tree, and keep leaves of species that have at least one gene
        //Then we add genes accordingly (all genes from the same species are just clustered
        //together at the bottom - we are kind of assuming one gene per species here).



        //TODO: handle errors
        unordered_map<string, vector<string> > genesPerSpecies;
        for (int i = 0; i < geneNames.size(); i++)
        {
            string spname = GetSpeciesName(geneNames[i]);
            Node* snode = speciesTree->GetTreeInfo()->GetNodeByLabel(spname);

            if (genesPerSpecies.find( snode->GetLabel() ) == genesPerSpecies.end())
            {
                genesPerSpecies[snode->GetLabel()] = vector<string>();
            }
            //threee lines below because C++ is C++, and makes copies of vectors whenever it gets a chance
            vector<string> tmp = genesPerSpecies[ snode->GetLabel() ];
            tmp.push_back(geneNames[i]);
            genesPerSpecies[ snode->GetLabel() ] = tmp;

        }

        Node* splca = GetLCASpecies();

        Node* scopy = new Node(false);
        scopy->CopyFrom(splca);

        vector<Node*> markedForDeletion;
        TreeIterator* it = scopy->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            if (!n->IsLeaf())
            {
                n->SetCustomField("event", "speciation");
            }
            else
            {
                //Node* spnode = speciesTree->GetTreeInfo()->GetNodeByLabel(n->GetLabel());
                if (genesPerSpecies.find(n->GetLabel()) == genesPerSpecies.end())
                {
                    markedForDeletion.push_back(n);   //will delete later, not nice to delete while iterating
                }
                else
                {
                    vector<string> genes = genesPerSpecies[n->GetLabel()];
                    if (genes.size() == 1)
                    {
                        //geneSpeciesMapping[n] = spnode;
                        n->SetLabel(genes[0]);
                    }
                    else
                    {

                        n->SetLabel("");
                        n->SetCustomField("event", "duplication");
                        for (int i = 0; i < genes.size(); i++)
                        {
                            Node* ch = n->AddChild();
                            ch->SetLabel(genes[i]);
                            //geneSpeciesMapping[ch] = spnode;
                        }
                    }
                }
            }
        }
        scopy->CloseIterator(it);


        //cleanup useless leaves
        for (int i = 0; i < markedForDeletion.size(); i++)
        {
            Node* n = markedForDeletion[i];
            n->GetParent()->RemoveChild(n);
            delete n;
        }
        scopy->DeleteSingleChildDescendants();

        unordered_map<Node*, Node*> lca_mapping;

        lca_mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(scopy, speciesTree, g_gene_label_separator, g_gene_label_species_index);

        GeneTree g(scopy, speciesTree, lca_mapping);
        return g;
    }

private:

    vector<string> geneNames;
    Node* speciesTree;
    Node* lcaSpecies;

};




bool cluster_compare_ascending(OrthologyCluster* c1, OrthologyCluster* c2)
{
    Node* s1 = c1->GetLCASpecies();
    Node* s2 = c2->GetLCASpecies();
    return (s1->HasAncestor(s2));
}


bool cluster_compare_descending(OrthologyCluster* c1, OrthologyCluster* c2)
{
    Node* s1 = c1->GetLCASpecies();
    Node* s2 = c2->GetLCASpecies();
    return (s2->HasAncestor(s1));
}


// The best paralogies cluster is the one that maximizes the average weight between
// genes pairs that are both from species below mycluster->GetLCASpecies().
// Normally we'd compute the average between all possible gene pairs, but this
// adds too much noise if clusters are big.
OrthologyCluster* GetBestParalogiesCluster(OrthologyCluster* mycluster, vector<OrthologyCluster*> clusters, unordered_map<string, double> weights, Node* speciesTree)
{
    Node* lca_species = mycluster->GetLCASpecies();

    OrthologyCluster* bestCluster = NULL;
    double bestAverageWeight = 0;

    for (int c = 0; c < clusters.size(); c++)
    {
        OrthologyCluster* cluster = clusters[c];
        double totalWeight = 0.0;
        int nbPairs = 0;

        for (int i = 0; i < mycluster->GetNbGenes(); i++)
        {
            string gene1 = mycluster->GetGeneName(i);

            for (int j = 0; j < cluster->GetNbGenes(); j++)
            {
                string gene2 = cluster->GetGeneName(j);
                string sname2 = GetSpeciesName(gene2);
                Node* s = speciesTree->GetTreeInfo()->GetNodeByLabel(sname2);

                if (s->HasAncestor(lca_species))
                {
                    totalWeight += GetWeight(weights, gene1, gene2);
                    nbPairs++;
                }

            }
        }


        double avg = totalWeight/((double)nbPairs);

        if (avg > bestAverageWeight)
        {
            bestAverageWeight = avg;
            bestCluster = cluster;
        }
    }

    return bestCluster;

}


// Computes all orthologs we would get for oc if we grafted the oc subtee above
// node_we_graft_above.  Returns the sum of the weights of orthologous pairs
// containing one member of oc, and one member not in oc.
// Assumes the gene tree is labeled by events.
double GetOrthologyWeightForGrafting(OrthologyCluster* oc, unordered_map<string, double> weights, GeneTree geneTree, Node* node_we_graft_above)
{
    Node* prevnode = node_we_graft_above;
    Node* curnode = node_we_graft_above->GetParent();

    double totalWeight = 0;

    while (curnode)
    {
        if (curnode->GetCustomField("event") == "speciation")
        {
            for (int i = 0; i < curnode->GetNbChildren(); i++)
            {
                Node* child = curnode->GetChild(i);
                if (child != prevnode)
                {
                    //all descendants of child are thus orthologs
                    vector<Node*> leaves = child->GetLeafVector();

                    for (int c_oc = 0; c_oc < oc->GetNbGenes(); c_oc++)
                    {
                        string ocgene = oc->GetGeneName(c_oc);

                        for (int c_leaf = 0; c_leaf < leaves.size(); c_leaf++)
                        {
                            string leafgene = leaves[c_leaf]->GetLabel();

                            string key = ocgene + ";;" + leafgene;
                            totalWeight += weights[key];
                        }
                    }
                }
            }
        }

        prevnode = curnode;
        curnode = curnode->GetParent();
    }

    return totalWeight;
}




int main(int argc, char *argv[])
{

    int verbose = 0;

    map<string, string> args;

    //BUILD DICTIONARY OF ARGS
    string prevArg = "";
    for (int i = 0; i < argc; i++)
    {
        if (string(argv[i]) == "-v")
        {
            verbose = 1;
            prevArg = "";
        }
        else
        {
            if (prevArg != "" && prevArg[0] == '-')
            {
                args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
            }

            prevArg = string(argv[i]);
        }
    }




    string scontent = "";

    if (args.find("s") != args.end())
    {
        scontent = Util::GetFileContent(args["s"]);
    }

    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);


    unordered_map< string, double > orthologyWeights;   //format for string key is GENEID1;;GENEID2
    vector<OrthologyCluster*> clusters;

    std::sort(clusters.begin(), clusters.end(), cluster_compare_descending);

    GeneTree gtree = clusters[0]->BuildGeneTree();

    //NOTE: this will modify the gene trees that are members of the OrthologyCluster
    //instance, as we are dealing with pointer references.  But, we don't really care here.
    for (int c = 1; c < clusters.size(); c++)
    {
        OrthologyCluster* cluster = clusters[c];
        Node* sp_cluster = cluster->GetLCASpecies();

        //graft at best spot.  We only consider uv such that v <= cluster and u > cluster
        Node* best_node_to_graft_above = NULL;
        double bestTotalWeight = 0;

        Node* tree = gtree.geneTreeRoot;

        TreeIterator* it = tree->GetPostOrderIterator();
        while (Node* node = it->next())
        {
            if (!node->IsRoot())
            {
                Node* sp_node = gtree.GetLCAMapping(node);
                Node* sp_above = gtree.GetLCAMapping(node->GetParent());

                if (( sp_node->HasAncestor(sp_cluster) || sp_node == sp_cluster)
                        && sp_cluster->HasAncestor(sp_above) && sp_cluster != sp_above)
                {
                    double weight = GetOrthologyWeightForGrafting(cluster, orthologyWeights, gtree, node);
                    if (weight > bestTotalWeight || best_node_to_graft_above == NULL)
                    {
                        best_node_to_graft_above = node;
                        bestTotalWeight = weight;
                    }
                }
            }
        }
        tree->CloseIterator(it);


        //We've found the best grafting spot.  Now apply it.
        GeneTree cluster_gtree = cluster->BuildGeneTree();
        if (best_node_to_graft_above)
        {
            gtree.GraftDuplicationGeneTree(best_node_to_graft_above, cluster_gtree);
        }
        else
        {
            gtree.GraftDuplicationGeneTree(tree, cluster_gtree);
        }

    }






    //wrap it up
    for (int i = 0; i < clusters.size(); i++)
    {
        delete clusters[i];
    }

    delete gtree.geneTreeRoot;
    delete speciesTree;



    return 0;

}

#endif
