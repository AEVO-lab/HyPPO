
#include "trees/node.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"

#include <iostream>
#include <string>
#include <map>
#include <time.h>
#include <algorithm>
#include "div/tinydir.h"

#include "supergenetreemaker.h"
#include "cluster.h"
#include "speciestreemaker.h"
#include "clusterfinder.h"

using namespace std;

string g_species_separator = "__";
int g_species_index = 0;

//TODO: this was copied from sgtutils
//Needs duplication node labels to end with __Dup
vector<string> GetOrthologsParalogs(Node* gtree)
{
    vector<string> relations;

    TreeIterator* it = gtree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            bool isdup = false;
            string relstr = "Orthologs";
			
            if (Util::EndsWith(n->GetLabel(), "Dup"))
            {
                isdup = true;
                relstr = "Paralogs";
				
            }




            for (int c = 0; c < n->GetNbChildren(); c++)
            {
                vector<Node*> c_leaves = n->GetChild(c)->GetLeafVector();
                for (int d = c + 1; d < n->GetNbChildren(); d++)
                {
                    vector<Node*> d_leaves = n->GetChild(d)->GetLeafVector();

                    for (int ic = 0; ic < c_leaves.size(); ic++)
                    {
                        for (int id = 0; id < d_leaves.size(); id++)
                        {
                            string gc = c_leaves[ic]->GetLabel();
                            string gd = d_leaves[id]->GetLabel();

                            string pz = gc + "\t" + gd + "\t" + relstr;
                            relations.push_back(pz);
                        }
                    }
                }
            }
        }
    }
    gtree->CloseIterator(it);

    return relations;
}

Node* FixSpeciesTree(Node* speciesTree)
{
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf() && n->GetLabel() != "")
        {
            n->SetLabel("i_" + n->GetLabel());
        }
    }
    speciesTree->CloseIterator(it);

    string newick = NewickLex::ToNewickString(speciesTree);

    delete speciesTree;


    return NewickLex::ParseNewickString(newick, true);
}





double GetWeight(unordered_map<string, double> &weights, string gene1, string gene2)
{
    if (gene1 == gene2)
        return 1.0;
    string key = gene1 + ";;" + gene2;

    if (weights.find(key) == weights.end())
    {
        cout<<"WARNING: weight not set for " + gene1 + " " + gene2<<endl;
		throw "Weight not set.";
        return 0.0;
    }

    double weight = weights[key];

    //AVOIDS LOTS OF PROBLEMS
    //if (weight == 0.0)
    //    weight = 0.00000000001;

    return weight;
}











bool cluster_compare_ascending(Cluster* c1, Cluster* c2)
{
    Node* s1 = c1->GetLCASpecies();
    Node* s2 = c2->GetLCASpecies();
    return (s1->HasAncestor(s2));
}


bool cluster_compare_descending(Cluster* c1, Cluster* c2)
{
    Node* s1 = c1->GetLCASpecies();
    Node* s2 = c2->GetLCASpecies();
    return (s2->HasAncestor(s1));
}



//SLOW SORTING BECAUSE STD SORT IS BUGGY
vector<Cluster*> SlowSort(vector<Cluster*> clusters)
{
    vector<Cluster*> sorted;

    while (clusters.size() > 0)
    {
        Cluster* curmin = clusters[0];
        int indexMin = 0;
        for (int c = 1; c < clusters.size(); c++)
        {
            Cluster* cluster = clusters[c];

            Node* smin = curmin->GetLCASpecies();
            Node* sclus = cluster->GetLCASpecies();

            if (sclus->HasAncestor(smin))
            {
                curmin = cluster;
                indexMin = c;
            }
        }

        sorted.push_back(curmin);
        clusters.erase(clusters.begin()+indexMin);
    }

    return sorted;
}


vector<Cluster*> ParseClustersFile(vector<string> filenames,Node* speciesTree)
{

    vector<string> lines;

    for (int i = 0; i < filenames.size(); i++)
    {
        vector<string> l_lines = Util::GetFileLines(filenames[i]);
        for (int l = 0; l < l_lines.size(); l++)
        {
            lines.push_back(l_lines[l]);
        }
    }

    vector<Cluster*> clusters;

    bool ingroup = false;
    vector<string> curGeneNames;
    for (int l = 0; l < lines.size(); l++)
    {
        string line = lines[l];

        if (line == "<GROUP>")
        {
            ingroup = true;
            curGeneNames.clear();
        }
        else
        {
            if (line == "</GROUP>")
            {
                if (ingroup)
                {
                    vector<string> curGeneNamesCopy = curGeneNames; //we make a copy of the vector to pass to the cluster
                    Cluster* cluster = new Cluster(curGeneNamesCopy, speciesTree);
                    clusters.push_back(cluster);
                    curGeneNames.clear();
                    ingroup = false;
                }
            }
            else if (line != "" && ingroup)
            {
                curGeneNames.push_back(line);
            }
        }
    }

    return clusters;
}



//format is one line = NAME1 NAME2 WEIGHT
unordered_map< string, double > ParseWeightsFile(string filename)
{
    unordered_map< string, double > weights;

    vector<string> lines = Util::GetFileLines(filename);

    for (int l = 0; l < lines.size(); l++)
    {
        vector<string> pz = Util::Split(lines[l], " ");

        if (pz.size() >= 3)
        {
            string s1 = pz[0] + ";;" + pz[1];
            string s2 = pz[1] + ";;" + pz[0];

            double w = Util::ToDouble(pz[2]);

            weights[s1] = w;
            weights[s2] = w;
        }
    }

    return weights;
}



Cluster* GetBestParalogiesCluster(Cluster* mycluster, vector<Cluster*> &clusters, unordered_map<string, double> &weights, Node* speciesTree, vector<Cluster*> &forbiddenClusters, int startIndex = 0, int verbose = 0)
{
    int* votes = new int[clusters.size()];
    for (int c = 0; c < clusters.size(); c++)
    {
        votes[c] = 0;
    }

    for (int i = 0; i < mycluster->GetNbGenes(); i++)
    {
        double maxweight = -999;
        int favCluster = startIndex;
        //for (int c = startIndex; c < clusters.size(); c++)
        for (int c = 0; c < clusters.size(); c++)
        {
            Cluster* otherCluster = clusters[c];

            bool isForbidden = (std::find(forbiddenClusters.begin(), forbiddenClusters.end(), otherCluster) != forbiddenClusters.end());

            if (!isForbidden && otherCluster != mycluster)
            {
                for (int j = 0; j < otherCluster->GetNbGenes(); j++)
                {

                    double w = GetWeight(weights, mycluster->GetGeneName(i), otherCluster->GetGeneName(j));
                    if (w > maxweight)
                    {
                        maxweight = w;
                        favCluster = c;
                    }

                }
            }
        }

        if (verbose > 1)
            cout<<mycluster->GetGeneName(i)<<" votes for "<<favCluster<<endl;
        votes[favCluster] = votes[favCluster] + 1;
    }
    Cluster* winnerCluster = NULL;
    int maxIndex = 0;

    if (verbose > 1)
    {
        cout<<"Votes are ";
        for (int i = 0; i < clusters.size(); i++)
        {
            cout<<votes[i]<<" ";
        }
        cout<<endl;
    }

    //for (int c = startIndex; c < clusters.size(); c++)
    for (int c = 0; c < clusters.size(); c++)
    {
        if (votes[c] > votes[maxIndex])
            maxIndex = c;
    }

    if (verbose > 1)
        cout<<"Winner is "<<maxIndex<<endl;
    winnerCluster = clusters[maxIndex];

    delete [] votes;

    return winnerCluster;

}



// The best paralogies cluster is the one that maximizes the average weight between
// genes pairs that are both from species below mycluster->GetLCASpecies().
Cluster* GetBestParalogiesCluster__2(Cluster* mycluster, vector<Cluster*> &clusters, unordered_map<string, double> &weights, Node* speciesTree, int startIndex = 0)
{
    Node* lca_species = mycluster->GetLCASpecies();

    Cluster* bestCluster = NULL;
    double bestAverageWeight = -1;

    for (int c = startIndex; c < clusters.size(); c++)
    {
        Cluster* cluster = clusters[c];

        if (cluster != mycluster)
        {
            double totalWeight = 0.0;
            int nbPairs = 0;

            for (int i = 0; i < mycluster->GetNbGenes(); i++)
            {
                string gene1 = mycluster->GetGeneName(i);

                for (int j = 0; j < cluster->GetNbGenes(); j++)
                {
                    string gene2 = cluster->GetGeneName(j);
                    string sname2 = cluster->GetSpeciesName(gene2);
                    Node* s = speciesTree->GetLeafByLabel(sname2);

                    if (s->HasAncestor(lca_species))
                    {
                        totalWeight += GetWeight(weights, gene1, gene2);
                        nbPairs++;
                    }

                }
            }


            if (nbPairs > 0)
            {
                double avg = totalWeight/((double)nbPairs);

                if (avg > bestAverageWeight)
                {
                    bestAverageWeight = avg;
                    bestCluster = cluster;
                }
            }
        }
    }

    return bestCluster;

}



/*vector<Cluster*> FindBadComponent(vector<Cluster*> &clusters, map<Cluster*, int> favorites)
{

}*/

void DDD(vector<Cluster*> clusters, unordered_map<string, double> &weights, Node* speciesTree, int startIndex = 0, int verbose = 0)
{



    /*map<Cluster*, int> favorites;


    for (int c = 0; c < clusters.size(); c++)
    {
        int fav = GetBestParalogiesCluster(clusters[c], clusters, weights, speciesTree, 0, verbose);
        favorites[c] = fav;
    }

    vector<Cluster*> lowClusters;
    lowClusters.push_back(clusters[0]);
    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* oc = clusters[c];
        if (oc->GetLCASpecies() == clusters[0]->GetLCASpecies())
        {
            lowClusters.push_back(oc);
        }
    }

    vector<Cluster*> badComponent = FindBadComponent(lowClusters, favorites);
    while (badComponent != NULL)
    {
        for (int b = 0; b < badComponent.size(); b++)
        {

        }
    }*/



}


Node* FindBestNodeToGraftAbove(Cluster* c_to_graft, Cluster* c_recv)
{
    Node* species_c_to_graft = c_to_graft->GetLCASpecies();
    GeneTree recv_gene_tree = c_recv->GetGeneTree();
    Node* best_node_to_graft_above = NULL;



    TreeIterator* it = recv_gene_tree.geneTreeRoot->GetPostOrderIterator();
    while (Node* node = it->next())
    {
        if (!node->IsRoot())
        {
            Node* sp_node = recv_gene_tree.GetLCAMapping(node);
            Node* sp_above = recv_gene_tree.GetLCAMapping(node->GetParent());



            if (( sp_node->HasAncestor(species_c_to_graft) || sp_node == species_c_to_graft)
                    && species_c_to_graft->HasAncestor(sp_above) && species_c_to_graft != sp_above)
            {
                best_node_to_graft_above = node;


            }
        }
    }
    recv_gene_tree.geneTreeRoot->CloseIterator(it);

    if (!best_node_to_graft_above)
        best_node_to_graft_above = recv_gene_tree.geneTreeRoot;

   return best_node_to_graft_above;
}



double GetBestWeight(Cluster* c1, Cluster* c2, unordered_map< string, double > &weights)
{
    double max = 0.0;
    for (int i = 0; i < c1->GetNbGenes(); i++)
    {
        for (int j = 0; j < c2->GetNbGenes(); j++)
        {
            double w = GetWeight(weights, c1->GetGeneName(i), c2->GetGeneName(j));

            if (w > max)
                max = w;
        }
    }
    return max;
}


bool HaveDisjointSpecies(Cluster* c1, Cluster* c2)
{
    for (int i = 0; i < c1->GetNbGenes(); i++)
    {
        vector<string> pz1 = Util::Split(c1->GetGeneName(i), g_species_separator);

        for (int j = 0; j < c2->GetNbGenes(); j++)
        {
            vector<string> pz2 = Util::Split(c2->GetGeneName(j), g_species_separator);

            if (pz1[g_species_index] == pz2[g_species_index])
                return false;
        }
    }

    return true;
}

map<Cluster*, Cluster*> ComputeFavoriteMergeableClusters(vector<Cluster*> clusters, unordered_map< string, double > &weights, Node* speciesTree, int verbose)
{
    map<Cluster*, Cluster*> favs;

    for (int i = 0; i < clusters.size(); i++)
    {
        Cluster* c1 = clusters[i];

        vector<Cluster*> forbidden;

        for (int j = 0; j < clusters.size(); j++)
        {
            Cluster* c2 = clusters[j];

            if (!HaveDisjointSpecies(c1, c2))
            {
                forbidden.push_back(c2);
            }

        }

        vector<Cluster*> empty;
        Cluster* favUnrestricted = GetBestParalogiesCluster(c1,clusters, weights, speciesTree, empty, 0, verbose);
        Cluster* favRestricted = GetBestParalogiesCluster(c1,clusters, weights, speciesTree, forbidden, 0, verbose);

        if (favUnrestricted == favRestricted)
        {
            favs[c1] = favUnrestricted;
        }

    }

    return favs;
}


vector<Cluster*> PreMergeClusters(int preMergeType, vector<Cluster*> clusters, unordered_map< string, double > &weights, Node* speciesTree, int verbose = 0)
{
    //TODO: delete hanging clusters
    map<Cluster*, Cluster*> favs = ComputeFavoriteMergeableClusters(clusters, weights, speciesTree, verbose);

    bool done = false;
    while (!done)
    {
        /*Cluster* pair1 = NULL;
        Cluster* pair2 = NULL;

        for (int i = 0; i < clusters.size(); i++)
        {
            if (!pair1)
            {
                Cluster* c1 = clusters[i];
                if (favs.find(c1) != favs.end() && favs[c1] != NULL)
                {
                    Cluster* c2 = favs[c1];

                    if (favs.find(c2) != favs.end() && favs[c2] != NULL)
                    {
                        //a nice pair: merge them
                        if (favs[c2] == c1)
                        {
                            pair1 = c1;
                            pair2 = c2;
                        }
                    }
                }
            }
        }



        if (!pair1)
        {
            done = true;
        }
        else
        {

            if (verbose > 0)
            {
                cout<<"Premerging"<<endl<<pair1->ToString()<<endl<<"with"<<endl<<pair2->ToString()<<endl;
            }
            vector<string> newnames;
            for (int g = 0; g < pair1->GetNbGenes(); g++)
            {
                newnames.push_back(pair1->GetGeneName(g));
            }
            for (int g = 0; g < pair2->GetNbGenes(); g++)
            {
                newnames.push_back(pair2->GetGeneName(g));
            }

            Cluster* newcluster = new Cluster(newnames, speciesTree);

            clusters.erase(std::remove(clusters.begin(), clusters.end(), pair1), clusters.end());
            clusters.erase(std::remove(clusters.begin(), clusters.end(), pair2), clusters.end());
            clusters.push_back(newcluster);
            clusters = SlowSort(clusters);
        }*/

        bool found = false;
        for (int i = 0; i < clusters.size(); i++)
        {
            if (!found)
            {
                Cluster* pair1 = clusters[i];
                if (favs.find(pair1) != favs.end())
                {
                    Cluster* pair2 = favs[pair1];

                    bool ok = true;

                    if (preMergeType == 2)
                    {
                        if (favs.find(pair2) == favs.end())
                        {
                            ok = false;
                        }
                        else
                        {
                            if (favs[pair2] != favs[pair1])
                                ok = false;
                        }
                    }


                    if (ok)
                    {
                        if (verbose > 0)
                        {
                            cout<<"Premerging"<<endl<<pair1->ToString()<<endl<<"with"<<endl<<pair2->ToString()<<endl;
                        }
                        vector<string> newnames;
                        for (int g = 0; g < pair1->GetNbGenes(); g++)
                        {
                            newnames.push_back(pair1->GetGeneName(g));
                        }
                        for (int g = 0; g < pair2->GetNbGenes(); g++)
                        {
                            newnames.push_back(pair2->GetGeneName(g));
                        }

                        Cluster* newcluster = new Cluster(newnames, speciesTree);

                        clusters.erase(std::remove(clusters.begin(), clusters.end(), pair1), clusters.end());
                        clusters.erase(std::remove(clusters.begin(), clusters.end(), pair2), clusters.end());
                        clusters.push_back(newcluster);
                        clusters = SlowSort(clusters);
                        found = true;
                    }
                }

            }
        }

        if (!found)
            done = true;

    }

    return clusters;
}



Cluster* ApplyGraftingOnClusters(Cluster* c1, Cluster* c2, int verbose)
{
    /*SuperGeneTreeMaker* sgt = new SuperGeneTreeMaker();

    Node* t1 = c1->GetGeneTree().geneTreeRoot;
    Node* t2 = c2->GetGeneTree().geneTreeRoot;


    vector<Node*> trees;
    trees.push_back(t1);
    trees.push_back(t2);

    vector<unordered_map<Node*, Node*> > maps;
    maps.push_back(c1->GetGeneTree().lca_mapping);
    maps.push_back(c2->GetGeneTree().lca_mapping);

cout<<"GONNA MERGE"<<endl;
cout<<NewickLex::ToNewickString(t1)<<endl;
cout<<NewickLex::ToNewickString(t2)<<endl;

    pair<Node*, int> res = sgt->GetSuperGeneTreeMinDL(trees, maps, c1->GetGeneTree().speciesTree, true);
cout<<"MERGE DONE"<<endl;

    Node* mergedTree = res.first;

    vector<string> newnames;
    for (int g = 0; g < c1->GetNbGenes(); g++)
    {
        newnames.push_back(c1->GetGeneName(g));
    }
    for (int g = 0; g < c2->GetNbGenes(); g++)
    {
        newnames.push_back(c2->GetGeneName(g));
    }

    //TODO: delte hanging pointer
    Cluster* newcluster = new Cluster(newnames, c1->GetGeneTree().speciesTree);

    newcluster->SetGeneTree(mergedTree, c1->GetGeneTree().speciesTree);

    TreeIterator* it = mergedTree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            unordered_map<Node*,Node*> tmp_map = newcluster->GetGeneTree().lca_mapping; //because c++
            if (GeneSpeciesTreeUtil::Instance()->IsNodeDup(n, tmp_map))
            {
                n->SetCustomField("event", "duplication");
                n->SetLabel("Dup");
            }
            else
            {
                n->SetCustomField("event", "speciation");
                n->SetLabel("Spec");
            }
        }
    }
    mergedTree->CloseIterator(it);

    return newcluster;*/

    Cluster* c_to_graft = c1;
    Cluster* c_recv = c2;
    if (c2->GetLCASpecies()->HasAncestor(c1->GetLCASpecies()))
    {
        c_to_graft = c2;
        c_recv = c1;
    }

    vector<string> newnames;
    for (int g = 0; g < c1->GetNbGenes(); g++)
    {
        newnames.push_back(c1->GetGeneName(g));
    }
    for (int g = 0; g < c2->GetNbGenes(); g++)
    {
        newnames.push_back(c2->GetGeneName(g));
    }

    //TODO: delte hanging pointer
    Cluster* newcluster = new Cluster(newnames, c1->GetGeneTree().speciesTree);




    Node* bestNode = FindBestNodeToGraftAbove(c_to_graft, c_recv);
    Node* newGeneTree = NULL;
    Node* newnode = NULL;

    if (verbose > 1)
    {
        cout<<endl<<endl<<c_to_graft->ToString()<<"( anc="<<c_to_graft->GetLCASpecies()->GetLabel()<<")"<<endl<<" will graft on "<<
              c_recv->ToString()<<"( anc="<<c_recv->GetLCASpecies()->GetLabel()<<")"<<
              endl<<" above "<<bestNode->GetLabel()<<endl;
    }

    if (!bestNode->GetParent())  //graft above root
    {

        newnode = new Node(false);
        newnode->AddSubTree(c_to_graft->GetGeneTree().geneTreeRoot);
        newnode->AddSubTree(c_recv->GetGeneTree().geneTreeRoot);

        newGeneTree = newnode;

    }
    else
    {
        newnode = bestNode->GraftOnParentEdge(c_to_graft->GetGeneTree().geneTreeRoot);
        newGeneTree = c_recv->GetGeneTree().geneTreeRoot;
    }

    newcluster->SetGeneTree(newGeneTree, c1->GetGeneTree().speciesTree);


    if (newcluster->GetGeneTree().IsDupUnderLCAMapping(newnode))
    {
        newnode->SetCustomField("event", "duplication");
        newnode->SetLabel("Dup");
    }
    else
    {
        newnode->SetCustomField("event", "speciation");
        newnode->SetLabel("Spec");
    }

    return newcluster;
}



vector<Cluster*> GetRecursiveInGrafters(Cluster* mycluster, map<Cluster*, vector<Cluster*> > &inGrafters)
{
    vector<Cluster*> out;

    for (int g = 0; g < inGrafters[mycluster].size(); g++)
    {
        out.push_back(inGrafters[mycluster][g]);

        vector<Cluster*> deeperGrafters = GetRecursiveInGrafters(inGrafters[mycluster][g], inGrafters);
        for (int d = 0; d < deeperGrafters.size(); d++)
        {
            out.push_back(deeperGrafters[d]);
        }
    }

    return out;
}


void ComputePctID(string filename, string outfile)
{
	vector<string> lines = Util::Split(Util::GetFileContent(filename), "\n");
	
	vector<string> seqnames;
	vector<string> seqs;
	
	string curseqname = "";
	string curseq = "";
	
	int maxlen = 0;
	
	for (int l = 0; l < lines.size(); l++)
	{
		string line = lines[l];
		line = Util::ReplaceAll(line, "\n", "");
		line = Util::ReplaceAll(line, "\r", "");
		if (line[0] == '>')
		{
			if (curseqname != "")
			{
				seqnames.push_back(curseqname);
				seqs.push_back(curseq);
				if (curseq.length() > maxlen)
			maxlen = curseq.length();
			}
			
			curseqname = Util::ReplaceAll(line, ">", "");
			curseq = "";
		}
		else if (curseqname != "")
		{
			curseq += Util::ReplaceAll(line, "-", "");
		}
	}
	
	
	if (curseqname != "")
	{
		seqnames.push_back(curseqname);
		seqs.push_back(curseq);
		if (curseq.length() > maxlen)
			maxlen = curseq.length();
	}
	
	/*cout<<"Loaded seqs"<<endl;
	for (int s = 0; s < seqnames.size(); s++)
	{
		cout<<seqnames[s]<<endl<<seqs[s]<<endl;
	}*/
	cout<<"Found "<<seqs.size()<<" sequences"<<endl;
	
	double** arr = new double*[maxlen];
	for(int i = 0; i < maxlen; ++i)
		arr[i] = new double[maxlen];
	
	double gapcost = 0.1;
	double mutcost = 0.1;
	
	string outstr = "";
	
	for (int a = 0; a < seqnames.size(); a++)
	{
		cout<<"Calculating scores for "<<seqnames[a]<<" ("<<a<<"/"<<seqnames.size()<<")"<<endl;
		for (int b = a; b < seqnames.size(); b++)
		{
			
			
			string s1 = seqs[a];
			string s2 = seqs[b];
			
			
			
			for (int i = 0; i < s1.length(); i++)
			{
				for (int j = 0; j < s2.length(); j++)
				{
					if (i == 0 || j == 0)
						arr[i][j] = 0;
					else
					{
						int dx = arr[i - 1][j - 1] - mutcost;
						if (s1[i] == s2[j])
							dx = arr[i - 1][j - 1] + 1;
						
						double max = arr[i - 1][j] - gapcost;
						if (arr[i][j - 1] - gapcost > max)
							max = arr[i][j - 1] - gapcost;
						if (dx > max)
							max = dx;
						arr[i][j] = max;
					}
				}
			}
			
			int minlen = s1.length();
			if (s2.length() < minlen)
				minlen = s2.length();
			
			double pctid = (arr[s1.length() - 1][s2.length() - 1]); ///((double)minlen);
			if (outstr != "")
				outstr += "\n";
			outstr += seqnames[a] + ";" + seqnames[b] + ";" + Util::ToString(pctid);
		}
	}
	
	Util::WriteFileContent(outfile, outstr);
	
	for(int i = 0; i < maxlen; ++i)
	{
		delete [] arr[i];
	}
	delete [] arr;
	
}



int main(int argc, char *argv[])
{

    int verbose = 0;

    int preMergeType = 0;

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


    //FOR TESTING PURPOSES
    /*string workdir = "C:/Users/Manuel/Desktop/tmp/";
    workdir = "W:/Users/Manuel/Desktop/tmp/";
    //args["s"] = workdir + "speciestree.nw";
    //args["c"] = workdir + "cplex-3best-di_clusters_08_2.0000000000.clusters";
    //args["w"] = workdir + "alignment_08_2.0000000000_TRUE.dist.edgelist";
    //args["o"] = workdir + "inferred.txt";
    args["m"] = "species_tree";
    args["c"] = workdir + "true_clusters_01.clusters;;" + workdir + "true_clusters_01.clusters;;" + workdir +
                "true_clusters_02.clusters;;" + workdir + "true_clusters_03.clusters;;" +
            workdir + "true_clusters_04.clusters;;" + workdir + "true_clusters_05.clusters;;" + workdir + "true_clusters_06.clusters;;" +
            workdir + "true_clusters_07.clusters;;" + workdir + "true_clusters_08.clusters;;" + workdir + "true_clusters_09.clusters;;";
    verbose = 0;*/



    /*string workdir = "W:/Users/Manuel/Desktop/tmp/";
    args["m"] = "find_clusters";
    args["w"] = workdir + "alignment_04_1.0000000000_TRUE.dist.edgelist";
    args["s"] = workdir + "speciestree_01.nw";*/


	if (args.find("spsep") != args.end())
	{
		g_species_separator = args["spsep"];
	}
	
	if (args.find("spindex") != args.end())
	{
		g_species_index = Util::ToInt(args["spindex"]);
	}
	
	Cluster::speciesSeparator = g_species_separator;
	Cluster::speciesIndex = g_species_index;
	
    vector<Cluster*> clusters;
    Node* speciesTree = NULL;
    unordered_map< string, double > weights;   //format for string key is GENEID1;;GENEID2

    string scontent = "";

    string outfile = "";

	
	if (args.find("m") != args.end())
	{
		if (args["m"] == "pctid")
		{
			ComputePctID(args["f"], args["o"]);
			return 0;
		}
		
		//./OCR/OCR -m getrelations -g ./data/swisstree/ST004/consensus_tree.nhx -s ./data/swisstree/speciestree.nhx
		if (args["m"] == "getrelations")
		{
			scontent = Util::GetFileContent(args["s"]);
			speciesTree = NewickLex::ParseNewickString(scontent, true);

			speciesTree = FixSpeciesTree(speciesTree);
			
			
			string gnewick = Util::GetFileContent(args["g"]);
			Node* geneTree = NewickLex::ParseNewickString(gnewick, false);
			
			unordered_map<Node*, Node*> lca_mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, g_species_separator, g_species_index);
			
			TreeIterator* it = geneTree->GetPostOrderIterator();
			while (Node* n = it->next())
			{
				if (!n->IsLeaf())
				{
					Node* sn = lca_mapping[n];
					Node* sn1 = lca_mapping[n->GetChild(0)];
					Node* sn2 = lca_mapping[n->GetChild(1)];
					
					
					if (sn == sn1 || sn == sn2)
					{
						n->SetLabel( n->GetLabel() + "__Dup" );
						
					}
				}
			}
			geneTree->CloseIterator(it);
			
			string outstring = "";

			vector<string> relations = GetOrthologsParalogs(geneTree);
			
			for (int r = 0; r < relations.size(); r++)
			{
				if (r != 0)
					outstring += ";;";
				outstring += relations[r];
			}
			outstring += "\n";
			
			cout<<outstring<<endl;
			
			delete speciesTree;
			delete geneTree;
			return 0;
		}
	}
	
	
    if (args.find("o") != args.end())
    {
        outfile = args["o"];
    }


    if (args.find("p") != args.end())
    {
        preMergeType = Util::ToInt(args["p"]);
    }

    if (args.find("s") != args.end())
    {
        scontent = Util::GetFileContent(args["s"]);

        if (scontent == "")
        {
            cout<<"SCONTENT IS EMPTY"<<endl;
            throw "scontent empty";
        }
        speciesTree = NewickLex::ParseNewickString(scontent, true);

        speciesTree = FixSpeciesTree(speciesTree);


        if (verbose)
            cout<<"Species tree loaded"<<endl;
    }
    else
    {
        if (args.find("m") == args.end())
        {
            cout<<"We need a species tree (use -s filename)";
            throw "No species tree";
        }
    }

    if (args.find("c") != args.end())
    {
        vector<string> filenames = Util::Split(args["c"], ";;");
        cout<<"Loading clusters file "<<args["c"]<<endl;
        clusters = ParseClustersFile(filenames, speciesTree);

        if (clusters.size() == 0)
        {
            cout<<"WARNING: no clusters were loaded."<<endl;
        }



        if (verbose)
        {
            cout<<"Clusters loaded ["<<clusters.size()<<"]"<<endl;


            //std::sort(clusters.begin(), clusters.end(), cluster_compare_ascending);


            for (int c = 0; c < clusters.size(); c++)
            {
                Cluster* oc = clusters[c];
                for (int i = 0; i < oc->GetNbGenes(); i++)
                {
                    cout<<oc->GetGeneName(i);
                    cout<<" ";
                }
                if (speciesTree)
                {
                    cout<<oc->GetLCASpecies()->GetLabel();
                }
                cout<<endl;
            }
        }
    }









    if (args.find("w") != args.end())
    {
        weights = ParseWeightsFile(args["w"]);

        if (verbose)
            cout<<"Weights loaded"<<endl;
    }






    if (args.find("m") != args.end() && args["m"] == "species_tree")
    {
        SpeciesTreeMaker maker;
        if (verbose > 0)
            maker.verbose = verbose;
        string newick = maker.BuildSpeciesTree(clusters);

        if (outfile == "")
        {
            cout<<newick;
        }
        else
        {
            Util::WriteFileContent(outfile, newick);
        }

        return 0;
    }




    if (args.find("m") != args.end() && args["m"] == "find_clusters")
    {


        vector<string> genes;
        set<string> genes_set;

        for (unordered_map<string, double>::iterator it = weights.begin(); it != weights.end(); it++)
        {
            string names = (*it).first;
            vector<string> names_v = Util::Split(names, ";;", false);
            if (names_v.size() == 2)
            {
                string name = names_v[0];
                if (genes_set.find(name) == genes_set.end())
                {
                    genes.push_back(name);
                    genes_set.insert(name);
                }
            }
        }

        ClusterFinder finder(genes, weights, speciesTree);
        vector<Cluster*> clusters = finder.FindClusters();

        string outstr = "";
        for (int c = 0; c < clusters.size(); c++)
        {
            if (outstr != "")
                outstr += "\n";

            outstr += "<GROUP>\n";

            for (int g = 0; g < clusters[c]->GetNbGenes(); g++)
            {
                outstr += clusters[c]->GetGeneName(g) + "\n";
            }

            outstr += "</GROUP>";
        }

        if (outfile == "")
        {
            cout<<outstr;
        }
        else
        {
            Util::WriteFileContent(outfile, outstr);
        }

        return 0;
    }






    if (verbose)
        cout<<"Clusters sorted"<<endl;


    clusters = SlowSort(clusters);

    if (preMergeType >= 1)
    {
        clusters = PreMergeClusters(preMergeType, clusters, weights, speciesTree, verbose);
    }

    //TODO: hanging POINTERS!
    vector<Cluster*> forbidden;
    while (clusters.size() > 1)
    {
        Cluster* c1 = clusters[0];
        Cluster* c2 = GetBestParalogiesCluster(c1, clusters, weights, speciesTree, forbidden, 0, verbose);

        Cluster* c_merged = ApplyGraftingOnClusters(c1, c2, verbose);
        clusters.erase(std::remove(clusters.begin(), clusters.end(), c1), clusters.end());
        clusters.erase(std::remove(clusters.begin(), clusters.end(), c2), clusters.end());
        clusters.push_back(c_merged);
        clusters = SlowSort(clusters);
    }


    Cluster* lastCluster = clusters[0];
    Node* finalGeneTree = lastCluster->GetGeneTree().geneTreeRoot;
    if (finalGeneTree->GetParent())
    {
        cout<<endl<<endl<<"ERROR: finalGeneTree is not the root.  Why?"<<endl<<endl;
    }

    string newick = NewickLex::ToNewickString(finalGeneTree);

    string outstring = "";

    outstring += "TREE=" + newick + "\n";

    vector<string> relations = GetOrthologsParalogs(finalGeneTree);
    outstring += "RELATIONS=";
    for (int r = 0; r < relations.size(); r++)
    {
        if (r != 0)
            outstring += ";;";
        outstring += relations[r];
    }
    outstring += "\n";


    if (outfile == "")
    {
        cout<<outstring;
    }
    else
    {
        Util::WriteFileContent(outfile, outstring);
    }

    //wrap it up
    /*for (int i = 0; i < clusters.size(); i++)
    {
        delete clusters[i];
    }*/

    delete speciesTree;



    return 0;

    /*map<Cluster*, vector<Cluster*> > inGrafters;
    unordered_map<Cluster*, Node*> bestGraftingSpots;
    unordered_map<Cluster*, Cluster*> bestClusters;

    for (int c = 0; c < clusters.size() - 1; c++)   //-1 because we won't graft last cluster
    {
        Cluster* cluster = clusters[c];


        if (verbose)
            cout<<"Finding best paralogy cluster for "<<cluster->GetGeneName(0)<<endl;

        vector<Cluster*> allForbidden = GetRecursiveInGrafters(cluster, inGrafters);
        Cluster* bestCluster = GetBestParalogiesCluster(cluster, clusters, orthologyWeights, speciesTree, allForbidden, 0, verbose);

        if (!bestCluster)  //This happens if cur cluster is minimal wrt LCAspecies
        {
            bestCluster = clusters[clusters.size() - 1];
        }

        cout<<cluster->GetGeneName(0)<<" --> "<<bestCluster->GetGeneName(0)<<endl;

        //Find best spot to graft on.  We only consider uv such that v <= cluster and u > cluster
        //Node* best_node_to_graft_above = NULL;

        if (verbose)
            cout<<"Building gene tree for "<<bestCluster->GetGeneName(0)<<endl;

        //Add cluster as an in-node to bestCluster
        if (inGrafters.find(bestCluster) == inGrafters.end())
        {
            inGrafters[bestCluster] = vector<Cluster*>();
        }
        vector<Cluster*> vtmp = inGrafters[bestCluster];
        vtmp.push_back(cluster);
        inGrafters[bestCluster] = vtmp;

        //GeneTree bestGeneTree = bestCluster->GetGeneTree();

        //for debug
        if (verbose > 0)
        {
            string tmp = NewickLex::ToNewickString(cluster->GetGeneTree().geneTreeRoot);
            cout<<"CLUSTER "<<c<<endl<<tmp<<endl;
        }




        bestClusters[cluster] = bestCluster;
        //bestGraftingSpots[cluster] = best_node_to_graft_above;  //NULL if root



    }



    //we have the best spots, now we just graft
    //NOTE: the reason we don't do the grafting while computing the best spot is that
    //      the 'best spot' computation would be affected by prior graftings
    for (int c = 0; c < clusters.size() - 1; c++)
    {
        Cluster* cluster = clusters[c];
        GeneTree mygeneTree = cluster->GetGeneTree();

        Cluster* bestCluster = bestClusters[cluster];
        //Node* bestNode = bestGraftingSpots[cluster];
        GeneTree bestGeneTree = bestCluster->GetGeneTree();

        Node* cluster_species = cluster->GetLCASpecies();

        Node* best_node_to_graft_above = NULL;

        //----------------------------------------------------------------------------------------------------------------
        //now find out where to graft it
        TreeIterator* it = bestGeneTree.geneTreeRoot->GetPostOrderIterator();
        while (Node* node = it->next())
        {
            if (!node->IsRoot())
            {
                Node* sp_node = bestGeneTree.GetLCAMapping(node);
                Node* sp_above = bestGeneTree.GetLCAMapping(node->GetParent());

                if (( sp_node->HasAncestor(cluster_species) || sp_node == cluster_species)
                        && cluster_species->HasAncestor(sp_above) && cluster_species != sp_above)
                {
                    best_node_to_graft_above = node;
                }
            }
        }
        bestGeneTree.geneTreeRoot->CloseIterator(it);

        if (!best_node_to_graft_above)
            best_node_to_graft_above = bestGeneTree.geneTreeRoot;


        if (bestCluster)
        {
            if (verbose > 0)
            {
                cout<<cluster->GetGeneName(0)<<" WILL GRAFT ON "<<bestCluster->GetGeneName(0)<<" AT "<<
                      best_node_to_graft_above->GetLabel()<<endl;
                cout<<endl;
            }
        }
        else        //not supposed to happen
        {

            if (verbose > 0)
            {
                cout<<cluster->GetGeneName(0)<<" has no best"<<endl;
            }

            //we default to the last cluster
            bestCluster = clusters[clusters.size() - 1];
            best_node_to_graft_above = bestCluster->GetGeneTree().geneTreeRoot;
        }
        //----------------------------------------------------------------



        Node* bestNode = best_node_to_graft_above;
        Node* newnode = NULL;
        if (!bestNode->GetParent())  //graft above root of bestCluster
        {



            if (verbose > 0)
                cout<<"c="<<c<<" "<<cluster->GetGeneName(0)<<" grafting above "<<bestCluster->GetGeneName(0)<<endl;

            newnode = new Node(false);
            newnode->AddSubTree(mygeneTree.geneTreeRoot);
            newnode->AddSubTree(bestGeneTree.geneTreeRoot);
			

			
            cluster->SetGeneTreeRoot(newnode);
            bestCluster->SetGeneTreeRoot(newnode);

        }
        else
        {
            newnode = bestNode->GraftOnParentEdge(mygeneTree.geneTreeRoot);
        }

        cluster->UpdateLCAMapping();
        bestCluster->UpdateLCAMapping();

        if (bestCluster->GetGeneTree().IsDupUnderLCAMapping(newnode))
        {
            newnode->SetCustomField("event", "duplication");
            newnode->SetLabel("Dup");
        }
        else
        {
            newnode->SetCustomField("event", "speciation");
            newnode->SetLabel("Spec");
        }

    }

    Cluster* lastCluster = clusters[clusters.size() - 1];
    Node* finalGeneTree = lastCluster->GetGeneTree().geneTreeRoot;
    if (finalGeneTree->GetParent())
    {
        cout<<endl<<endl<<"ERROR: finalGeneTree is not the root.  Why?"<<endl<<endl;
    }

    string newick = NewickLex::ToNewickString(finalGeneTree);

    string outstring = "";

    outstring += "TREE=" + newick + "\n";

    vector<string> relations = GetOrthologsParalogs(finalGeneTree);
    outstring += "RELATIONS=";
    for (int r = 0; r < relations.size(); r++)
    {
        if (r != 0)
            outstring += ";;";
        outstring += relations[r];
    }
    outstring += "\n";


    if (outfile == "")
    {
        cout<<outstring;
    }
    else
    {
        Util::WriteFileContent(outfile, outstring);
    }

    //wrap it up
    for (int i = 0; i < clusters.size(); i++)
    {
        delete clusters[i];
    }

    delete speciesTree;



    return 0;*/

}

