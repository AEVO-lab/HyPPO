#include "clusterfinder.h"

ClusterFinder::ClusterFinder(vector<string> &_genes, unordered_map< string, double > &_weights, Node* _speciesTree)
{

    this->genes = _genes;
    this->weights = _weights;
    this->speciesTree = _speciesTree;

    //we work with clusters in this algo.  To start with, each gene is a cluster, and they grow.
    for (int g = 0; g < genes.size(); g++)
    {
        vector<string> v;
        v.push_back(genes[g]);
        Cluster* cl = new Cluster(v, speciesTree);
        this->clusters.push_back(cl);
    }
}


ClusterFinder::~ClusterFinder()
{
    for (int c = 0; c < clusters.size(); c++)
    {
        delete clusters[c];
    }
    clusters.clear();
}

double ClusterFinder::GetWeight(string gene1, string gene2)
{
    if (gene1 == gene2)
        return 1.0;
    string key = gene1 + ";;" + gene2;

    if (weights.find(key) == weights.end())
    {
        cout<<"WARNING: weight not set for " + gene1 + " " + gene2<<endl;
        return 0.0;
    }

    double weight = weights[key];

    return weight;
}


//returns list of cluster indexes in cls that have a best hit with gene_name in species.  A vector is returned in case of ties.
vector<int> ClusterFinder::GetClustersWithBestHit(string gene_name, string species, vector<Cluster*> cls)
{
    vector<int> ret;
    double maxWeight = 0;

    for (int c = 0; c < cls.size(); c++)
    {
        Cluster* cl = cls[c];

        bool hasGene = false;	//TODO: not used
        for (int g = 0; g < cl->GetNbGenes(); g++)
        {
            if (cl->GetGeneName(g) == gene_name)
            {
                hasGene = true;
            }
            else
            {
                if (cl->GetSpeciesName(cl->GetGeneName(g)) == species)
                {
                    double w = GetWeight(cl->GetGeneName(g), gene_name);
                    if (w > maxWeight)
                    {
                        ret.clear();
                        ret.push_back(c);
                        maxWeight = w;
                    }
                    else if (w == maxWeight)
                    {
                        ret.push_back(c);
                    }
                }
            }
        }
    }

    return ret;
}


//Let maxScore be the max score of a gene in cl and a gene outside.  Returns all clusters that achieve this maxScore.
vector<int> ClusterFinder::GetClustersWithMaxScore(Cluster* cl, vector<Cluster*> cls)
{
    double bestscore = 0.0;
    vector<int> ret;

    for (int c = 0; c < cls.size(); c++)
    {
        Cluster* cl2 = cls[c];

        double w = GetMaxScore(cl, cl2);


        if (w > bestscore)
        {
            bestscore = w;
            ret.clear();
            ret.push_back(c);
        }
        else if (w == bestscore)
        {
            ret.push_back(c);
        }

    }

    return ret;
}



double ClusterFinder::GetMaxScore(Cluster* cl1, Cluster* cl2)
{
    double max = 0;
    for (int g1 = 0; g1 < cl1->GetNbGenes(); g1++)
    {
        for (int g2 = 0; g2 < cl2->GetNbGenes(); g2++)
        {
            double w = this->GetWeight(cl1->GetGeneName(g1), cl2->GetGeneName(g2));

            if (w > max)
            {
                max = w;
            }
        }
    }

    return max;
}



vector<int> ClusterFinder::GetClustersWithBestOutsiderAgreement(Cluster* cl, vector<Cluster*> cls, Node* sp_of_cl)
{
    vector<int> ret;
    if (sp_of_cl->GetParent()->IsRoot() || sp_of_cl->GetParent()->GetParent()->IsRoot())
    {
        return ret;
    }
    else
    {
        Node* s_gp = sp_of_cl->GetParent()->GetParent();
        Node* sibl = s_gp->GetChild(0);
        if (sibl == sp_of_cl->GetParent())
            sibl = s_gp->GetChild(1);

        vector<Cluster*> cl_sibl = this->GetClustersUnder(sibl);

        /*TreeIterator* it = sibl->GetPostOrderIterator(true);
        while (Node* s = it->next())
        {
            string s_lbl = s->GetLabel();


        }
        sibl->CloseIterator(it);*/

        vector<double> averages;
        for (int c = 0; c < cl_sibl.size(); c++)
        {
            averages.push_back(GetAverageWeight(cl, cl_sibl[c]));
        }

        double min_sum_of_diff = 9999.0;

        for (int c = 0; c < cls.size(); c++)
        {
            double sum_of_diff = 0.0;

            for (int d = 0; d < cl_sibl.size(); d++)
            {
                double avg = GetAverageWeight(cls[c], cl_sibl[d]);
                sum_of_diff += abs(averages[d] - avg);
            }

            if (sum_of_diff < min_sum_of_diff)
            {
                min_sum_of_diff = sum_of_diff;
                ret.clear();
                ret.push_back(c);
            }
            else if (sum_of_diff == min_sum_of_diff)
            {
                ret.push_back(c);
            }

        }

    }

    return ret;
}




double ClusterFinder::GetAverageWeight(Cluster* c1, Cluster* c2)
{
    double ttl = 0.0;
    for (int g1 = 0; g1 < c1->GetNbGenes(); g1++)
    {
        for (int g2 = 0; g2 < c2->GetNbGenes(); g2++)
        {
            ttl += this->GetWeight(c1->GetGeneName(g1), c2->GetGeneName(g2));
        }
    }

    double avg = ttl / double(c1->GetNbGenes() * c2->GetNbGenes());

    return avg;
}


vector<int> ClusterFinder::GetClustersWithMaxAverage(Cluster* cl, vector<Cluster*> cls)
{
    double bestavg = 0.0;
    vector<int> ret;

    for (int c = 0; c < cls.size(); c++)
    {
        Cluster* cl2 = cls[c];

        double avg = GetAverageWeight(cl, cl2);

        if (avg > bestavg)
        {
            bestavg = avg;
            ret.clear();
            ret.push_back(c);
        }
        else if (avg == bestavg)
        {
            ret.push_back(c);
        }
    }

    return ret;
}


int ClusterFinder::GetBestHitCluster(Cluster* cl, vector<Cluster*> cls, Node* sp_of_cl, Node* sp_out)
{


    vector<int> v1 = this->GetClustersWithMaxScore(cl, cls);
    vector<int> v2 = this->GetClustersWithMaxAverage(cl, cls);


    vector<int> inter = Util::GetVectorIntersection(v1, v2);
    if (inter.size() == 1)
    {
        return inter[0];
    }
    else
    {
        /*if (sp_of_cl->GetNbLeaves() <= 3 || sp_out->GetNbLeaves() <= 3)
            return v1[0];
        else
            return v2[0];*/

        vector<int> v3 = GetClustersWithBestOutsiderAgreement(cl, cls, sp_of_cl);

        vector<int> inter_inter = Util::GetVectorIntersection(inter, v3);

        if (inter_inter.size() > 0)
            return inter_inter[0];

        vector<int> inter1 = Util::GetVectorIntersection(v1, v3);
        if (inter1.size() > 0)
            return inter1[0];

        vector<int> inter2 = Util::GetVectorIntersection(v2, v3);
        if (inter2.size() > 0)
            return inter2[0];

        //if all else fails
        //cout<<"NO CLEAR BEST CLUSTER FOUND!"<<endl;

        if (sp_of_cl->GetNbLeaves() <= 3 || sp_out->GetNbLeaves() <= 3)
            return v1[0];
        else
            return v2[0];
        return v1[0];
    }



    /*vector<int> counts(cls.size());
    vector<int> maxpossible(cls.size());
    for (int i = 0; i < counts.size(); i++)
    {
        counts[i] = 0;
        maxpossible[i] = 0;
    }

    TreeIterator* it = sp_out->GetPostOrderIterator(true);
    while (Node* s = it->next())
    {
        //each gene casts a vote on its preferred cluster for s
        //majority wins
        vector<int> s_counts(cls.size());

        for (int i = 0; i < s_counts.size(); i++)
            s_counts[i] = 0; //not sure if needed

        for (int g = 0; g < cl->GetNbGenes(); g++)
        {
            int c_index = GetClusterWithBestHit(cl->GetGeneName(g), s->GetLabel(), cls);
            s_counts[c_index] = s_counts[c_index] + 1;
        }

        pair<int, int> p = Util::GetMaxInVector(s_counts);

        int best_cluster_for_s = p.first;   //this line is for clarity
        counts[best_cluster_for_s] = counts[best_cluster_for_s] + 1;

        //each cluster containing s had an opportunity to be the best.
        //here we record those that did
        for (int c = 0; c < cls.size(); c++)
        {
            if (cls[c]->HasSpecies(s->GetLabel()))
            {
                maxpossible[c] = maxpossible[c] + 1;
            }
        }
    }
    sp_out->CloseIterator(it);

    int curbest = -1;
    double curmax = -1;
    for (int c = 0; c < cls.size(); c++)
    {
        double score = (double)counts[c] / (double)maxpossible[c];
        if (score > curmax)
        {
            curbest = c;
            curmax = score;
        }
    }


    return curbest;*/

}


vector<Cluster*> ClusterFinder::GetClustersUnder(Node* sp)
{
    vector<Cluster*> ret;
    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];

        if (cl->GetLCASpecies()->HasAncestor(sp))
        {
            ret.push_back(cl);
        }
    }
    return ret;
}


vector<Cluster*> ClusterFinder::FindClusters()
{
    if (!this->speciesTree)
    {
        cout<<"Error: no species tree provided."<<endl;
        throw "No species tree provided";
    }



    TreeIterator* it = speciesTree->GetPostOrderIterator();

    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            Node* n1 = n->GetChild(0);
            Node* n2 = n->GetChild(1);

            vector<Cluster*> cls1 = GetClustersUnder(n1);
            vector<Cluster*> cls2 = GetClustersUnder(n2);

            if (cls1.size() > 0 && cls2.size() > 0)
            {
                vector<pair< Cluster*, Cluster* > > bbh_list;

                //find BBH's
                for (int c1 = 0; c1 < cls1.size(); c1++)
                {
                    int best_hit = GetBestHitCluster(cls1[c1], cls2, n1, n2);

                    if (best_hit >= 0)
                    {
                        int rev_best_hit = GetBestHitCluster(cls2[best_hit], cls1, n2, n1);

                        //we have a bbh!
                        if (rev_best_hit == c1)
                        {
                            pair<Cluster*, Cluster*> p;
                            p.first = cls1[c1];
                            p.second = cls2[best_hit];
                            bbh_list.push_back(p);
                        }
                    }
                    else
                    {
                        cout<<"WARNING: no best hit found."<<endl;
                    }
                }


                //handle bbh's
                for (int b = 0; b < bbh_list.size(); b++)
                {
                    Cluster* cl1 = bbh_list[b].first;
                    Cluster* cl2 = bbh_list[b].second;

                    MergeClusters(cl1, cl2);
                }

            }
        }
    }
    speciesTree->CloseIterator(it);

    return clusters;
}



void ClusterFinder::MergeClusters(Cluster* cl1, Cluster* cl2)
{
    //merge them by creating a new cluster and deleting the old ones
    vector<string> new_genes;
    for (int g = 0; g < cl1->GetNbGenes(); g++)
        new_genes.push_back(cl1->GetGeneName(g));
    for (int g = 0; g < cl2->GetNbGenes(); g++)
        new_genes.push_back(cl2->GetGeneName(g));
    Cluster* new_cluster = new Cluster(new_genes, speciesTree);
    clusters.erase( std::remove( clusters.begin(), clusters.end(), cl1 ), clusters.end() );
    clusters.erase( std::remove( clusters.begin(), clusters.end(), cl2 ), clusters.end() );
    delete cl1;
    delete cl2;
    clusters.push_back(new_cluster);
}



/*vector<Cluster*> ClusterFinder::FindClustersWithoutSpeciesTree()
{
    bool done = false;
    while (!done)
    {
        double max = -1;
        int max_c = 0;
        int max_d = 0;

        for (int c = 0; c < clusters.size(); c++)
        {
            for (int d = c + 1; d < clusters.size(); d++)
            {

                vector<string> inter = clusters[c]->GetSpeciesIntersection(clusters[d]->GetSpeciesSet());

                if (inter.size() == 0)
                {
                    double w = GetMaxScore(clusters[c], clusters[d]);
                }
            }
        }

        if (max < 0)
            done = true;
        else
        {

        }
    }

}*/
