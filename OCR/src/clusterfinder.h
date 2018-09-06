#ifndef CLUSTERFINDER_H
#define CLUSTERFINDER_H


#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <unordered_map>
#include "cluster.h"
#include "div/longcounter.h"

using namespace std;


class ClusterFinder
{
public:
    ClusterFinder(vector<string> &_genes, unordered_map<string, double> &_weights, Node* _speciesTree);
    ~ClusterFinder();

    vector<Cluster*> FindClusters();
    //vector<Cluster*> FindClustersWithoutSpeciesTree();

    vector<int> GetClustersWithBestHit(string gene_name, string species, vector<Cluster *> cls);
    vector<int> GetClustersWithMaxAverage(Cluster* cl, vector<Cluster*> cls);
    vector<int> GetClustersWithMaxScore(Cluster *cl, vector<Cluster*> cls);

    vector<int> GetClustersWithBestOutsiderAgreement(Cluster* cl, vector<Cluster*> cls, Node* sp_of_cl);

    double GetMaxScore(Cluster* cl1, Cluster* cl2);
    double GetAverageWeight(Cluster* c1, Cluster* c2);

    double GetWeight(string gene1, string gene2);



    vector<Cluster*> GetClustersUnder(Node* sp);

    int GetBestHitCluster(Cluster* cl, vector<Cluster*> cls, Node *sp_of_cl, Node *sp_out);

private:
    vector<Cluster*> clusters;
    vector<string> genes;
    unordered_map< string, double > weights;
    Node* speciesTree;

    void MergeClusters(Cluster* cl1, Cluster* cl2);

};

#endif // CLUSTERFINDER_H
