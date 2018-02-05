#include "speciestreemaker.h"




SpeciesTreeMaker::SpeciesTreeMaker()
{
    verbose = 0;
    applyCladesFilter = true;
}



void SpeciesTreeMaker::BuildPossibleClades(vector<Cluster*> &clusters, set<string> &allSpecies)
{
    //TODO: make all this faster
    possibleClades.clear();

    if (verbose > 0)
        cout<<"Adding direct clades"<<endl;

    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];

        set<string> clsp = cl->GetSpeciesSet();

        if (!IsPossibleClade(clsp))
        {
            possibleClades.insert(clsp);
            set<string> compl_clsp;
            compl_clsp = Util::GetSetComplement(clsp, allSpecies);
            possibleClades.insert(compl_clsp);
        }
    }

    if (verbose > 0)
        cout<<"Adding order clades (done part 0/3)"<<endl;


    //experimental stuff: order the species by order of appearance (s1,s2,...,sk), and
    //any contiguous interval is a possible clade
    //This includes, in particular, any single species.
    vector< pair<string, int> > speciesCounts;
    for (set<string>::iterator it = allSpecies.begin(); it != allSpecies.end(); it++)
    {
        string sp = (*it);
        //count how many have sp
        int cnt = 0;
        for (int c = 0; c < clusters.size(); c++)
        {
            Cluster* cl = clusters[c];
            if (cl->HasSpecies(sp))
            {
                cnt++;
            }
        }
        pair<string, int> p;
        p.first = sp;
        p.second = cnt;
        speciesCounts.push_back(p);
    }

    if (verbose > 0)
        cout<<"Adding order clades (done part 1/3)"<<endl;

    std::sort(speciesCounts.begin(), speciesCounts.end(), compare_pair_second<>());

    if (verbose > 0)
        cout<<"Adding order clades (done part 2/3)"<<endl;


    for (int i = 0; i < speciesCounts.size(); i++)
    {
        for (int j = i; j < speciesCounts.size(); j++)
        {
            set<string> curset;

            for (int k = i; k <= j; k++)
            {
                curset.insert(speciesCounts[k].first);
            }

            if (!IsPossibleClade(curset))
            {
                possibleClades.insert(curset);
                set<string> compl_curset = Util::GetSetComplement(curset, allSpecies);
                possibleClades.insert(compl_curset);
            }

        }
    }

    if (verbose > 0)
        cout<<"Done 3/3"<<endl;

}




bool SpeciesTreeMaker::IsPossibleClade(set<string> clade)
{
    string key = GetSpeciesSetKey(clade);
    //TODO: could be faster
    for (set< set<string> >::iterator it = possibleClades.begin(); it != possibleClades.end(); it++)
    {
        set<string> clade = (*it);

        if (key == GetSpeciesSetKey(clade))
            return true;
    }

    return false;
}


string SpeciesTreeMaker::BuildSpeciesTree(vector<Cluster*> &clusters)
{
    //return val = newick of best species tree

    this->speciesSetCache.clear();

    map<string, set<string> > blocks = Cluster::GetBlocks(clusters);


    if (verbose > 0)
    {
        cout<<"THE BLOCKS="<<endl;
        for (map<string, set<string> >::iterator it = blocks.begin(); it != blocks.end(); it++)
        {
            cout<<Util::ToString((*it).second)<<endl;
        }
		cout<<"/BLOCKS"<<endl;
    }

    vector<Cluster*> compressed_clusters;




    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];
        string new_gene_name = Util::ToString(c);
        Cluster* compressed = cl->GetCompressedCluster(blocks, new_gene_name);
        compressed_clusters.push_back(compressed);

        if (verbose > 0)
        {
            cout<<compressed->ToString()<<endl<<endl<<endl;
        }
    }

    //build set of all species
    set<string> species = Cluster::GetAllSpecies(compressed_clusters);

    if (verbose > 0)
    {
        cout<<"Species set"<<endl;
        for (set<string>::iterator it = species.begin(); it != species.end(); it++)
        {
            cout<<(*it)<<endl;
        }
        cout<<endl;
        cout<<"NBSPECIES="<<species.size()<<endl;
    }


    if (this->applyCladesFilter)
    {
        cout<<"BUILDING POSSIBLE CLADES"<<endl;
        BuildPossibleClades(compressed_clusters, species);
        if (verbose > 0)
        {
            for (set< set<string> >::iterator it = possibleClades.begin(); it != possibleClades.end(); it++)
            {
                cout<<GetSpeciesSetKey(*it)<<endl;
            }
        }
        cout<<"DONE, "<<possibleClades.size()<<" CLADES"<<endl;

    }


    pair<int, string> vals = EvaluateSpeciesSet(species, compressed_clusters);


    for (int c = 0; c < compressed_clusters.size(); c++)
    {
        delete compressed_clusters[c];
    }

    return vals.second + ";";
}


/*vector< pair< set<string>, set<string> > > SpeciesTreeMaker::GetAllSplits(set<string> &species)
{

}*/


string SpeciesTreeMaker::GetSpeciesSetKey(set<string> species)
{
    string key = "";
    for (set<string>::iterator it = species.begin(); it != species.end(); it++)
    {
        if ( key != "")
            key += ";;;;";
        key += (*it);
    }
    return key;
}

set<string> SpeciesTreeMaker::GetKeySpecies(string key)
{
    vector<string> v = Util::Split(key, ";;;;");
    set<string> ret;
    ret.insert(v.begin(), v.end());

    return ret;
}


//return val = <cost, newick>
pair<int, string> SpeciesTreeMaker::EvaluateSpeciesSet(set<string> &species, vector<Cluster*> &clusters)
{

    string speciesSetKey = this->GetSpeciesSetKey(species);
    if (this->speciesSetCache.find( speciesSetKey ) != this->speciesSetCache.end() )
    {
        //cout<<"USING CACHE FOR "<<speciesSetKey<<endl;

        return this->speciesSetCache[speciesSetKey];
    }

    vector<string> species_vector;  //we'll need a guaranteed fixed order later on
    for (set<string>::iterator it = species.begin(); it != species.end(); it++)
    {
        species_vector.push_back(*it);
    }

    //compute intersection between cluster and species.  intersections[i] for clusters[i]
    //allSupersetOfSpecies is true if each cluster is either a superset of species, or has none in common
    vector< set<string> > intersections;
    bool allSupersetOfSpecies = true;
    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];
        set<string> inter = cl->GetSpeciesIntersection(species);

        intersections.push_back(inter);

        if (inter.size() > 0 && inter.size() < species.size())
            allSupersetOfSpecies = false;

    }

    if (allSupersetOfSpecies)
    {
        //any tree will do, so of course we make a caterpillar
        string newick = NewickLex::GetCaterpillarNewick(species);
        pair<int, string> retval;
        retval.first = 0;
        retval.second = newick;

        this->speciesSetCache[speciesSetKey] = retval;

        return retval;
    }
    else
    {
        int min = 99999999;
        string minnewick = "";



        if (this->applyCladesFilter && species.size() > 4)  //even if there's a filter, we allow trying every 4-set since it's not that long
        {
            int nbtried = 0;
            for (set< set<string> >::iterator it = possibleClades.begin(); it != possibleClades.end(); it++)
            {
                set<string> set1 = (*it);

                if (set1.size() < species.size() && set1.size() > 0 && Util::SetContains(species, set1))
                {
                    for (set< set<string> >::iterator it2 = it; it2 != possibleClades.end(); it2++)
                    {

                        set<string> set2 = (*it2);

                        if (Util::SetContains(species, set2) && Util::GetSetIntersection(set1, set2).size() == 0 && (set1.size() + set2.size() == species.size()))
                        {

                            /*cout<<"SP="<<Util::ToString(species)<<endl
                                <<"S1="<<Util::ToString(set1)<<endl
                                <<"S2="<<Util::ToString(set2)<<endl<<endl;*/
                            pair<int, string> vals = EvaluateSplit( set1, set2, clusters, intersections );

                            if (vals.first < min)
                            {
                                min = vals.first;
                                minnewick = vals.second;
                            }

                            nbtried++;

                        }

                    }
                }


                /*set<string> s1 = Util::GetSetIntersection(*it, species);
                set<string> s2 = Util::GetSetComplement(s1, species);

                if (s1.size() > 0 && s2.size() > 0)
                {
                    pair<int, string> vals = EvaluateSplit( s1, s2, clusters, intersections );

                    if (vals.first < min)
                    {
                        min = vals.first;
                        minnewick = vals.second;
                    }
                }*/
            }


            //TODO : check if not already tried
            set<string> clade = GetMostCommonSubclade(clusters, species);
            set<string> complement = Util::GetSetComplement(clade, species);

            pair<int, string> vals = EvaluateSplit( clade, complement, clusters, intersections );

            if (vals.first < min)
            {
                min = vals.first;
                minnewick = vals.second;
            }

        }
        else
        {
            //try every possible split
            LongCounter counter(species.size(), 2);
            counter.Increment();    //skip all zeros
            bool metLast = false;




            int cpt = 0;
            while (!metLast)
            {

                if (counter.GetFromLeft(0) == 1)
                    metLast = true;

                //cout<<"Counter = "<<counter.ToString()<<endl;
                set<string> s1, s2;

                for (int i = 0; i < species_vector.size(); i++)
                {
                    if (counter.GetFromRight(i) == 0)
                    {
                        s1.insert(species_vector[i]);
                    }
                    else
                    {
                        s2.insert(species_vector[i]);
                    }
                }


                pair<int, string> vals = EvaluateSplit( s1, s2, clusters, intersections );

                if (vals.first < min)
                {
                    min = vals.first;
                    minnewick = vals.second;
                }

                counter.Increment();
                cpt++;
            }
        }




        pair<int, string> retval;
        retval.first = min;
        retval.second = minnewick;

        this->speciesSetCache[speciesSetKey] = retval;

        return retval;
    }
}


set<string> SpeciesTreeMaker::GetMostCommonSubclade(vector<Cluster*> &clusters, set<string> &species)
{
    map<string, int> counts;

    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];

        set<string> inter = cl->GetSpeciesIntersection(species);

        if (inter.size() < species.size())
        {
            string key = GetSpeciesSetKey(inter);

            if (counts.find(key) == counts.end())
                counts[key] = 1;
            else
                counts[key] = counts[key] + 1;
        }
    }

    set<string> ret;


    if (counts.size() == 0)  //return anything
    {
        ret.insert(*species.begin());
    }
    else
    {
        string maxstr = "";
        int curmax = 0;
        for (map<string, int>::iterator it = counts.begin(); it != counts.end(); it++)
        {
            if ((*it).second > curmax)
                maxstr = (*it).first;
        }

        ret = GetKeySpecies(maxstr);
    }

    return ret;
}


//return val = <cost, newick>
pair<int, string> SpeciesTreeMaker::EvaluateSplit(set<string> &s1, set<string> &s2, vector<Cluster*> &clusters, vector<set<string> > &intersections)
{

    if (verbose > 0)
    {
        cout<<"EVAL SPLIT "<<endl;
        cout<<Util::ToString(s1)<<endl<<Util::ToString(s2);
        cout<<endl;
    }

    pair<int, string> x1 = EvaluateSpeciesSet(s1, clusters);
    pair<int, string> x2 = EvaluateSpeciesSet(s2, clusters);


    int costWithin = 0;
    for (int c = 0; c < clusters.size(); c++)
    {
        Cluster* cl = clusters[c];

        set<string> inter1 = cl->GetSpeciesIntersection(s1);
        set<string> inter2 = cl->GetSpeciesIntersection(s2);

        //if cl doesn't have anything outside, it won't lose anything
        //so here we check that cl has something outside
        if (inter1.size() + inter2.size() < cl->GetSpeciesSet().size())
        {
            //a loss occurs when cl misses all of s1 or all of s2 but not both
            if ((inter1.size() == 0 || inter2.size() == 0) && inter1.size() + inter2.size() > 0)
            {
                costWithin += 1;    //python increment

                if (verbose > 1)
                {
                    cout<<"MUST PAY:"<<endl<<cl->ToString()<<endl;
                }
            }
        }
    }

    if (verbose > 1)
    {
        cout<<"COST WITHIN="<<costWithin<<endl<<endl;
    }


    pair<int, string> retval;
    retval.first = x1.first + x2.first + costWithin;
    retval.second = "(" + x1.second + "," + x2.second + ")";
    return retval;
}
