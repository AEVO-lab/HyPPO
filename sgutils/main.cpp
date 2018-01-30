#include "trees/node.h"
#include "trees/autester.h"

#include <string>
#include <map>
#include "div/util.h"
#include "div/tinydir.h"
#include "trees/newicklex.h"
#include "trees/node.h"
#include "trees/genespeciestreeutil.h"
#include <cstdio>
#include <stdlib.h>

using namespace std;


//assumes format [genename]__[clusterid]__[whatever]
string GetTrueClustersString(Node* gtree)
{
    unordered_map<string, string > clusters;

    vector<Node*> leaves = gtree->GetLeafVector();

    for (int i = 0; i < leaves.size(); i++)
    {
        Node* leaf = leaves[i];

        string lbl = leaf->GetLabel();
        vector<string> pz = Util::Split(lbl, "__");

        string locusid = pz[1];

        if (clusters.find(locusid) == clusters.end())
        {
            clusters[locusid] = lbl;
        }
        else
        {
            clusters[locusid] = clusters[locusid] + " " + lbl;
        }
    }

    int nb = 0;
    string outstr = "";
    for (unordered_map<string, string>::iterator it = clusters.begin(); it != clusters.end(); it++)
    {
        if (nb != 0)
            outstr += ";;";
        string key = it->first;

        outstr += clusters[key];

        nb++;
    }

    return outstr;
}




void DoRandomNNI(Node* tree, float sizeThreshold)
{
	vector<Node*> nodes = tree->GetPostOrderedNodes();
	
	int nbLeaves = tree->GetNbLeaves();
	
	bool done = false;
	
	while (!done)
	{
		int r = rand() % nodes.size();
		int ch1 = rand() % 2;
		int ch2 = 1 - ch1;
		
		Node* n = nodes[r];
		if (!n->IsLeaf() && ((float)(n->GetNbLeaves())/(float)(nbLeaves) >= sizeThreshold ))
		{
			Node* x = n->GetChild(ch1);
			if (!x->IsLeaf())
			{
				//Node* a = x->GetChild(0);
				Node* b = x->GetChild(1);
				Node* c = n->GetChild(ch2);
				
				x->RemoveChild(b);
				n->RemoveChild(c);
				n->AddSubTree(b);
				x->AddSubTree(c);
				done = true;
			}
		}
		
	}
}


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




string EditINDELibleControlFile(string control_filename, vector< vector<Node*> > multipliedGeneTrees, vector<double> multiples)
{
    string outstring = "";

    vector<string> lines = Util::GetFileLines(control_filename);


    string treestr = "", partitionsstr = "", evolvestr = "";

    for (int locus = 0; locus < multipliedGeneTrees.size(); locus++)
    {
        vector<Node*> gtrees = multipliedGeneTrees[locus];
        for (int multiple_index = 0; multiple_index < gtrees.size(); multiple_index++)
        {
            Node* gtree = gtrees[multiple_index];
            double multiple = multiples[multiple_index];


            string id = (locus + 1 < 10 ? "0" : "") + Util::ToString(locus + 1) + "_" + Util::ToString(multiple);

            treestr += "[TREE] T" + id + " " +
                        NewickLex::ToNewickString(gtree, true, false) + "\n";

            partitionsstr += "[PARTITIONS] simple_" + id +
                                " [T" + id + " simple_unlinked_" + Util::ToString(locus + 1) + " 1000]" + "\n";

            evolvestr += "simple_" + id + " 1 alignment_" + id + "\n";

        }
    }


    bool stopCopying = false;
    for (int l = 0; l < lines.size(); l++)
    {
        //TODO: clean this up
        string line = lines[l];
        if (Util::StartsWith(line, "[TREE]"))
        {
            stopCopying = true;
        }
        else if (Util::StartsWith(line, "[PARTITIONS]"))
        {
            stopCopying = true;
        }
        else if (Util::StartsWith(line, "[EVOLVE]"))
        {
            stopCopying = true;
        }
        else
        {
            if (!stopCopying)
                outstring += line + "\n";
        }
    }

    outstring += treestr + "\n" + partitionsstr + "\n" +
                "[EVOLVE] " + evolvestr + "\n";

    return outstring;
}


int FindNewLocusChild(Node* gtree, Node* dupnode)
{
    //if g is the new locus, all locus index under g aren't found outside of g
    vector<Node*> all_leaves = gtree->GetLeafVector();
    set<string> left_leaves = dupnode->GetChild(0)->GetLeafLabels();
    set<string> right_leaves = dupnode->GetChild(1)->GetLeafLabels();

    set<string> loci_outside;

    for (int i = 0; i < all_leaves.size(); i++)
    {
        Node* leaf = all_leaves[i];
        if (left_leaves.find(leaf->GetLabel()) == left_leaves.end() && right_leaves.find(leaf->GetLabel()) == right_leaves.end())
        {
            vector<string> pz = Util::Split(leaf->GetLabel(), "__");
            string locus = pz[1];

            if (loci_outside.find(locus) == loci_outside.end())
                loci_outside.insert(locus);
        }
    }

    vector<Node*> v_left = dupnode->GetChild(0)->GetLeafVector();
    for (int i = 0; i < v_left.size(); i++)
    {
        Node* leaf = v_left[i];
        vector<string> pz = Util::Split(leaf->GetLabel(), "__");
        string locus = pz[1];

        //locus creation must've occurred at other child
        if (loci_outside.find(locus) != loci_outside.end())
            return 1;
    }

    //not found outside, child 0 must be new locus
    return 0;
}



void ExtractOrthologsParalogsFromGigaTrees(string infile)
{
    string incontent = Util::GetFileContent(infile);
    Node* gtree = NewickLex::ParseNewickString(incontent);
    TreeIterator* it = gtree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            if (n->GetLabel().find("Ev=1") != std::string::npos)
            {
                n->SetLabel("Dup");
            }
        }
    }
    gtree->CloseIterator(it);

    vector<string> relations = GetOrthologsParalogs(gtree);
    string outstring = "";
    for (int r = 0; r < relations.size(); r++)
    {
        if (r != 0)
            outstring += ";;";
        outstring += relations[r];
    }

    Util::WriteFileContent(Util::ReplaceAll(infile, ".tree", ".relations"), outstring);

}



int main(int argc, char *argv[])
{


    //This reads all the gene trees from a single simphy directory, and
    //creates a file containing the species tree, each gene tree with multiplied
    //dup edges, and the list of relations
    //freopen("W:/ZN/test.txt","w",stdout);


    int verbose = 0;
    map<string, string> args;

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

    //for testing
    //args["m"] = "ortho_from_gigatree";
    //args["f"] = "C:/Users/Manuel/Desktop/giga.tree";


    string mode = "";

    if (args.find("m") != args.end())
    {
        mode = args["m"];
		cout<<"Mode set to "<<mode<<endl;
    }

    //SPECIAL MODE = parse a file from giga, and output its orthologs/paralogs
    if (mode == "ortho_from_gigatree")
    {
        if (args.find("f") != args.end())
        {
            ExtractOrthologsParalogsFromGigaTrees(args["f"]);
        }
        else
        {
            cout<<"You need an input file (use -f)"<<endl;
        }
        return 0;
    }
	
	
	if (mode == "random_nni")
	{
		int nbnni = 1;
		
		if (args.find("n") != args.end())
		{
			nbnni = Util::ToInt(args["n"]);
		}
		
		
		
		if (args.find("f") != args.end())
        {
			string treestr = Util::ReplaceAll(Util::GetFileContent(args["f"]), "\n", "");
            Node* tree = NewickLex::ParseNewickString(treestr, false);
			
			//cout<<"In tree = " <<treestr<<endl;
			for (int i = 0; i < nbnni; i++)
			{
				DoRandomNNI(tree, 0.2);
			}
		
			if (args.find("o") != args.end())
			{
				Util::WriteFileContent(args["o"], Util::ReplaceAll(NewickLex::ToNewickString(tree), " ", ""));
			}
			else
			{
				cout<<Util::ReplaceAll(NewickLex::ToNewickString(tree), " ", "")<<endl;
			}
        }
        else
        {
            cout<<"You need an input file (use -f)"<<endl;
        }
		
		return 0;
	}


	////////////////////////////////////////////////////////////////////////////////////////
	/* Mode restrict_species_tree : just restrict a species tree to a predefined set */
	////////////////////////////////////////////////////////////////////////////////////////
	if (mode == "restrict_species_tree")
	{
		Node* stree = NewickLex::ParseNewickString(args["s"] , false);
		string species_restriction = args["l"];
		
		vector<string> sz = Util::Split(species_restriction, ",");
		vector<int> founds;
		for (int i = 0; i < sz.size(); i++)
		{
			founds.push_back(0);
		}
		
		set<Node*> leavesToKeep;
		TreeIterator* it = stree->GetPostOrderIterator();
		while (Node* n = it->next())
		{
			if (n->IsLeaf())
			{
				for (int i = 0; i < sz.size(); i++)
				{
					if (sz[i] == n->GetLabel())
					{
						leavesToKeep.insert(n);
						founds[i] = founds[i] + 1;
						break;
					}
				}
			}
			else
			{
				n->SetLabel("");
			}
			
		}
		stree->CloseIterator(it);
		
		
		for (int i = 0; i < sz.size(); i++)
		{
			if (founds[i] == 0)
			{
				cout<<"Warning: species "<<sz[i]<<" not found in species tree"<<endl;
			}
		}
		
		 Node::RestrictToLeafset(stree, leavesToKeep);
		 Util::WriteFileContent(args["o"], Util::ReplaceAll(NewickLex::ToNewickString(stree), " ", ""));
		 
		 delete stree;
		 return 0;
		 
	}
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////

    //string baseDir = "C:/Users/Manuel/Dropbox/projects/simgraphs/simphy/DL0000005/01";

    string baseDir = "";

    if (args.find("d") != args.end())
    {
        baseDir = args["d"];
        baseDir = Util::ReplaceAll(baseDir, "\\", "/");
    }

    if (baseDir == "")
    {
        cout<<"Please specify indir with -d"<<endl;
        return 0;
    }

    if (baseDir[baseDir.length() - 1] != '/')
        baseDir = baseDir + "/";





    //DEFAULT MODE = go into a simphy simulation directory,
    //create a mult subdirectory,
    //make a summary file containing each locus tree converted to a gene tree
    //with duplication branch lengths multiplied by 1, 1.5,2, 4,6,8,10,50
    //then make a control file for indelible

    /*
     * We don't care about the locus / gene tree differentiation.
     * However, we need the branch lengths from the gene trees, and the event (Dup/Spec)
     * from the locus tree.  So we merge the two, creating a new tree with
     * the node names and events of the locus tree, and the br len of the gene tree.
     * The node type is embedded in the newick format.
     */

    string outstring = "";

    string speciesTreeStr = Util::ReplaceAll(Util::GetFileContent(baseDir + "s_tree.trees"), "\n", "");

    Node* sp_tree = NewickLex::ParseNewickString(speciesTreeStr, true);

    outstring += "SPECIESTREE=" + speciesTreeStr + "\n";

    vector<string> locusTrees = Util::GetFileLines(baseDir + "l_trees.trees");

    vector<string> geneTrees;
    for (int i = 1; i <= locusTrees.size(); i++)
    {


        string fname = baseDir + "g_trees" +
                        (i < 10 ? "0" : "") +
                       Util::ToString(i) + ".trees";
        geneTrees.push_back( Util::ReplaceAll(Util::GetFileContent(fname), "\n", "") );

    }


    vector<double> multiples;
    multiples.push_back(1.0);
    multiples.push_back(1.5);
    multiples.push_back(2.0);
    multiples.push_back(4.0);
    multiples.push_back(6.0);
    multiples.push_back(8.0);
    multiples.push_back(10.0);
    multiples.push_back(50.0);

    vector< vector<Node*> > multipliedGeneTrees;

    for (int i = 0; i < locusTrees.size(); i++)
    {
        //cout<<"Parsing tree "<<i<<endl;


        int nbdups = 0;

        string padded_i = ((i + 1) < 10 ? "0" : "") + Util::ToString(i + 1);

        outstring += "LOCUS=" + padded_i + "\n";

        //compute node type (Dup or Spec)
        vector<string> locmap_lines = Util::GetFileLines(baseDir + padded_i + ".mapsl");
        map<string, string> loctypes;
        for (int ll = 0; ll < locmap_lines.size(); ll++)
        {
            vector<string> pz = Util::Split( Util::ReplaceAll(locmap_lines[ll], "'", ""), "\t");
            if (pz.size() >= 3)
            {
                loctypes[pz[0]] = pz[2];
                //THIS MIGHT COUNT DUPS THAT APPEAR NOWHERE IN THE GENE TREE (due to losses after dup)
                //if (pz[2] == "Dup")
                //    nbdups++;
            }
        }

        //compute gene to locus map
        vector<string> genemap_lines = Util::GetFileLines(baseDir + padded_i + "l1g.maplg");
        map<string, string> geneloci;
        for (int ll = 0; ll < genemap_lines.size(); ll++)
        {
            vector<string> pz = Util::Split(  Util::ReplaceAll(genemap_lines[ll], "'", ""), "\t");
            if (pz.size() >= 2)
            {
                geneloci[pz[0]] = pz[1];

            }
        }




        Node* gtree = NewickLex::ParseNewickString(geneTrees[i]);

        TreeIterator* it = gtree->GetPostOrderIterator();
        while (Node* n = it->next())
        {

            string loc = geneloci[n->GetLabel()];
            string type = loctypes[loc];

            //cout<<n->GetLabel()<<" "<<loc<<" "<<type<<endl;


            n->SetLabel( Util::ReplaceAll(loc, "_", "__") + "__" + type );

            if (type == "Dup")
                nbdups++;
        }
        gtree->CloseIterator(it);

        //for some reason, root never has an id nor a type.  We find the type the hard way using lca mapping
        {
            string type = "";

            unordered_map<Node*, Node*> mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(gtree, sp_tree, "__", 0);
            if (mapping[gtree] == mapping[gtree->GetChild(0)] || mapping[gtree] == mapping[gtree->GetChild(1)])
            {
                type = "Dup";
                nbdups++;
            }
            else
                type = "Spec";


            gtree->SetLabel(type);
        }



        outstring += "NBDUPS=" + Util::ToString(nbdups) + "\n";

        vector<string> relations = GetOrthologsParalogs(gtree);
        outstring += "RELATIONS=";
        for (int r = 0; r < relations.size(); r++)
        {
            if (r != 0)
                outstring += ";;";
            outstring += relations[r];
        }
        outstring += "\n";


        outstring += "CLUSTERS=" + GetTrueClustersString(gtree) + "\n";

        /*
         * Here we multiply branch lengths of duplications to simulate the DAD model
         */
        vector<Node*> local_multipliedGeneTrees;

        //cout<<"GENETREE1="<<NewickLex::ToNewickString(gtree, true)<<endl;
        //local_multipliedGeneTrees.push_back(gtree);



        for (int m = 0; m < multiples.size(); m++)
        {
            Node* copy = new Node(false);
            copy->CopyFrom(gtree);

            TreeIterator* it = copy->GetPostOrderIterator();
            while (Node* n = it->next())
            {
                if (Util::EndsWith(n->GetLabel(), "Dup"))
                {

                    //TODO: choose randomly
                    int newlocus_index = FindNewLocusChild(copy, n);
                    double len = n->GetChild(newlocus_index)->GetBranchLength();
                    n->GetChild(newlocus_index)->SetBranchLength( len * multiples[m] );
                }

            }
            copy->CloseIterator(it);

            local_multipliedGeneTrees.push_back(copy);

            outstring += "GENETREE" + padded_i + "_" + Util::ToString(multiples[m]) + "="  + NewickLex::ToNewickString(copy, true) + "\n";
        }



        multipliedGeneTrees.push_back(local_multipliedGeneTrees);



        /*delete gtree;
        for (int t = 0; t < multipliedGeneTrees.size(); t++)
        {
            delete multipliedGeneTrees[t];
        }*/
    }


    if (!Util::DirectoryExists(baseDir + "mult"))
        Util::CreateDirectory(baseDir + "mult");

    cout<<"Writing " + baseDir + "mult/summary.txt"<<endl;
    Util::WriteFileContent(baseDir + "mult/summary.txt", outstring);

    string control = EditINDELibleControlFile(baseDir + "control.txt", multipliedGeneTrees, multiples);

    cout<<"Writing " + baseDir + "mult/control.txt"<<endl;
    Util::WriteFileContent(baseDir + "mult/control.txt", control);


    cout<<"Writing " + baseDir + "mult/speciestree_nobranch.nw"<<endl;
    Util::WriteFileContent(baseDir + "mult/speciestree_nobranch.nw", NewickLex::ToNewickString(sp_tree, false, false));

    delete sp_tree;

    return 0;
}
