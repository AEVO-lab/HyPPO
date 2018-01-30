#ifndef AUTESTER_H
#define AUTESTER_H

#include <string>
#include <map>
#include <unordered_map>
#include "trees/node.h"
#include "trees/newicklex.h"

#include <QProcess>
#include <QString>

#include <iostream>
#include <fstream>

using namespace std;

//TODO : TODO : TODO !!!!!!!!!!!!!!!!!!!!!!!! QT DEPENDENCY HERE !

class AUTester
{
public:
    AUTester();

    string GetFastaFetcherCommand(Node* tree, string &outfile, string labelSplitter, int geneIDIndex, int speciesIndex);
    void FetchFasta(Node* genetree, string &outfile, string labelSplitter, int geneIDIndex, int speciesIndex);
    void AlignFile(string &fastaFile, string &outfile, string &repairedOutfile);
    void ExecutePhyML(string alignedFile, string treesFile);

    /**
      CONSEL outputs one line of stats per tree
      1st tree line is the original tree, 2nd is the corrected tree
      **/
    void ExecuteConsel(string &phyMLFile, string &resultsFile);
    void OutputTrees(Node* originalTree, Node* correctedTree, string &outfile, string &treeID, string labelSplitter, int geneIDIndex, int speciesIndex);
    void RunAUTest(Node* beforeTree, Node* afterTree, string treeID, string labelSplitter, int geneIDIndex, int speciesIndex,
                   string bothTreesDir, string fastaDir, string alignedDir, string conselOutDir);

};


#endif // AUTESTER_H
