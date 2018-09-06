#include "autester.h"

AUTester::AUTester()
{

}






string AUTester::GetFastaFetcherCommand(Node* tree, string &outfile, string labelSplitter, int geneIDIndex, int speciesIndex)
{
    //WE ASSUME GENE LABEL FORMAT IS species_geneid

    //the getfasta command invokes a perl script that fetches all nucleotide sequences of the tree genes
    //this command will be dropped if user chose to use peptide sequences instead
    string commandFasta = "getfasta_b.perl -s \"";

    TreeIterator* it = tree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            vector<string> sg = Util::Split(n->GetLabel(), labelSplitter);

            commandFasta += sg[speciesIndex] + ":" + sg[geneIDIndex] + ";";
            //commandFasta += EnsemblUtil::Instance()->GetSpeciesEnsemblName(info->taxa, GLOB_COMPARA_DBNAME) + ":" + info->geneStableID + ";";

        }
    }
    tree->CloseIterator(it);


    commandFasta += "\"";

    if (outfile != "")
        commandFasta += " > " + outfile;

    commandFasta+= "\n";

    //cout<<fromStdString(commandFasta);

    return commandFasta;

}


void AUTester::FetchFasta(Node* genetree, string &outfile, string labelSplitter, int geneIDIndex, int speciesIndex)
{

    //fasta file might exist
    ifstream ifile(outfile);
    if (ifile) {
        cout<<outfile<<" exists";
        ifile.close();
        return;
    }





    string tstr = "";
    string command = GetFastaFetcherCommand(genetree, tstr, labelSplitter, geneIDIndex, speciesIndex);

    QProcess *fastaProcess = new QProcess();

    fastaProcess->setStandardOutputFile(QString::fromStdString(outfile));


    fastaProcess->start(QString::fromStdString(command));
    fastaProcess->waitForFinished(-1);


    cout<<"DONE FASTA-ING WITH CODE "<<fastaProcess->error();

    delete fastaProcess;
}


void AUTester::AlignFile(string &fastaFile, string &outfile, string &repairedOutfile)
{
    QProcess *clustalProcess = new QProcess();


    string command = "clustalw2 -INFILE=" + fastaFile + " -OUTFILE=" + outfile + " -OUTPUT=NEXUS";
    cout<<"ALIGNING ";

    clustalProcess->start(QString::fromStdString(command));
    clustalProcess->waitForFinished(-1);


    cout<<"DONE ALIGNING WITH CODE "<<clustalProcess->error();

    delete clustalProcess;


    //phyml dislikes the clustal nexus format...nexrepair is a script that prepares stuff for phyml
    QProcess *repairProcess = new QProcess();

    string in = outfile;

    command = "nexrepair.py -i \"" + in + "\" -o \"" + repairedOutfile + "\"";
    cout<<"REPAIRING";

    repairProcess->start(QString::fromStdString(command));
    repairProcess->waitForFinished(-1);


    cout<<"DONE REPAIRING WITH CODE "<<repairProcess->error();

    delete repairProcess;

}





void AUTester::ExecutePhyML(string alignedFile, string treesFile)
{
    ifstream alfile(alignedFile);
    if (!alfile) {
        cout<<alignedFile<<" DOESN'T exists FOR phyml !";
        return;
    }
    alfile.close();

    string command = "phyml -i " + alignedFile +
                     " -u " + treesFile + " -n 1 -o lr --print_site_lnl --no_memory_check";

    QProcess *phymlProcess = new QProcess();


    cout<<"PHYML-ing "<<command;


    phymlProcess->setProcessChannelMode(QProcess::MergedChannels);

    phymlProcess->start(QString::fromStdString(command));


    phymlProcess->waitForReadyRead();

    //in case of error, phyml blocks up and asks for input.  Here we give phyml his damn input
    //so it won't clog up the whole toolchain
    while (phymlProcess->canReadLine())
    {
        QByteArray bline = phymlProcess->readLine();
        QString line(bline);

        if (line.toUpper().contains("TYPE ANY KEY"))
        {
            phymlProcess->putChar('13');
        }
    }
    phymlProcess->waitForFinished(-1);


    cout<<"DONE PHYML-ing WITH CODE "<<phymlProcess->error();

    delete phymlProcess;
}


/**
  CONSEL outputs one line of stats per tree
  1st tree line is the original tree, 2nd is the corrected tree
  **/
void AUTester::ExecuteConsel(string &phyMLFile, string &resultsFile)
{
    string command = "/u/lafonman/consel/consel/bin/makermt --phyml " + phyMLFile;

    QProcess *mkmtProcess = new QProcess();


    cout<<"MAKERMT";

    mkmtProcess->start(QString::fromStdString(command));
    mkmtProcess->waitForFinished(-1);


    cout<<"DONE MAKERMT WITH CODE "<<mkmtProcess->error();

    delete mkmtProcess;




    string noext = Util::ReplaceAll(phyMLFile, ".txt", "");
    command = "/u/lafonman/consel/consel/bin/consel " + noext;

    QProcess *conselProcess = new QProcess();


    cout<<"CONSEL";

    conselProcess->start(QString::fromStdString(command));
    conselProcess->waitForFinished(-1);


    cout<<"DONE MAKERMT WITH CODE "<<conselProcess->error();

    delete conselProcess;



    command = "/u/lafonman/consel/consel/bin/catpv " + noext + ".pv -s 1";

    QProcess *catpvProcess = new QProcess();


    cout<<"CATPV";

    catpvProcess->setStandardOutputFile(QString::fromStdString(resultsFile));
    catpvProcess->start(QString::fromStdString(command));
    catpvProcess->waitForFinished(-1);



    cout<<"DONE CATPV WITH CODE "<<catpvProcess->error();

    delete catpvProcess;

}






//outputs before/after trees for phyml/consel testing
void AUTester::OutputTrees(Node* originalTree, Node* correctedTree, string &outfile, string &treeID, string labelSplitter, int geneIDIndex, int speciesIndex)
{
    //WE ASSUMED GENE LABEL FORMAT IS species_geneid

    //make a copy since we change labels
    Node* treeBefore = new Node(false);
    treeBefore->CopyFrom(originalTree);

    Node* treeAfter = new Node(false);
    treeAfter->CopyFrom(correctedTree);

    TreeIterator* it = treeBefore->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            vector<string> sg = Util::Split(n->GetLabel(), labelSplitter);
            n->SetLabel(sg[geneIDIndex]);
        }
        else
        {
            n->SetLabel("");
        }
    }
    treeBefore->CloseIterator(it);

    it = treeAfter->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            vector<string> sg = Util::Split(n->GetLabel(), labelSplitter);
            n->SetLabel(sg[geneIDIndex]);
        }
        else
        {
            n->SetLabel("");
        }
    }
    treeBefore->CloseIterator(it);


    ofstream myfile;
    myfile.open (outfile);
    myfile << NewickLex::ToNewickString(treeBefore) << "\n" << NewickLex::ToNewickString(treeAfter);
    myfile.close();



    delete treeBefore;
    delete treeAfter;

}

void AUTester::RunAUTest(Node* beforeTree, Node* afterTree, string treeID, string labelSplitter, int geneIDIndex, int speciesIndex,
               string bothTreesDir, string fastaDir, string alignedDir, string conselOutDir)
{

    string conselOut = conselOutDir + treeID + ".consel";
    ifstream resultsfile(conselOut);
    if (resultsfile) {
        cout<<"RESULTS " + conselOut + " exists";
        resultsfile.close();
    }
    else
    {
        string treesFile = bothTreesDir + treeID + ".trees";

        cout<<"Outputting "<<treesFile<<".tree";

        OutputTrees(beforeTree, afterTree, treesFile, treeID, labelSplitter, geneIDIndex, speciesIndex);

        string fastaFile = fastaDir + treeID + ".fasta";
        FetchFasta(beforeTree, fastaFile, labelSplitter, geneIDIndex, speciesIndex);


        string aligned = alignedDir + treeID + ".nxs";
        string repaired = alignedDir + treeID + "-rep.nxs";

        AlignFile(fastaFile, aligned, repaired);


        ExecutePhyML(repaired, treesFile);

        string conselIn = repaired + "_phyml_lk.txt";

        ExecuteConsel(conselIn, conselOut);
    }

}
