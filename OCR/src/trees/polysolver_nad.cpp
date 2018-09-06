#include "polysolver_nad.h"


SADNADGraph::SADNADGraph()
{

}

void SADNADGraph::AddNode(Node* n)
{
    SADNADNode* sad = new SADNADNode();
    this->nodes[n] = sad;

    //n is in its own cc until proven otherwise
    set<Node*> cc;
    cc.insert(n);
    //this->ADComponents[n] = cc;
    sad->ADComponent = n;
}

void SADNADGraph::AddEdge(Node* n1, Node* n2, string edgeType)
{
    SADNADNode* s1 = nodes[n1];
    SADNADNode* s2 = nodes[n2];
    if (edgeType == "AD")
    {
        s1->AD_Neighbors.insert(n2);
        s2->AD_Neighbors.insert(n1);

        //merge the 2 cc's if we need to
        if (s1->ADComponent != s2->ADComponent)
        {
            /*set<Node*> cc1 = GetADComponentOf(n1);
            set<Node*> cc2 = GetADComponentOf(n2);

            for (set<Node*>::iterator it = cc2.begin(); it != cc2.end(); it++)
            {
                cc1.insert((*it));
                //nodes[(*it)]->ADComponent = s1->ADComponent;
            }
            cc2.clear();
            this->ADComponents.erase(s2->ADComponent);
            s2->ADComponent = s1->ADComponent;*/

            //AD component of s1 globs the others...yes, globs
            for (map<Node*, SADNADNode*>::iterator it = nodes.begin(); it != nodes.end(); it++)
            {
                Node* n = (*it).first;
                SADNADNode* sad = (*it).second;

                if (n != n2 && sad->ADComponent == s2->ADComponent)
                {
                    sad->ADComponent = s1->ADComponent;
                }
            }
            s2->ADComponent = s1->ADComponent;
        }
    }
    else if (edgeType == "S")
    {
        s1->S_Neighbors.insert(n2);
        s2->S_Neighbors.insert(n1);
    }
    else if (edgeType == "NAD")
    {
        s1->NAD_Neighbors.insert(n2);
        s2->NAD_Neighbors.insert(n1);
    }
}

/*set<Node*> SADNADGraph::GetADComponentOf(Node* n)
{
    SADNADNode* s1 = nodes[n];
    return this->ADComponents[s1->ADComponent];
}*/

Node* SADNADGraph::GetADComponentRepresentantOf(Node* n)
{
    SADNADNode* s1 = nodes[n];
    return s1->ADComponent;
}

string SADNADGraph::GetEdgeType(Node* n1, Node* n2, unordered_map<Node*, Node*> &lcaMapping)
{
    Node* s1 = lcaMapping[n1];
    Node* s2 = lcaMapping[n2];
    if (!s1->HasAncestor(s2) && !s2->HasAncestor(s1))
    {
        return "S";
    }
    else
    {
        bool hasOne = false;
        unordered_set<Node*> n1Species = GeneSpeciesTreeUtil::Instance()->GetGeneTreeSpecies(n1, lcaMapping);
        TreeIterator* it = n2->GetPostOrderIterator(true);
        while (Node* n2leaf = it->next())
        {
            if (n1Species.find(lcaMapping[n2leaf]) != n1Species.end())
            {
                hasOne = true;
                break;
            }
        }
        n2->CloseIterator(it);

        if (hasOne)
            return "AD";
        else
            return "NAD";

    }
}



vector<pair<Node *, Node *> > SADNADGraph::GetUsefulSpeciationEdges()
{
    set<pair<Node *, Node *> > retset;   //what a waste...just to remove symmetry
    vector<pair<Node *, Node *> > ret;

    for (map<Node*, SADNADNode*>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        SADNADNode* sad = (*it).second;
        Node* n1 = (*it).first;

        for (set<Node*>::iterator sadit = sad->S_Neighbors.begin(); sadit != sad->S_Neighbors.end(); sadit++)
        {
            Node* n2 = (*sadit);

            if (!this->HaveSameADComponent(n1, n2))
            {
                pair<Node*, Node*> symedge = make_pair(n2, n1);

                if (retset.find(symedge) == retset.end())
                {
                    ret.push_back( make_pair(n1, n2) );
                    retset.insert( make_pair(n1, n2) );
                }
            }
        }
    }

    return ret;
}


int SADNADGraph::GetNbADComponents()
{
    set<Node*> adreps;
    for (map<Node*, SADNADNode*>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        Node* n = (*it).first;
        Node* rep = this->GetADComponentRepresentantOf(n);

        if (adreps.find(rep) == adreps.end())
            adreps.insert(rep);
    }

    return adreps.size();
}

int SADNADGraph::GetNbADComponents(set<pair<Node*, Node*> > &additionalSEdges)
{
    //this fakes another sadnadgraph where each AD edge and each additionalEdge is inserted as an AD edge
    //Probably suboptimal but easy to implement
    SADNADGraph graph;

    set<pair<Node*, Node*> > uniqueADEdges;

    for (map<Node*, SADNADNode*>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        Node* n = (*it).first;
        SADNADNode* sad = (*it).second;
        graph.AddNode(n);

        for (set<Node*>::iterator adit = sad->AD_Neighbors.begin(); adit != sad->AD_Neighbors.end(); adit++)
        {
            pair<Node*, Node*> erev = make_pair((*adit), n);
            if (uniqueADEdges.find(erev) == uniqueADEdges.end())
            {
                uniqueADEdges.insert(make_pair(n, (*adit)));
            }
        }
    }

    for (set<pair<Node*, Node*> >::iterator it = uniqueADEdges.begin(); it != uniqueADEdges.end(); it++)
    {
        pair<Node*, Node*> e = (*it);
        graph.AddEdge(e.first, e.second, "AD");
    }

    for (set<pair<Node*, Node*> >::iterator it = additionalSEdges.begin(); it != additionalSEdges.end(); it++)
    {
        pair<Node*, Node*> e = (*it);
        graph.AddEdge(e.first, e.second, "AD");
    }

    return graph.GetNbADComponents();

}



bool SADNADGraph::HasSEdge(Node* n1, Node* n2)
{
    SADNADNode* sad1 = this->nodes[n1];

    return (sad1->S_Neighbors.find(n2) != sad1->S_Neighbors.end());
}

void SADNADGraph::BuildGraph(vector<Node*> nodes, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping)
{
    //TODO : this is suboptimal in terms of complexity.  Each GetEdgeType takes something like linear time,
    //       summing up to n^3 complexity.  But, can we do better anyway ?
    for (int i = 0; i < nodes.size(); i++)
    {
        this->AddNode(nodes[i]);

        for (int j = 0; j < i; j++)
        {
            string edgeType = this->GetEdgeType(nodes[i], nodes[j], lcaMapping);
            this->AddEdge(nodes[i], nodes[j], edgeType);
        }
    }
}



void SADNADGraph::MergeNodes(Node* n1, Node* n2, Node* newNode, unordered_map<Node*, Node*> &lcaMapping)
{
    delete this->nodes[n1];
    delete this->nodes[n2];
    this->nodes.erase(n1);
    this->nodes.erase(n2);
    this->AddNode(newNode);

    for (map<Node*, SADNADNode*>::iterator it = this->nodes.begin(); it != this->nodes.end(); it++)
    {
        Node* n = (*it).first;
        SADNADNode* sad = (*it).second;

        if (n != newNode)
        {
            string type1 = sad->GetEdgeTypeWith(n1);
            string type2 = sad->GetEdgeTypeWith(n2);

            //infer the new edge relationship between n and newNode
            if (type1 == "AD" || type2 == "AD")
            {
                this->AddEdge(n, newNode, "AD");
            }
            else if (type1 == "NAD" || type2 == "NAD")
            {
                this->AddEdge(n, newNode, "NAD");
            }
            else
            {
                if (lcaMapping[n]->HasAncestor(lcaMapping[newNode]))
                    this->AddEdge(n, newNode, "NAD");
                else
                    this->AddEdge(n, newNode, "S");
            }


            sad->RemoveNeighbor(n1);
            sad->RemoveNeighbor(n2);
        }
    }
}


bool SADNADGraph::HaveSameADComponent(Node* n1, Node* n2)
{
    SADNADNode* s1 = nodes[n1];
    SADNADNode* s2 = nodes[n2];

    return s1->ADComponent == s2->ADComponent;
}


SADNADGraph::~SADNADGraph()
{
    map<Node*, SADNADNode*>::iterator it = nodes.begin();
    while (it != nodes.end())
    {
        delete (*it).second;
        it++;
    }
    nodes.clear();
}



void SADNADGraph::PrintGraph()
{
    string str = "";
    for (map<Node*, SADNADNode*>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        Node* n = (*it).first;
        SADNADNode* sad = (*it).second;

        str += "NODE : " + n->GetLabel();

        str += "\nS : ";
        for (set<Node*>::iterator it2 = sad->S_Neighbors.begin(); it2 != sad->S_Neighbors.end(); it2++)
        {
            Node* n2 = (*it2);
            str += " " + n2->GetLabel();
        }

        str += "\nNAD : ";
        for (set<Node*>::iterator it2 = sad->NAD_Neighbors.begin(); it2 != sad->NAD_Neighbors.end(); it2++)
        {
            Node* n2 = (*it2);
            str += " " + n2->GetLabel();
        }

        str += "\nAD : ";
        for (set<Node*>::iterator it2 = sad->AD_Neighbors.begin(); it2 != sad->AD_Neighbors.end(); it2++)
        {
            Node* n2 = (*it2);
            str += " " + n2->GetLabel();
        }

        str += "\n\n";
    }

    cout<<str;
}


PolySolverNAD::PolySolverNAD()
{
}

PolySolverCorrectionInfo PolySolverNAD::CorrectHighestNAD(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping)
{



    unordered_map<Node*, Node*> oldlcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, geneLeavesSpeciesMapping);
    unordered_map<Node*, Node*> lcaMapping;

    Node* geneTreeCopy = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(geneTree, oldlcaMapping, lcaMapping);

    //GeneSpeciesTreeUtil::Instance()->PrintMapping(geneTreeCopy, lcaMapping);

    TreeIterator* it = geneTreeCopy->GetPreOrderIterator();
    while (Node* n = it->next())
    {
        if (!n->IsLeaf())
        {
            //first check if it's a duplication, lca mapping classic rule
            if (lcaMapping[n->GetChild(0)] == lcaMapping[n] || lcaMapping[n->GetChild(1)] == lcaMapping[n])
            {
                if (!GeneSpeciesTreeUtil::Instance()->HaveCommonSpecies(n->GetChild(0), n->GetChild(1), lcaMapping))
                {
                    vector<string> leafLabels;
                    vector<Node*> n_leaves = n->GetLeafVector();
                    for (int i = 0; i < n_leaves.size(); i++)
                    {
                        leafLabels.push_back(n_leaves[i]->GetLabel());
                    }

                    //there it is ! the highest NAD
                    pair<Node*, unordered_map<Node*, Node*> > polytomizedWithMapping = PolytomizeNAD(n, speciesTree, lcaMapping);
                    Node* polytomized = polytomizedWithMapping.first;


                    //HERE we do some not so clean stuff...because whatever we do, we'll exit this function
                    geneTreeCopy->CloseIterator(it);

                    //replace the subtree that just got polytomized
                    if (!n->IsRoot())
                    {
                        Node* parent = n->GetParent();
                        parent->RemoveChild(n);
                        parent->AddSubTree(polytomized);
                        delete n;
                    }
                    else
                    {
                        delete geneTreeCopy;
                        geneTreeCopy = polytomized;
                    }

                    //cout<<"COPY AFTER = "<<NewickLex::ToNewickString(geneTreeCopy)<<endl;


                    PolySolverCorrectionInfo info;
                    info.nadCladeGenes = leafLabels;
                    info.firstPolySize = polytomized->GetChild(0)->GetNbChildren();
                    info.secondPolySize = polytomized->GetChild(1)->GetNbChildren();
                    this->SolvePolytomy(polytomized->GetChild(0), speciesTree, polytomizedWithMapping.second);
                    this->SolvePolytomy(polytomized->GetChild(1), speciesTree, polytomizedWithMapping.second);

                    //cout<<"CORRECTED = "<<NewickLex::ToNewickString(geneTreeCopy)<<endl;

                    geneTreeCopy->DeleteSingleChildDescendants();

                    info.correction = geneTreeCopy;

                    return info;
                }
            }
        }
    }
    geneTreeCopy->CloseIterator(it);

    //if we got here, we found no NAD, and since we didn't do anything we return NULL
    delete geneTreeCopy;

    PolySolverCorrectionInfo info;
    info.correction = NULL;
    return info;
}



PolySolverCorrectionInfo PolySolverNAD::CorrectNodeByMultifurcation(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping, Node* n)
{
    //TODO : code copied from above

    unordered_map<Node*, Node*> oldlcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, geneLeavesSpeciesMapping);
    unordered_map<Node*, Node*> lcaMapping;

    //here we'll copy the original gene tree and manage to find the node of interest in this copy
    string prevLabel = n->GetLabel();
    string tempLabel = "temp-label-no-one-else-should-use";
    n->SetLabel(tempLabel);
    Node* geneTreeCopy = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(geneTree, oldlcaMapping, lcaMapping);
    n->SetLabel(prevLabel);

    //find the node of interest
    Node* node_to_correct = NULL;
    TreeIterator* it = geneTreeCopy->GetPostOrderIterator();
    while (Node* ncopy = it->next())
    {
        if (ncopy->GetLabel() == tempLabel)
        {
            node_to_correct = ncopy;
            node_to_correct->SetLabel(prevLabel);
            break;
        }

    }
    geneTreeCopy->CloseIterator(it);

    vector<string> leafLabels;
    vector<Node*> n_leaves = node_to_correct->GetLeafVector();
    for (int i = 0; i < n_leaves.size(); i++)
    {
        leafLabels.push_back(n_leaves[i]->GetLabel());
    }

    pair<Node*, unordered_map<Node*, Node*> > polytomizedWithMapping = PolytomizeNAD(node_to_correct, speciesTree, lcaMapping);
    Node* polytomized = polytomizedWithMapping.first;


    //replace the subtree that just got polytomized
    if (!node_to_correct->IsRoot())
    {
        Node* parent = node_to_correct->GetParent();
        parent->RemoveChild(node_to_correct);
        parent->AddSubTree(polytomized);
        delete node_to_correct;
    }
    else
    {
        delete geneTreeCopy;
        geneTreeCopy = polytomized;
    }


    PolySolverCorrectionInfo info;
    info.nadCladeGenes = leafLabels;
    info.firstPolySize = polytomized->GetChild(0)->GetNbChildren();
    info.secondPolySize = polytomized->GetChild(1)->GetNbChildren();
    this->SolvePolytomy(polytomized->GetChild(0), speciesTree, polytomizedWithMapping.second);
    this->SolvePolytomy(polytomized->GetChild(1), speciesTree, polytomizedWithMapping.second);


    geneTreeCopy->DeleteSingleChildDescendants();

    info.correction = geneTreeCopy;

    return info;
}


pair<Node*, unordered_map<Node*, Node*> > PolySolverNAD::PolytomizeNAD(Node* nadNode, Node* speciesTree, unordered_map<Node*, Node*> lcaMapping)
{
    set<Node*> leftSubtrees, rightSubtrees;
    Node* s = lcaMapping[nadNode];

    //TODO : there should be a way not to iterate uselessly into a taken subtree (preorder traversal that we stop)
    TreeIterator* it = nadNode->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n != nadNode)
        {
            //here we maximal subtree either on the left or right
            if (lcaMapping[n] != s && lcaMapping[n->GetParent()] == s)
            {
                if (lcaMapping[n]->HasAncestor(s->GetChild(0)))
                {
                    leftSubtrees.insert(n);
                }
                else    //if (lcaMapping[n]->HasAncestor(s->GetChild(1))) should be the only possibility here
                {
                    rightSubtrees.insert(n);
                }
            }
        }
    }
    nadNode->CloseIterator(it);


    Node* newShizzle = new Node(false);
    Node* left = newShizzle->AddChild();
    Node* right = newShizzle->AddChild();
    unordered_map<Node*, Node*> newMapping;

    for (set<Node*>::iterator itLeft = leftSubtrees.begin(); itLeft != leftSubtrees.end(); itLeft++)
    {
        Node* copy = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping((*itLeft), lcaMapping, newMapping);
        left->AddSubTree(copy);
    }

    for (set<Node*>::iterator itRight = rightSubtrees.begin(); itRight != rightSubtrees.end(); itRight++)
    {
        Node* copy = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping((*itRight), lcaMapping, newMapping);
        right->AddSubTree(copy);
    }

    newMapping[newShizzle] = s;
    if (left->GetNbChildren() > 1)
    {
        newMapping[left] = GeneSpeciesTreeUtil::Instance()->GetSingleNodeLCAMapping(left, speciesTree, newMapping);
    }
    if (right->GetNbChildren() > 1)
    {
        newMapping[right] = GeneSpeciesTreeUtil::Instance()->GetSingleNodeLCAMapping(right, speciesTree, newMapping);
    }

    newShizzle->DeleteSingleChildDescendants();

    return make_pair(newShizzle, newMapping);
}


Node* PolySolverNAD::SolvePolytomies(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping)
{
    //TODO : THIS HASN'T BEEN TESTED AFTER SOME MODIFICATIONS !!
    Node* geneTreeCopy;

    unordered_map<Node*, Node*> mappingCopy;
    geneTreeCopy = GeneSpeciesTreeUtil::Instance()->CopyTreeWithNodeMapping(geneTree, geneLeavesSpeciesMapping, mappingCopy);

    unordered_map<Node*, Node*> lcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTreeCopy, speciesTree, mappingCopy);

    TreeIterator* it = geneTreeCopy->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        vector<Node*> curLeaves;
        for (int i = 0; i < g->GetNbChildren(); i++)
        {
            curLeaves.push_back(g->GetChild(i));
        }
        this->SolvePolytomy(curLeaves, speciesTree, lcaMapping);
    }

    return geneTreeCopy;
}

Node* PolySolverNAD::SolvePolytomy(Node* polytomyRoot, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping)
{
    vector<Node*> curLeaves;
    for (int i = 0; i < polytomyRoot->GetNbChildren(); i++)
    {
        curLeaves.push_back(polytomyRoot->GetChild(i));
    }

    if (curLeaves.size() == 0)
        return polytomyRoot;
    return this->SolvePolytomy(curLeaves, speciesTree, lcaMapping);
}


Node* PolySolverNAD::SolvePolytomy(vector<Node*> leaves, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping)
{


    SADNADGraph graph;
    graph.BuildGraph(leaves, speciesTree, lcaMapping);



    //successively apply lowest useful speciations
    TreeIterator* sit = speciesTree->GetPostOrderIterator();
    while (Node* s = sit->next())
    {
        //gLeft : the nodes of leaves with a species on the left of s
        //gRight : the nodes of leaves with a species on the right of s
        //gTaken : the nodes already used in a speciation
        set<Node*> gLeftGuys, gRightGuys, gTaken;
        if (!s->IsLeaf())
        {
            Node* sLeft = s->GetChild(0);
            Node* sRight = s->GetChild(1);

            //This could be optimized, but code is simpler and more easily modifiable by just beginning by building sLeft and sRight
            const map<Node*, SADNADNode*> graphNodes = graph.GetNodes();
            for (map<Node*, SADNADNode*>::const_iterator graphIt = graphNodes.begin(); graphIt != graphNodes.end(); graphIt++)
            {
                Node* g = (*graphIt).first;
                if (!lcaMapping[g])
                {
                    cout<<"CRITICAL PROBLEM "<<g->GetLabel();
                    cout<<endl;
                }

                if (lcaMapping[g] == sLeft || lcaMapping[g]->HasAncestor(sLeft))
                {
                    gLeftGuys.insert(g);
                }
                else if (lcaMapping[g] == sRight || lcaMapping[g]->HasAncestor(sRight))
                {
                    gRightGuys.insert(g);
                }

            }


            //make speciations GREEDILY
            /*for (set<Node*>::iterator itLeft = gLeftGuys.begin(); itLeft != gLeftGuys.end(); itLeft++)
            {
                bool wasJoined = false;
                Node* gLeft = (*itLeft);
                for (set<Node*>::iterator itRight = gRightGuys.begin(); itRight != gRightGuys.end(); itRight++)
                {
                    if (!wasJoined)
                    {
                        Node* gRight = (*itRight);

                        //useful spec x y : x y not in same component, y not already used
                        if ( gTaken.find(gRight) == gTaken.end() && !graph.HaveSameADComponent(gLeft, gRight))
                        {
                            Node* newNode = gLeft->InsertParentWith(gRight);
                            lcaMapping[newNode] = s;
                            graph.MergeNodes(gLeft, gRight, newNode, lcaMapping);
                            gTaken.insert(gRight);
                            wasJoined = true;
                        }
                    }
                }
            }*/

            this->SpeciateCleverly(graph, gLeftGuys, gRightGuys, s, lcaMapping);
            gLeftGuys.clear();
            gRightGuys.clear();
        }
    }
    speciesTree->CloseIterator(sit);


    last_nb_ad_components = graph.GetNbADComponents();

    //OK then...from here, we've made all the useful speciations we could.
    //Next, we merge AD nodes
    //TODO : here, we do this the ugly way
    bool thereIsAnAD = true;
    while (thereIsAnAD)
    {
        thereIsAnAD = false;

        Node* chosenNode = NULL;
        Node* chosenFriend = NULL;
        const map<Node*, SADNADNode*> graphNodes = graph.GetNodes();
        for (map<Node*, SADNADNode*>::const_iterator it = graphNodes.begin(); it != graphNodes.end(); it++)
        {
            Node* n = (*it).first;
            SADNADNode* sad = (*it).second;

            if (sad->AD_Neighbors.size() > 0)
            {
                chosenFriend = (*sad->AD_Neighbors.begin());
                chosenNode = n;
                thereIsAnAD = true;
            }
        }

        if (thereIsAnAD)
        {
            Node* newNode = chosenNode->InsertParentWith(chosenFriend);
            lcaMapping[newNode] = lcaMapping[chosenNode]->FindLCAWith(lcaMapping[chosenFriend]);
            graph.MergeNodes(chosenNode, chosenFriend, newNode, lcaMapping);
        }

    }


    //TODO : the NAD
    //And then, if we get here, we don't have a choice but to create NADs
    while (graph.GetNodes().size() > 1)
    {
        const map<Node*, SADNADNode*> graphNodes = graph.GetNodes();
        map<Node*, SADNADNode*>::const_iterator it = graphNodes.begin();

        Node* n1 = (*it).first;
        it++;
        Node* n2 = (*it).first;
        Node* newNode = n1->InsertParentWith(n2);
        lcaMapping[newNode] = lcaMapping[n1]->FindLCAWith(lcaMapping[n2]);
        graph.MergeNodes(n1, n2, newNode, lcaMapping);
    }




    return (*graph.GetNodes().begin()).first;
}






void PolySolverNAD::SpeciateCleverly(SADNADGraph &graph, set<Node*> &gLeftGuys, set<Node*> &gRightGuys, Node *curSpecies, unordered_map<Node *, Node *> &lcaMapping)
{

    //cout<<"STEP !   LEFT="<<gLeftGuys.size()<<" RIGHT="<<gRightGuys.size()<<endl;
    vector<vector<Node*> > compsLeft = GetSortedLocalADComponents(graph, gLeftGuys, gRightGuys);
    vector<vector<Node*> > compsRight = GetSortedLocalADComponents(graph, gRightGuys, gLeftGuys);

    //if (gLeftGuys.size() > 0)
    //    cout<<"LEFTSIZE="<<gLeftGuys.size()<<" COMPSLEFTSIZE="<<compsLeft.size();
    vector<Node*> leftOrdering;
    for (int i = 0; i < compsLeft.size(); i++)
    {
        for (int j = 0; j < compsLeft[i].size(); j++)
        {
            leftOrdering.push_back(compsLeft[i][j]);
            //cout<<":"<<compsLeft[i][j]->GetLabel() + " ";
        }
    }

    //if (gLeftGuys.size() > 0)
    //    cout<<endl;
    //if (gRightGuys.size() > 0)
    //    cout<<"RIGHTSIZE="<<gRightGuys.size()<<" ";
    vector<Node*> rightOrdering;
    for (int i = 0; i < compsRight.size(); i++)
    {
        for (int j = 0; j < compsRight[i].size(); j++)
        {
            rightOrdering.push_back(compsRight[i][j]);
            //cout<<":"<<compsRight[i][j]->GetLabel() + " ";
        }
    }

    //if (gRightGuys.size() > 0)
    //    cout<<endl;


    set<Node*> gTaken;

    for (int l = 0; l < leftOrdering.size(); l++)
    {
        bool wasJoined = false;
        Node* gLeft = leftOrdering[l];
        for (int r = 0; r < rightOrdering.size() && !wasJoined; r++)
        {
            Node* gRight = rightOrdering[r];
            if (gTaken.find(gLeft) == gTaken.end() && gTaken.find(gRight) == gTaken.end())
            {
                if (!graph.HaveSameADComponent(gLeft, gRight))
                {
                    //cout<<"JOINING "<<gLeft->GetLabel()<<" "<<gRight->GetLabel()<<endl;
                    Node* newNode = gLeft->InsertParentWith(gRight);
                    newNode->SetLabel(gLeft->GetLabel() + "__" + gRight->GetLabel());
                    lcaMapping[newNode] = curSpecies;
                    graph.MergeNodes(gLeft, gRight, newNode, lcaMapping);
                    gTaken.insert(gRight);
                    gTaken.insert(gLeft);
                    wasJoined = true;
                }
            }
        }
    }

}

struct vector_size_sorter
{
    inline bool operator() (const vector<Node*>& struct1, const vector<Node*>& struct2)
    {
        return (struct1.size() > struct2.size());
    }
};

vector<vector<Node*> > PolySolverNAD::GetSortedLocalADComponents(SADNADGraph &graph, set<Node*> guys, set<Node*> otherGuys)
{
    //TODO : I was sick when I wrote the fake node addition thing.  It could probably be a lot better...
    map<Node*, vector<Node*> > localADComps;


    for (set<Node*>::iterator leftIt = guys.begin(); leftIt != guys.end(); leftIt++)
    {
        Node* g = (*leftIt);

        vector<Node*> comp;
        Node* gCompRep = graph.GetADComponentRepresentantOf(g);
        if (localADComps.find(gCompRep) != localADComps.end())
            comp = localADComps[gCompRep];
        comp.push_back(g);

        localADComps[gCompRep] = comp;
    }


    //for those AD-comps with a bridge, we need to add one fake node
    //that might break out ties in the sort and it is important
    for (map<Node*, vector<Node*> >::iterator it = localADComps.begin(); it != localADComps.end(); it++)
    {
        vector<Node*> curComp = (*it).second;
        bool hasBridge = false;

        for (int i = 0; i < curComp.size() && !hasBridge; i++)
        {
            Node* n1 = curComp[i];
            for (set<Node*>::iterator otherIt = otherGuys.begin(); otherIt != otherGuys.end() && !hasBridge; otherIt++)
            {
                Node* n2 = (*otherIt);

                if (graph.HaveSameADComponent(n1, n2))
                    hasBridge = true;
            }
        }

        if (hasBridge)
        {
            Node* nfake = new Node(false);
            //cout<<"BRIDGE ADDED TO CC OF "<<curComp[0]->GetLabel()<<endl;
            nfake->SetLabel("fake");
            curComp.push_back(nfake);
            localADComps[(*it).first] = curComp;
        }
    }


    //sort the local comps by size
    vector<vector<Node*> > v_localADComps;
    for (map<Node*, vector<Node*> >::iterator lit = localADComps.begin(); lit != localADComps.end(); lit++)
    {
        v_localADComps.push_back((*lit).second);
        //cout<<"COMP OF "<<(*lit).second[0]->GetLabel()<<"="<<(*lit).second.size()<<endl;
    }

    std::sort(v_localADComps.begin(), v_localADComps.end(), vector_size_sorter());




    vector<vector<Node*> > v_final_ad_comps;
    //then, delete any fakeass node we added
    for (int i = 0; i < v_localADComps.size(); i++)
    {
        vector<Node*> curComp = v_localADComps[i];

        //cout<<i<<" : "<<curComp[0]->GetLabel()<<" SIZE="<<curComp.size()<<endl;

        Node* last = curComp[curComp.size() - 1];
        if (last->GetLabel() == "fake")
        {
            curComp.pop_back();
            delete last;

            //localADComps[(*it).first] = curComp;
        }
        v_final_ad_comps.push_back(curComp);

    }


    return v_final_ad_comps;
}



pair<Node*, Node*> PolySolverNAD::GetRandomPolytomy(int k, int verbose)
{
    Node* speciesTree = new Node(false);

    double s_size_factor = 2.5 * (double)(rand() % 1000)/1000.0 + 0.5;  //between 0.5 and 3

    for (int i = 0; i < s_size_factor*k; i++)
    {
        Node* c = speciesTree->AddChild();
        c->SetLabel("S" + Util::ToString(i));
    }

    speciesTree->BinarizeRandomly();

    //get an ordering of the internal nodes...this will let us pick one at random
    vector<Node*> internalNodes;

    TreeIterator* it = speciesTree->GetPostOrderIterator(false);
    while (Node* s = it->next())
    {
        if (!s->IsLeaf())
        {
            internalNodes.push_back(s);
        }

    }
    speciesTree->CloseIterator(it);

    //generate k gene subtrees
    unordered_map<Node*, Node*> lcaMapping;
    vector<Node*> forest;
    map<Node*, Node*> geneLeftSpecies;
    map<Node*, Node*> geneRightSpecies;

    for (int i = 0; i < k; i++)
    {
        Node* g = new Node(false);
        g->SetLabel("G" + Util::ToString(i));

        //pick an lca for g at random
        Node* lca = internalNodes[rand() % internalNodes.size()];
        lca->SetLabel(lca->GetLabel() + "_" + Util::ToString(i));
        lcaMapping[g] = lca;

        //add something left and right to enforce s(g) = lca
        //by adding a species specific to g on both sides
        bool done = false;
        TreeIterator* itLeft = lca->GetChild(0)->GetPostOrderIterator();
        while (Node* s = itLeft->next())
        {
            if (!done)
            {
                string slbl = s->GetLabel();
                if (slbl[0] == 'S') //got an original species leaf
                {
                    Node* sg = s->AddChild();
                    sg->SetLabel("XL" + Util::ToString(i));

                    Node* gs = g->AddChild();
                    gs->SetLabel("XL" + Util::ToString(i));

                    lcaMapping[gs] = sg;
                    done = true;

                    geneLeftSpecies[g] = s;
                }
            }
        }
        lca->CloseIterator(itLeft);


        done = false;
        TreeIterator* itRight = lca->GetChild(1)->GetPostOrderIterator();
        while (Node* s = itRight->next())
        {
            if (!done)
            {
                string slbl = s->GetLabel();
                if (slbl[0] == 'S') //got an original species leaf
                {
                    Node* sg = s->AddChild();
                    sg->SetLabel("XR" + Util::ToString(i));

                    Node* gs = g->AddChild();
                    gs->SetLabel("XR" + Util::ToString(i));

                    lcaMapping[gs] = sg;
                    done = true;

                    geneRightSpecies[g] = s;
                }
            }
        }
        lca->CloseIterator(itRight);

        forest.push_back(g);
    }

    int AD_prob = rand() % 50 + 25; //between 25-75% chances of having a dup

    //ok, we have a forest.  Now, everything is either S or NAD (no species are shared since we created one specific to each gene)
    //so here we add a couple AD
    for (int i = 0; i < forest.size(); i++)
    {
        Node* g1 = forest[i];
        Node* s1 = lcaMapping[g1];
        for (int j = i + 1; j < forest.size(); j++)
        {
            Node* g2 = forest[j];
            Node* s2 = lcaMapping[g2];

            //they're related...make them AD if we're lucky enough
            if (s1->HasAncestor(s2) || s2->HasAncestor(s1))
            {
                int r = rand() % 100;

                //add a species near the g1left species s.t. g1 and g2 will share a gene of this species
                if (r < AD_prob)
                {
                    Node* s_to_add_to = geneLeftSpecies[g1];
                    if (!s1->HasAncestor(s2))
                        s_to_add_to = geneLeftSpecies[g2];

                    Node* dspecies = s_to_add_to->AddChild();
                    dspecies->SetLabel("AD_" + g1->GetLabel() + "_" + g2->GetLabel());

                    Node* newg1 = g1->AddChild();
                    newg1->SetLabel(dspecies->GetLabel());
                    lcaMapping[newg1] = dspecies;

                    Node* newg2 = g2->AddChild();
                    newg2->SetLabel(dspecies->GetLabel());
                    lcaMapping[newg2] = dspecies;
                }
            }
        }
    }



    //if everything was done correctly, binarizing S
    speciesTree->BinarizeRandomly();
    speciesTree->DeleteSingleChildDescendants();

    string sstr = NewickLex::ToNewickString(speciesTree);
    if (verbose > 0)
        cout<<"S="<<sstr<<endl;

    Node* poly = new Node(false);
    for (int i = 0; i < forest.size(); i++)
    {
        forest[i]->BinarizeRandomly();


        poly->AddSubTree(forest[i]);
    }

    string gstr = NewickLex::ToNewickString(poly);
    if (verbose > 0)
        cout<<"G="<<"="<<gstr<<endl;


    //we have to recreate the species tree, or later on lca mapping will get messed up FOR UNKNOWN REASONS !
    string spNewick = NewickLex::ToNewickString(speciesTree);
    delete speciesTree;

    speciesTree = NewickLex::ParseNewickString(spNewick, true);

    lcaMapping.clear();

    return make_pair(poly, speciesTree);
}




void PolySolverNAD::PerformRandomTest(int k, int verbose)
{

    //here's a good case for you :
    /*
S=((((S8, S0), XL0), XR0)_0, (((XL2, XL4)S1, (((AD_G1_G2, (XL3, ((XR2, AD_G3_G4), XL1)))S5, (S6, S3)), (XR1, XR3)S4)_1_3)_2, XR4)_4);
G=((XL0, XR0)G0, (AD_G1_G2, (XL1, XR1))G1, (XR2, (AD_G1_G2, XL2))G2, (XL3, (XR3, AD_G3_G4))G3, (XL4, (XR4, AD_G3_G4))G4);
SOLVED=((((XL0, XR0)G0, (AD_G1_G2, (XL1, XR1))G1), (XR2, (AD_G1_G2, XL2))G2), ((XL4, (XR4, AD_G3_G4))G4, (XL3, (XR3, AD_G3_G4))G3));


(((XL0, (XL3, XL1))S11, (XR0, XR1)S6)_0_1, (((XL4, ((AD_G4_G8, (AD_G3_G8, XL8)), XR3))S15, ((XL2, ((AD_G2_G4, (XL7, XL5)), (XR4, AD_G2_G3)))S4, (XR5, (XR7, XR2))S13)_2_5_7)_4, (((((XL6, XR8)S10, ((S7, S9), (S0, S16))), (S3, XR6))_6, ((S2, (S8, S5)), S1)), (S14, S17)))_8)_3;
((XL0, XR0)G0, (XL1, XR1)G1, (AD_G2_G4, (AD_G2_G3, (XL2, XR2)))G2, (AD_G2_G3, (AD_G3_G8, (XL3, XR3)))G3, (XR4, ((XL4, AD_G2_G4), AD_G4_G8))G4, (XL5, XR5)G5, (XL6, XR6)G6, (XL7, XR7)G7, (XL8, ((XR8, AD_G3_G8), AD_G4_G8))G8);
(((((XR4, ((XL4, AD_G2_G4), AD_G4_G8))G4, ((XL8, ((XR8, AD_G3_G8), AD_G4_G8))G8, (AD_G2_G3, (AD_G3_G8, (XL3, XR3)))G3)), (AD_G2_G4, (AD_G2_G3, (XL2, XR2)))G2), ((XL1, XR1)G1, (XL7, XR7)G7)), ((XL0, XR0)G0, ((XL5, XR5)G5, (XL6, XR6)G6)));
      */
    pair<Node*, Node*> trees = GetRandomPolytomy(k, verbose);

    Node* poly = trees.first;
    Node* speciesTree = trees.second;

    string gstr = NewickLex::ToNewickString(poly);
    string sstr = NewickLex::ToNewickString(speciesTree);

    vector<Node*> forest = poly->GetChildrenVector();

    unordered_map<Node*, Node*> newMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(poly, speciesTree, "-", 0);


    //note : later on, delete poly also deletes solved
    Node* solved = this->SolvePolytomy(forest, speciesTree, newMapping);

    string solvedStr = NewickLex::ToNewickString(solved);

    if (verbose > 0)
        cout<<"Solved="<<solvedStr<<endl;


    //here's how we count the NADs we created : count the NADs in the forest, count the total NADs, the difference is our fault
    //int forestNads = GeneSpeciesTreeUtil::Instance()->GetNADNodes(poly, speciesTree, newMapping).size();
    /*for (int i = 0; i < forest.size(); i++)
    {
        vector<Node*> fnads = GeneSpeciesTreeUtil::Instance()->GetNADNodes(forest[i], speciesTree, newMapping);
        forestNads += (int)fnads.size();
        fnads.clear();
    }*/

    //vector<Node*> allNads = GeneSpeciesTreeUtil::Instance()->GetNADNodes(solved, speciesTree, newMapping);

    //int mynads = ((int)allNads.size() - forestNads);


    int mynads = this->last_nb_ad_components - 1;

    //cout<<"MYNADS = "<<mynads<<" ("<<this->last_nb_ad_components<<" - 1)"<<endl;

    //find the actual best solution
    SADNADGraph graph;
    graph.BuildGraph(forest, speciesTree, newMapping);

    if (verbose > 1)
        graph.PrintGraph();
    vector<pair<Node*, Node*> > availableSpecs = graph.GetUsefulSpeciationEdges();
    set<pair<Node*, Node*> > chosenSpecs;

    if (verbose > 0)
    {
        cout<<"USEFUL EDGES SIZE="<<availableSpecs.size()<<endl;
        cout<<"INIT AD COMPONENTS="<<graph.GetNbADComponents()<<endl;

        if (verbose > 1)
        {
            const map<Node*, SADNADNode*> nodes = graph.GetNodes();
            for (map<Node*, SADNADNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
            {
                Node* n = (*it).first;

                Node* rep = graph.GetADComponentRepresentantOf(n);

                cout<<n->GetLabel()<<" "<<rep->GetLabel()<<endl;
            }
        }
    }


    bestSolSoFar = 99999;
    bestChoiceSoFar.clear();
    int best = FindBestSolution(graph, availableSpecs, chosenSpecs); //, (k == 7) ? 1 : 0);

    cout<<"BEST - MYNADS = "<<best<<" - "<<mynads<<" = "<<(best - mynads)<<endl;

    if (best < mynads)
    {
        cout<<"S, G, Solved : "<<endl<<sstr<<endl<<gstr<<endl<<solvedStr<<endl;
        //graph.PrintGraph();

        for (set<pair<Node*, Node*> >::iterator it = bestChoiceSoFar.begin(); it != bestChoiceSoFar.end(); it++)
        {
            cout<<(*it).first->GetLabel()<<" "<<(*it).second->GetLabel()<<endl;
        }

        cout<<"------------------------"<<endl;
    }

    delete speciesTree;
    delete poly;


    newMapping.clear();

}






int PolySolverNAD::FindBestSolution(SADNADGraph &graph, vector<pair<Node*, Node*> > availableSpecs, set<pair<Node*, Node*> > chosenSpecs, int verbose)
{


    if (availableSpecs.size() == 0)
    {
        //return the # of AD components
        int nb = graph.GetNbADComponents(chosenSpecs) - 1;

        if (verbose > 0)
        {

            cout<<" RETURNING "<<nb<<endl;

        }

        if (nb < bestSolSoFar)
        {
            bestSolSoFar = nb;
            bestChoiceSoFar = chosenSpecs;
        }

        return nb;
    }

    int bestWith = 0;
    int bestWithout = 0;

    //2 cases : next available spec is chosen, or not

    //EASY CASE : NOT CHOSEN
    vector<pair<Node*, Node*> > localAvailableSpecsWithout(availableSpecs);
    localAvailableSpecsWithout.erase(localAvailableSpecsWithout.begin());
    bestWithout = FindBestSolution(graph, localAvailableSpecsWithout, chosenSpecs, verbose);


    //TOUGH CASE : IT'S CHOSEN.  Then update what's available or not, and chosen or not
    set<pair<Node*, Node*> > localChosen(chosenSpecs);
    pair<Node*, Node*> nextSpec = availableSpecs[0];
    vector<pair<Node*, Node*> > localAvailableSpecs(availableSpecs);
    localChosen.insert(nextSpec);
    localAvailableSpecs.erase(localAvailableSpecs.begin());

    if (verbose > 0)
    {
        cout<<endl<<"CHOOSING AN EDGE "<<endl<<"AVAIL="<<availableSpecs.size()<<" CHOSEN="<<chosenSpecs.size();
        cout<<" Choice="<<nextSpec.first->GetLabel()<<" "<<nextSpec.second->GetLabel();

        for (vector<pair<Node*, Node*> >::iterator vit = availableSpecs.begin(); vit != availableSpecs.end(); vit++)
        {
            cout<<(*vit).first->GetLabel()<<" "<<(*vit).second->GetLabel()<<endl;
        }

        cout<<"CHOSEN="<<endl;

        for (set<pair<Node*, Node*> >::iterator vit = chosenSpecs.begin(); vit != chosenSpecs.end(); vit++)
        {
            cout<<(*vit).first->GetLabel()<<" "<<(*vit).second->GetLabel()<<endl;
        }
    }

    //automatically include every spec in the same clique as nextSpec
    //what nextSpec does is that it "combines" two cliques
    //first, get the 2 cliques the nextSpec vertices are in
    Node* v1 = nextSpec.first;
    Node* v2 = nextSpec.second;
    set<Node*> cliqueV1;
    cliqueV1.insert(v1);
    set<Node*> cliqueV2;
    cliqueV2.insert(v2);

    for (set<pair<Node*, Node*> >::iterator it = localChosen.begin(); it != localChosen.end(); it++)
    {
        pair<Node*, Node*> e = (*it);

        if (e != nextSpec)
        {
            if (e.first == v1)
                cliqueV1.insert(e.second);
            if (e.second == v1)
                cliqueV1.insert(e.first);

            if (e.first == v2)
                cliqueV2.insert(e.second);
            if (e.second == v2)
                cliqueV2.insert(e.first);
        }
    }


    //now, each cliqueV1 and cliqueV2 should be complete to eachother (otherwise this code is not working)
    //The spec edges that link these 2 cliques are chosen "by default", so no need to make them available
    vector<pair<Node*, Node*> >::iterator avit = localAvailableSpecs.begin();
    while (avit != localAvailableSpecs.end())
    {
        pair<Node*, Node*> e = (*avit);

        //NOTE : next spec shouldn't be in localAvailableSpecs
        if ((cliqueV1.find(e.first) != cliqueV1.end() && cliqueV2.find(e.second) != cliqueV2.end()) ||
            (cliqueV2.find(e.first) != cliqueV2.end() && cliqueV1.find(e.second) != cliqueV1.end()))
        {
            avit = localAvailableSpecs.erase(avit);
        }
        else
        {
            avit++;
        }
    }

    //we must add all these clique-joiners edges to the chosen ones.  They are not necessarely all
    //in the availables, so we must do it manually
    for (set<Node*>::iterator c1it = cliqueV1.begin(); c1it != cliqueV1.end(); c1it++)
    {
        Node* n1 = (*c1it);
        for (set<Node*>::iterator c2it = cliqueV2.begin(); c2it != cliqueV2.end(); c2it++)
        {
            Node* n2 = (*c2it);
            localChosen.insert( make_pair(n1, n2) );
        }
    }

    //then, make all speciations that wouldn't grow the clique properly unavailable
    //If W = cliqueV1 union cliqueV2, a node x  is bad w.r.t. W if it has a non-neighbor in W
    //or if x has a bad neighbor in chosen
    //A speciation in available cannot be chosen if
    //- it has one endpoint in W and the other in the bad set


    set<Node*> badNodes;  //first we find the set of such x that are bad
    const map<Node*, SADNADNode*> nodes = graph.GetNodes();
    for (map<Node*, SADNADNode*>::const_iterator itNodes = nodes.begin(); itNodes != nodes.end(); itNodes++)
    {
        Node* n = (*itNodes).first;
        SADNADNode* sad = (*itNodes).second;

        bool itIsBad = false;
        //try to find a non-neighbor of sad in W.  The AD+NAD edges are the complement of the S edges, by the way
        for (set<Node*>::iterator itSad = sad->AD_Neighbors.begin(); itSad != sad->AD_Neighbors.end() && !itIsBad; itSad++)
        {
            if (cliqueV1.find((*itSad)) != cliqueV1.end() || cliqueV2.find((*itSad)) != cliqueV2.end())
            {
                itIsBad = true;
            }
        }
        for (set<Node*>::iterator itSad = sad->NAD_Neighbors.begin(); itSad != sad->NAD_Neighbors.end() && !itIsBad; itSad++)
        {
            if (cliqueV1.find((*itSad)) != cliqueV1.end() || cliqueV2.find((*itSad)) != cliqueV2.end())
            {
                itIsBad = true;
            }
        }

        if (itIsBad)
        {
            badNodes.insert(n);

            //also add every clique nbr of n to the unchoosable
            for (set<pair<Node*, Node*> >::iterator vit = localChosen.begin(); vit != localChosen.end(); vit++)
            {
                if ((*vit).first == n)
                    badNodes.insert((*vit).second);
                if ((*vit).second == n)
                    badNodes.insert((*vit).first);
            }
        }
    }

    if (verbose > 0)
        cout<<" badNodes size="<<badNodes.size();

    //make unavailable all specs with one endpoint in badNodes and one in cliqueNbrs
    //TODO : this could be incorporated in the previous loop...here we keep it simple
    avit = localAvailableSpecs.begin();
    while (avit != localAvailableSpecs.end())
    {
        pair<Node*, Node*> e = (*avit);

        if ((badNodes.find(e.first) != badNodes.end() && cliqueV1.find(e.second) != cliqueV1.end()) ||
            (badNodes.find(e.second) != badNodes.end() && cliqueV1.find(e.first) != cliqueV1.end()) ||
            (badNodes.find(e.first) != badNodes.end() && cliqueV2.find(e.second) != cliqueV2.end()) ||
            (badNodes.find(e.second) != badNodes.end() && cliqueV2.find(e.first) != cliqueV2.end()))
        {
            avit = localAvailableSpecs.erase(avit);
        }
        else
        {
            avit++;
        }
    }

    if (verbose > 0)
    {
        cout<<endl<<"AfterAvail="<<localAvailableSpecs.size()<<endl;
    }

    bestWith = FindBestSolution(graph, localAvailableSpecs, localChosen, verbose);


    return min(bestWithout, bestWith);
}
