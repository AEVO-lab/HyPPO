
#include "newicklex.h"




Node* NewickLex::ParseNewickString(string& str, bool maintainTreeInfo)
{
    Node* root = new Node(maintainTreeInfo);
    int pos = str.find_last_of(')');
    int lastcolonpos = str.find_last_of(';');


    if (pos != string::npos)
    {
        string label = "";
        if (lastcolonpos > pos)
        {
            label = str.substr(pos + 1, lastcolonpos - pos - 1);
        }
        else
        {
            label = Util::RTrim(str.substr(pos + 1));
        }
        root->SetLabel(label);

        ReadNodeChildren(str, pos - 1, root);
    }
    else
    {
        //this should be the root label
        root->SetLabel(str.substr(0, str.length() - 1));
    }

    return root;
}


string NewickLex::ToNewickString(Node* root, bool addBranchLengthToLabel, bool addInternalNodesLabel)
{
    string str;
    WriteNodeChildren(str, root, addBranchLengthToLabel, addInternalNodesLabel);
    str += ";";
    return str;
}



int NewickLex::ReadNodeChildren(string &str, int revstartpos, Node* curNode)
{
    bool ok = true;
    int pos = revstartpos;

    //cout<<"POS="<<pos<<endl;

    while (ok && pos >= 0)
    {
        int closepos = str.find_last_of(')', pos);
        int openpos = str.find_last_of('(', pos);
        int commapos = str.find_last_of(',', pos);

        int maxpos = max(max(closepos, openpos), commapos);

        if (maxpos == string::npos)
        {
            ok = false;
            pos = -1;
        }
        else
        {

            string lbl = "";

            lbl = Util::Trim(str.substr(maxpos + 1, pos - maxpos));

            if (maxpos == closepos)
            {
                //cout<<"Close="<<maxpos<<"      LBL="<<lbl<<endl;
                Node* newNode = curNode->InsertChild(0);

                ParseLabel(newNode, lbl);
                //newNode->SetLabel(lbl);
                pos = ReadNodeChildren(str, maxpos - 1, newNode);

                while (str[pos] != ',' && str[pos] != '(')
                    pos--;

                if (str[pos] == ',')
                    pos--;
            }
            else if (maxpos == commapos)
            {
                //cout<<"Comma="<<maxpos<<"      LBL="<<lbl<<endl;

                int ptcomma = str.find_first_of(',', maxpos + 1);
                int ptopen = str.find_first_of('(', maxpos + 1);
                int ptclose = str.find_first_of(')', maxpos + 1);

                if ((ptcomma < ptopen && ptcomma < ptclose) ||
                    (ptopen == string::npos || ptopen > ptclose))
                {
                    Node* newNode = curNode->InsertChild(0);
                    ParseLabel(newNode, lbl);
                    //newNode->SetLabel(lbl);
                }
                pos = maxpos - 1;
            }
            else if (maxpos == openpos)
            {
                //cout<<"Open="<<maxpos<<"      LBL="<<lbl<<endl;

                //EDIT ML AUG 2012 : an opening parenthese creates a node only if followed by a ,
                // THIS WAS DONE DURING A PHASE OF TIREDNESS
                // IF SOMETHING IS BUGGY, IT'S PROBABLY AROUND HERE
                //NOTE ML AUG 2013 : it does seem to be holding up, even after extensive use
                int ptcomma = str.find_first_of(',', maxpos + 1);
                int ptopen = str.find_first_of('(', maxpos + 1);
                int ptclose = str.find_first_of(')', maxpos + 1);

                if (ptcomma != string::npos &&
                        (ptcomma < ptclose || ptclose == string::npos) &&
                        (ptcomma < ptopen || ptopen == string::npos))
                {
                    Node* newNode = curNode->InsertChild(0);
                    ParseLabel(newNode, lbl);
                    //newNode->SetLabel(lbl);
                }
                pos = maxpos - 1;
                ok = false;
            }
        }

    }

    return pos;
}

void NewickLex::WriteNodeChildren(string &str, Node* curNode, bool addBranchLengthToLabel, bool addInternalNodesLabel)
{
    if (curNode->IsLeaf())
    {

        str += curNode->GetLabel();

        if (addBranchLengthToLabel && !curNode->IsRoot())
            str += ":" + Util::ToString(curNode->GetBranchLength());
    }
    else
    {
        str += "(";
        for (int i = 0; i < curNode->GetNbChildren(); i++)
        {
            if (i != 0)
                str += ", ";

            Node* child = curNode->GetChild(i);
            WriteNodeChildren(str, child, addBranchLengthToLabel, addInternalNodesLabel);
        }

        str += ")";

        if (addInternalNodesLabel)
            str += curNode->GetLabel();

        if (addBranchLengthToLabel && !curNode->IsRoot())
            str += ":" + Util::ToString(curNode->GetBranchLength());
    }


}




void NewickLex::ParseLabel(Node* node, string label)
{
    string sublabel = "";
    int pos = label.find("[");
    if (pos != string::npos)
    {
        sublabel = label.substr(pos);
        label = label.substr(0, pos - 1);
    }

    vector<string> sz = Util::Split(label, ":");

    if (sz.size() == 1)
    {
        //If label is just a double, we consider it's a branch length
        /*if (Util::IsDouble(sz[0]))
        {
            node->SetBranchLength(Util::ToDouble(sz[0]));
        }
        else*/
        {
            node->SetLabel(label);
        }
    }
    else
    {
        node->SetLabel(sz[0]);
        node->SetBranchLength(Util::ToDouble(sz[1]));
    }

    node->SetLabel(node->GetLabel() + sublabel);
}




string NewickLex::GetCaterpillarNewick(set<string> labels)
{
    string newick = "";

    set<string>::iterator it;
    int cpt = 0;
    for (it = labels.begin(); it != labels.end(); it++)
    {
        if (cpt == 0)
            newick = *it;
        else if (cpt == 1)
            newick = "(" + newick + "," + *it + ")";
        else
        {
            newick = "(" + newick + "," + *it + ")";
        }

        cpt += 1;
    }

    return newick;
}
