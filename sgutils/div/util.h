#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <iostream>
#include <fstream>
#include "define.h"
#include <algorithm>
#include <cctype>
#include <cwctype>
#include <sstream>
#include <vector>
#include <iomanip>

#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

using namespace std;

inline bool caseInsCharCompareN(char a, char b) {
   return(toupper(a) == toupper(b));
}


/**
  Various useful static methods
  **/
class Util
{
public:

    /**
      Converts n in base 2, returning the string of ones and zeros
      **/
    static string UInt64ToBinary(uint64 n)
    {
        string str;
        uint64 buf = 1;

        for (int i = 63; i >= 0; i--)
        {
            if (n & (buf << (uint64)(i)))
                str += "1";
            else
                str += "0";
        }

        return str;
    }

    /**
      Simply outputs n as a binary string, with an optional message beforehand
      **/
    static void DumpUInt64Bin(uint64 n,string msg = "")
    {
        cout<<msg<<UInt64ToBinary(n)<<endl;
    }


    /**
      Trims s on the right from any single character in delimiters
      **/
    static string RTrim(string s, string delimiters = " \f\n\r\t\v")
    {
        if (s.length() == 0)
            return s;

      return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
    }


    /**
      Trims s on the left from any single character in delimiters
      **/
    static string LTrim(string s, string delimiters = " \f\n\r\t\v" )
    {
      if (s.length() == 0)
        return s;

      return s.substr( s.find_first_not_of( delimiters ) );
    }

    /**
      Trims s on the left and right from any single character in delimiters
      **/
    static string Trim(std::string s, string delimiters = " \f\n\r\t\v" )
    {
      return LTrim( RTrim( s, delimiters ), delimiters );
    }


    /**
      Simply outputs str on stdout
      **/
    static void DebugOut(string str)
    {
        cout<<str<<endl;
    }

    /**
      Returns true if s1 = s2, false otherwise
      **/
    static bool Streq(const string& s1, const string& s2) {
       return((s1.size( ) == s2.size( )) &&
              equal(s1.begin( ), s1.end( ), s2.begin( ), caseInsCharCompareN));
    }


    /**
      Replaces all occurences of find by rep in str.
      **/
    static string ReplaceAll( const string& str, const string& find, const string& rep ) {
        if ( str.empty() || find.empty() || find == rep || str.find(find) == string::npos ) {
            return str;
        }
        ostringstream build_it;
        size_t i = 0;
        for ( size_t pos; ( pos = str.find( find, i ) ) != string::npos; ) {
            build_it.write( &str[i], pos - i );
            build_it << rep;
            i = pos + find.size();
        }
        if ( i != str.size() ) {
            build_it.write( &str[i], str.size() - i );
        }
        return build_it.str();
    }

    static string ToString(int v)
    {
        stringstream ss;
        ss << v;
        return ss.str();
    }

    static string ToString(double v)
    {
        stringstream ss;

        ss.precision(10);

        ss  << fixed  << v;
        return ss.str();
    }

    static double ToDouble(string s)
    {
        double d;
        stringstream ss;
        ss.str(s.c_str());
        ss >> d;

        return d;
    }


    static bool IsDouble(string s)
    {
        //TODO : this might fail
        double d = 0;
        istringstream ss;
        ss.str(s.c_str());
        ss >> d;

        return (!ss.fail() && ss.eof());
    }


    static int ToInt(string s)
    {
        int i;
        stringstream ss;
        ss.str(s.c_str());
        ss >> i;

        return i;

    }


    /**
      Splits str by the splitter, returns a vector of all obtained strings
      */
    static vector<string> Split(string str, string splitter, bool includeEmpty = true)
    {
        vector<string> v;

        if (splitter.length() > str.length())
        {
            v.push_back(str);
            return v;
        }

        int current = 0;
        int next = -1 * splitter.length();


        do
        {
          current = next + splitter.length();
          next = str.find( splitter, current );

          if (next == string::npos)
          {
              if (includeEmpty || (current < str.length()))
              {
                v.push_back(str.substr( current ));
              }
          }
          else
          {
              if (includeEmpty || (next - current > 0))
              {
                v.push_back(str.substr( current, next - current ));
              }
          }
        }
        while (next != string::npos);

        return v;
    }


    /**
      Inserts splitter after each nbchars characters in str
      **/
    static string SplitByLength(string str, int nbchars, string splitter = "\n")
    {
        string out = "";

        int pos = 0;

        while (pos < str.length())
        {
            if (str.length() > pos + nbchars)
            {
                out += str.substr(pos, nbchars) + splitter;

            }
            else
            {
                out += str.substr(pos);
            }
            pos += nbchars;
        }

        return out;
    }

    static string ToLower(string str)
    {
        transform(str.begin(), str.end(),str.begin(), ::tolower );
        return str;
    }

    static string ToUpper(string str)
    {
        transform(str.begin(), str.end(),str.begin(), ::toupper );
        return str;
    }

    /**
      Just doubles every apostrophe in s, then wraps s with apostrophes
      **/
    static string DBEscape(string s)
    {
        return "'" + Util::ReplaceAll(s, "'", "''") + "'";
    }

    /**
      Prepares a list of values to be used in a query of the type
      "WHERE some_field IN (sz[0], sz[1], ...)"
      Usage :
      @code
      string q = "SELECT * FROM table WHERE field IN (" + Util::ToInstr(my_vector) + ")";
      @endcode
      **/
    static string ToInstr(vector<string> sz)
    {
        string instr = "";

        for (int i = 0; i < sz.size(); i++)
        {
            if (instr != "")
                instr +=",";
            instr += Util::DBEscape(sz[i]);
        }

        return instr;
    }



    static string GetSubstringBefore(string s, string separator)
    {
        int pos = s.find_first_of(separator);

        if (pos != string::npos)
            return s.substr(0, pos);

        return s;
    }

    static string GetSubstringAfter(string s, string separator)
    {
        int pos = s.find_last_of(separator);

        if (pos != string::npos)
            return s.substr(pos + 1);

        return s;
    }


    static bool FileExists(string filename)
    {
      ifstream ifile(filename.c_str());
      if (ifile)
          return true;
      else
        return false;
    }

    static vector<string> GetFileLines(string filename)
    {
        string superstr = Util::GetFileContent(filename);

        vector<string> unfilteredLines = Util::Split(superstr, "\n");

        vector<string> lines;
        for (int i = 0; i < unfilteredLines.size(); i++)
        {
            if (unfilteredLines[i] != "")
                lines.push_back(unfilteredLines[i]);
        }

        return lines;

    }

    static string GetFileContent(string filename)
    {
        std::ifstream sifs(filename);
        std::string spcontent( (std::istreambuf_iterator<char>(sifs) ),
                             (std::istreambuf_iterator<char>()    ));
        sifs.close();

        return spcontent;
    }


    static void WriteFileContent(string filename, string content, bool append = false)
    {
        ofstream outfile;
        if (!append)
            outfile.open (filename);
        else
            outfile.open (filename, ios_base::app | ios_base::out);
        outfile<<content;
        outfile.close();
    }

    static string GetPathFilename(string fullpath)
    {
        return Util::GetSubstringAfter(fullpath, "/");
    }


    static string GetFileLine(string filename, int lineIndex)
    {
        ifstream file(filename);
        string line;
        int line_number = 0;

        string myLine = "";

        bool wereDone = false;
        while (std::getline(file, line) && !wereDone)
        {
            if (line_number == lineIndex)
            {
                 myLine = line;
                 wereDone = true;
            }
            line_number++;
        }

        return myLine;
    }


    static bool StartsWith (string fullString, string begin) {
        if (fullString.length() >= begin.length()) {
            return (0 == fullString.compare (0, begin.length(), begin));
        } else {
            return false;
        }
    }


    static bool EndsWith (string fullString, string ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
        } else {
            return false;
        }
    }



    static bool DirectoryExists(const std::string& path)
    {
    #if defined(_WIN32)
        struct _stat info;
        if (_stat(path.c_str(), &info) != 0)
        {
            return false;
        }
        return (info.st_mode & _S_IFDIR) != 0;
    #else
        struct stat info;
        if (stat(path.c_str(), &info) != 0)
        {
            return false;
        }
        return (info.st_mode & S_IFDIR) != 0;
    #endif
    }

    static bool CreateDirectory(const std::string& path)
    {
    #if defined(_WIN32)
        int ret = _mkdir(path.c_str());
    #else
        mode_t mode = 0755;
        int ret = mkdir(path.c_str(), mode);
    #endif
        if (ret == 0)
            return true;

        switch (errno)
        {
        case ENOENT:
            // parent didn't exist, try to create it
            {
                int pos = path.find_last_of('/');
                if (pos == std::string::npos)
    #if defined(_WIN32)
                    pos = path.find_last_of('\\');
                if (pos == std::string::npos)
    #endif
                    return false;
                if (!CreateDirectory( path.substr(0, pos) ))
                    return false;
            }
            // now, try to create again
    #if defined(_WIN32)
            return 0 == _mkdir(path.c_str());
    #else
            return 0 == mkdir(path.c_str(), mode);
    #endif

        case EEXIST:
            // done!
            return DirectoryExists(path);

        default:
            return false;
        }
    }



};

#endif // UTIL_H
