#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <streambuf>

using namespace std;

void Split(string const& to_split, char delim, vector<string>& tokens)
{
    tokens.clear();
    string item;
    stringstream ss(to_split);
    while (getline(ss, item, delim))
    {
        tokens.push_back(item);
    }
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        cerr << "Not enough parameters" << endl;
        return 2;
    }

    ifstream chromosome_file(argv[1]);
    if (!chromosome_file.good())
    {
        cerr << "Bad chromosome file" << endl;
        return 1;
    }

    ifstream cpg_islands(argv[2]);
    if (!cpg_islands.good())
    {
        cerr << "Bad cpg islands file" << endl;
        return 1;
    }

    ofstream cpg_training_set("cpg_training_set");
    ofstream non_cpg_training_set("non_cpg_training_set");

    string chromosome;
    string line;
    while (chromosome_file)
    {
        getline(chromosome_file, line);
        if (line[0] != '>')
        {
            chromosome += line;
        }
    }
    chromosome_file.close();
    transform(chromosome.begin(), chromosome.end(), chromosome.begin(), ::tolower);

    string island;
    string non_island;
    vector<string> tokens;
    size_t island_begin;
    size_t island_end;
    size_t prev_end = 0;
    while (cpg_islands)
    {
        getline(cpg_islands, line);

        if (line.empty())
        {
            continue;
        }

        Split(line, '\t', tokens);
        island_begin = atoi(tokens[1].c_str());
        island_end = atoi(tokens[2].c_str());

        non_island = chromosome.substr(prev_end, island_begin - prev_end);
        island = chromosome.substr(island_begin, island_end - island_begin + 1);
        prev_end = island_end + 1;

        non_cpg_training_set << non_island << 'n' << endl;
        cpg_training_set << island << 'n' << endl;
    }

    non_island = chromosome.substr(prev_end);
    non_cpg_training_set << non_island << 'n' << endl;

    cpg_islands.close();
    non_cpg_training_set.close();
    cpg_training_set.close();
    return 0;
}
