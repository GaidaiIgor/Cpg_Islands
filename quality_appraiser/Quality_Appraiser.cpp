#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstring>

using namespace std;
typedef unsigned int uint;

ifstream tested_cpg;
ifstream true_cpg;

void Split(string const& to_split, char delim, vector<string>& tokens)
{
    tokens.clear();
    string item;
    stringstream splitted_stream(to_split);
    while (getline(splitted_stream, item, delim))
    {
        tokens.push_back(item);
    }
}

void Count_Hits()
{
    string tested_line;
    string true_line;
    vector<string> tokens;
    size_t tested_begin;
    size_t tested_end;
    size_t true_begin;
    size_t true_end;
    size_t false_islands = 0;
    size_t hits = 0;
    size_t missed_islands = 0;
    size_t total_islands = 0;
    size_t total_suggestions = 0;
    bool update_true = true;
    bool update_tested = true;
    while (tested_cpg || true_cpg)
    {
        if (update_tested)
        {
            update_tested = false;
            getline(tested_cpg, tested_line);
            Split(tested_line, '\t', tokens);
            tested_begin = atoi(tokens[1].c_str());
            tested_end = atoi(tokens[2].c_str());
            ++total_suggestions;
        }

        if (update_true)
        {
            update_true = false;
            getline(true_cpg, true_line);
            Split(true_line, '\t', tokens);
            true_begin = atoi(tokens[1].c_str());
            true_end = atoi(tokens[2].c_str());
            ++total_islands;
        }

        if (true_line.empty() && tested_line.empty())
        {
            break;
        }

        if (true_line.empty())
        {
            ++false_islands;
            update_tested = true;
            continue;
        }

        if (tested_line.empty())
        {
            ++missed_islands;
            update_true = true;
            continue;
        }

        if (tested_end < true_begin)
        {
            update_tested = true;
            ++false_islands;
        }

        if ((tested_end >= true_begin && tested_end <= true_end) || (tested_begin >= true_begin && tested_begin <= true_end)
                || (tested_end >= true_end && tested_begin <= true_begin))
        {
            ++hits;
            update_true = true;
            update_tested = true;
        }

        if (tested_begin > true_end)
        {
            ++missed_islands;
            update_true = true;
        }
    }

    cout << "Total islands: " << total_islands << endl;
    cout << "Total suggestions: " << total_suggestions << endl;
    cout << "Hits: " << hits << " (" << setprecision(2) << hits/(double)total_islands*100 << "%)" << endl;
    cout << "Missed: " << missed_islands << " (" << setprecision(2) << missed_islands/(double)total_islands*100 << "%)" << endl;
    cout << "False islands: " << false_islands << " (" << setprecision(2) << false_islands/(double)total_suggestions*100 << "%)" << endl;
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        cout << "Not enough arguments" << endl;
        return 1;
    }

    tested_cpg.open(argv[1]);
    true_cpg.open(argv[2]);

    if (!tested_cpg.good())
    {
        cout << "Bad tested cpg file" << endl;
        return 2;
    }
    if (!true_cpg.good())
    {
        cout << "Bad true cpg file" << endl;
        return 2;
    }

    if (!strcmp(argv[3], "--hits"))
    {
        Count_Hits();
    }

    tested_cpg.clear();
    true_cpg.clear();
    tested_cpg.seekg(0, tested_cpg.beg);
    true_cpg.seekg(0, true_cpg.beg);

    uint total_cpg_length = 0;
    uint total_tested_length = 0;

    string tested_line;
    string true_line;
    uint tested_start_interval = 0;
    uint tested_end_interval = 0;
    uint prev_tested_end = 0;
    uint true_start_interval = 0;
    uint true_end_interval = 0;
    uint tested_not_in_true = 0;
    uint common = 0;
    uint true_not_in_tested = 0;
    vector<string> tested_line_tokens;
    vector<string> true_line_tokens;
    bool update_true_intervals = true;
    bool update_test_intervals = true;

    while (tested_cpg || true_cpg)
    {
        if (update_test_intervals)
        {
            prev_tested_end = tested_end_interval;
            update_test_intervals = false;
            getline(tested_cpg, tested_line);
            Split(tested_line, '\t', tested_line_tokens);
            tested_start_interval = atoi(tested_line_tokens[1].c_str());
            tested_end_interval = atoi(tested_line_tokens[2].c_str());

            if (!tested_line.empty())
            {
                total_tested_length += tested_end_interval - tested_start_interval + 1;
            }
        }

        if (update_true_intervals)
        {
            update_true_intervals = false;
            getline(true_cpg, true_line);
            Split(true_line, '\t', true_line_tokens);
            true_start_interval = atoi(true_line_tokens[1].c_str());
            true_end_interval = atoi(true_line_tokens[2].c_str());

            if (!true_line.empty())
            {
                total_cpg_length += true_end_interval - true_start_interval + 1;
            }
        }

        if (true_line.empty() && tested_line.empty())
        {
            break;
        }

        if (true_line.empty())
        {
            tested_not_in_true += tested_end_interval - tested_start_interval + 1;
            update_test_intervals = true;
            continue;
        }

        if (tested_line.empty())
        {
            true_not_in_tested += true_end_interval - true_start_interval + 1;
            update_true_intervals = true;
            continue;
        }

        //  ####  ####
        // #############
        if (tested_start_interval <= true_end_interval && prev_tested_end > true_start_interval)
        {
            true_not_in_tested += tested_start_interval - prev_tested_end - 1;
        }

        //   ####      #####
        // ##########
        if (tested_start_interval > true_end_interval && prev_tested_end < true_end_interval && prev_tested_end >= true_start_interval)
        {
            true_not_in_tested += true_end_interval - prev_tested_end;
        }

        //  #####               #######
        //         ########
        if (tested_start_interval > true_end_interval && prev_tested_end < true_start_interval)
        {
            true_not_in_tested += true_end_interval - true_start_interval + 1;
            update_true_intervals = true;
        }

        // ###   ###
        //     #########
        if (prev_tested_end < true_start_interval && tested_start_interval > true_start_interval && tested_start_interval < true_end_interval)
        {
            true_not_in_tested += tested_start_interval - true_start_interval;
        }

        //tested_not_in_true

        // #####
        //         #######
        if (tested_end_interval < true_start_interval)
        {
            tested_not_in_true += tested_end_interval - tested_start_interval + 1;
            update_test_intervals = true;
        }

        // #####
        //    #######
        if (tested_end_interval >= true_start_interval && tested_start_interval < true_start_interval && tested_end_interval <= true_end_interval)
        {
            tested_not_in_true += true_start_interval - tested_start_interval;
            common += tested_end_interval - true_start_interval + 1;
            update_test_intervals = true;
        }

        // ##########
        //    ###
        if (tested_start_interval < true_start_interval && tested_end_interval > true_end_interval)
        {
            tested_not_in_true += true_start_interval - tested_start_interval;
            common += true_end_interval - true_start_interval + 1;
            tested_start_interval = true_end_interval + 1;
            update_true_intervals = true;
            continue;
        }

        //     ###
        // #########
        if (tested_end_interval <= true_end_interval && tested_start_interval >= true_start_interval)
        {
            common += tested_end_interval - tested_start_interval + 1;
            update_test_intervals = true;
        }

        //           #######
        // ######
        if (tested_start_interval > true_end_interval && tested_start_interval > true_end_interval)
        {
            update_true_intervals = true;
        }

        //    ######
        // ######
        if (tested_end_interval > true_end_interval && tested_start_interval <= true_end_interval)
        {
            common += true_end_interval - tested_start_interval + 1;
            update_true_intervals = true;
            tested_start_interval = true_end_interval + 1;
        }
    }

    cout << "Total true length: " << total_cpg_length << endl;
    cout << "Total tested length: " << total_tested_length << endl;
    cout << "Total tested not in true: " << tested_not_in_true << " (" << setprecision(2) << tested_not_in_true/(double)total_tested_length*100 << "%)" << endl;
    cout << "Total common1: " << common << " (" << setprecision(2) << common/(double)total_tested_length*100 << "%)" << endl;
    cout << "Total common2: " << common << " (" << setprecision(2) << common/(double)total_cpg_length*100 << "%)" << endl;
    cout << "Total true not in tested: " << true_not_in_tested << " (" << setprecision(2) << true_not_in_tested/(double)total_cpg_length*100 << "%)" << endl;

    return 0;
}
