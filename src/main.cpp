#include "precompiled.h"
#include "HMMlib/hmm_table.hpp"
#include "HMMlib/hmm_vector.hpp"
#include "HMMlib/hmm.hpp"

using hmmlib::HMM;
using hmmlib::HMMMatrix;
using hmmlib::HMMVector;

#define DEBUG

void Sequence_Stats(string sequence)
{
    uint cg_counter = 0;
    uint c_counter = 0;
    uint g_counter = 0;
    uint cpg_counter = 0;
    uint i = 0;
    for (; i < sequence.length() - 1; ++i)
    {
        if (sequence[i] == 'c')
        {
            c_counter += 1;
            cg_counter += 1;
        }

        if (sequence[i] == 'g')
        {
            g_counter += 1;
            cg_counter += 1;
        }

        if (sequence[i] == 'c' && sequence[i + 1] == 'g')
        {
            cpg_counter += 1;
        }
    }

    if (sequence[i] == 'c')
    {
        c_counter += 1;
        cg_counter += 1;
    }

    if (sequence[i] == 'g')
    {
        g_counter += 1;
        cg_counter += 1;
    }

    cerr << "Sequence length: " << sequence.length() << endl;
    cerr << "CG percenage: " << std::setprecision(3) << cg_counter/(double)sequence.length()*100 << "%" << endl;
    cerr << "Observed to expected CpG ratio: " << std::setprecision(3)
         << cpg_counter/(double)(c_counter*g_counter)*sequence.length()*100 << "%" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Print_Cpg(sequence hidden_sequence)
{
    uint start_cpg = 0;
    bool is_open_cpg = false;
    for (ull i = 0; i < hidden_sequence.size(); ++i)
    {
        if (hidden_sequence[i] < 4)
        {
            if (!is_open_cpg)
            {
                start_cpg = i;
                is_open_cpg = true;
            }
        }
        else
        {
            if (is_open_cpg)
            {
                cout << "chr1\t" << start_cpg << "\t" << i << endl;
                is_open_cpg = false;
            }
        }
    }
}

void Set_HMM_Parameters(shared_ptr< HMMVector<double> > initial_probabilities_sptr,
                       shared_ptr< HMMMatrix<double> > transition_probabilities_sptr,
                       shared_ptr< HMMMatrix<double> > emission_probabilities_sptr,
                        ifstream& parameters_txt)
{
    HMMVector<double>& initial_probabilities = *initial_probabilities_sptr;
    HMMMatrix<double>& transition_probabilities = *transition_probabilities_sptr;
    HMMMatrix<double>& emission_probabilities = *emission_probabilities_sptr;

    string line;

    getline(parameters_txt, line);
    // initial probabilities
    double next_value;
    for (ushort i = 0; i < initial_probabilities.get_size(); ++i)
    {
        parameters_txt >> next_value;
        initial_probabilities(i) = next_value;
    }

    getline(parameters_txt, line);
    getline(parameters_txt, line);
    //transition probabilities
    for (ushort i = 0; i < initial_probabilities.get_size(); ++i)
    {
        for (ushort j = 0; j < initial_probabilities.get_size(); ++j)
        {
            parameters_txt >> next_value;
            transition_probabilities(i, j) = next_value;
        }
    }

    getline(parameters_txt, line);
    getline(parameters_txt, line);
    //emission probabilities
    for (ushort i = 0; i < emission_probabilities.get_no_rows(); ++i)
    {
        for (ushort j = 0; j < initial_probabilities.get_size(); ++j)
        {
            parameters_txt >> next_value;
            emission_probabilities(i, j) = next_value;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    string mode = "predict";

    bool parameters_expectation = false;
    ifstream parameters_file;
    //command line arguments handling
    for (int arg = 1; arg < argc; ++arg)
    {
        if (parameters_expectation)
        {
            parameters_expectation = false;
            parameters_file.open(argv[arg]);
        }

        if (!strcmp(argv[arg], "--help") || !strcmp(argv[arg], "-h"))
        {
            cerr << "This program searches for CpG islands in genome" << endl;
            cerr << "Usage: cat fasta file with genome and pass output to program" << endl;
            cerr << "Example: cat genome.fa | ./HMM_CpG" << endl;
            return 0;
        }

        if (!strcmp(argv[arg], "--version") || !strcmp(argv[arg], "-v"))
        {
            cerr << "1.0.0" << endl;
            return 0;
        }

        if (!strcmp(argv[arg], "--train") || !strcmp(argv[arg], "-t"))
        {
            mode = "train";
        }

        if (!strcmp(argv[arg], "--parameters") || !strcmp(argv[arg], "-p"))
        {
            parameters_expectation = true;
        }
    }

    if (!parameters_file.good() || !parameters_file.is_open())
    {
        cerr << "Can't read file" << endl;
        return 1;
    }

    vector<int> nucleotides_mapping(127, -1);
    nucleotides_mapping['a'] = 0;
    nucleotides_mapping['c'] = 1;
    nucleotides_mapping['g'] = 2;
    nucleotides_mapping['t'] = 3;

//    vector<char> state_mapping(8);
//    state_mapping[0] = 'A';
//    state_mapping[1] = 'C';
//    state_mapping[2] = 'G';
//    state_mapping[3] = 'T';
//    state_mapping[4] = 'a';
//    state_mapping[5] = 'c';
//    state_mapping[6] = 'g';
//    state_mapping[7] = 't';

#ifdef DEBUG
    cerr << "DEBUG" << endl;
    freopen("chr1.fa", "rt", stdin);
#endif

    //Sequence_Stats(observed);

    const uint number_of_states = 8;
    const uint alphabet_size = 4;

    shared_ptr< HMMVector<double> > initial_probabilities(new HMMVector<double>(number_of_states));
    shared_ptr< HMMMatrix<double> > transition_probabilities(new HMMMatrix<double>(number_of_states, number_of_states));
    shared_ptr< HMMMatrix<double> > emission_probabilities(new HMMMatrix<double>(alphabet_size, number_of_states));

    Set_HMM_Parameters(initial_probabilities, transition_probabilities, emission_probabilities, parameters_file);

    parameters_file.close();

    HMM<double>* hmm_ptr;
    try
    {
        hmm_ptr = new HMM<double>(initial_probabilities, transition_probabilities, emission_probabilities);
    }
    catch(const char* ex)
    {
        cerr << ex << endl;
        return 1;
    }

    shared_ptr< HMM<double> > hmm_sptr(hmm_ptr);
    HMM<double> hmm = *hmm_sptr;

    string observed;

    string line;
    while (!cin.eof())
    {
        cin >> line;
        observed += line;
    }

    transform(observed.begin(), observed.end(), observed.begin(), ::tolower);

    sequence observed_sequence;
    observed_sequence.reserve(observed.length());

    for (ull i = 0; i < observed.length(); ++i)
    {
        if (observed[i] != 'a' && observed[i] != 'c' && observed[i] != 'g' && observed[i] != 't' && observed[i] != 'n')
        {
            cerr << "Unexpected character. Possible characters are a, c, g, t, n (in any register)" << endl;
            throw "Unexpected character. Possible characters are a, c, g, t, n (in any register)";
        }

        if (observed[i] != 'n')
        {
            observed_sequence.push_back(nucleotides_mapping[observed[i]]);
        }
        else
        {
            if (observed_sequence.size() > 0)
            {
                if (mode == "predict")
                {
                    sequence hidden_sequence;
                    hmm.Predict(observed_sequence, hidden_sequence);
                    Print_Cpg(hidden_sequence);
                }
                else
                {
                    hmm.Train(observed_sequence);
                }

                observed_sequence.clear();
            }
        }
    }

    //save new parameters
    if (mode == "train")
    {
        hmm.Save_Parameters();
    }

    return 0;
}
