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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Train(sequence& observed_sequence, HMM<double>& hmm)
{
    HMMMatrix<double> forward_dynamic(observed_sequence.size(), hmm.get_no_states());
    HMMVector<double> scales(observed_sequence.size());
    cerr << "Running forward" << endl;
    hmm.forward(observed_sequence, scales, forward_dynamic);

    cerr << "Running backward" << endl;
    HMMMatrix<double> backward_dynamic(observed_sequence.size(), hmm.get_no_states());
    hmm.backward(observed_sequence, scales, backward_dynamic);

    cerr << "Running Baum-Welch" << endl;
    shared_ptr< HMMVector<double> > new_initial_probabilities_sptr(new HMMVector<double>(hmm.get_no_states()));
    shared_ptr< HMMMatrix<double> > new_transition_probabilities_sptr(new HMMMatrix<double>(hmm.get_no_states(),
                                                                                                  hmm.get_no_states()));
    shared_ptr< HMMMatrix<double> > new_emission_probabilities_sptr(new HMMMatrix<double>(hmm.get_alphabet_size(),
                                                                                                hmm.get_no_states()));
    hmm.baum_welch(observed_sequence, forward_dynamic, backward_dynamic, scales, *new_initial_probabilities_sptr,
                   *new_transition_probabilities_sptr, *new_emission_probabilities_sptr);

    hmm.Set_Initial_Probabilities(new_initial_probabilities_sptr);
    hmm.Set_Transitions_Probabilities(new_transition_probabilities_sptr);
    hmm.Set_Emission_Probabilities(new_emission_probabilities_sptr);

    cerr << "Iteration complete" << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Set_HMM_Parameters(shared_ptr< HMMVector<double> > initial_probabilities_sptr,
                       shared_ptr< HMMMatrix<double> > transition_probabilities_sptr,
                       shared_ptr< HMMMatrix<double> > emission_probabilities_sptr)
{
    HMMVector<double>& initial_probabilities = *initial_probabilities_sptr;
    HMMMatrix<double>& transition_probabilities = *transition_probabilities_sptr;
    HMMMatrix<double>& emission_probabilities = *emission_probabilities_sptr;

    ifstream parameters_txt("parameters.txt");

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

    parameters_txt.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Predict(sequence& observed_sequence, HMM<double>& hmm)
{
    double log_likelihood;
    cerr << "Running viterbi" << endl;
    sequence hidden_sequence;
    hidden_sequence.resize(observed_sequence.size());
    log_likelihood = hmm.viterbi(observed_sequence, hidden_sequence);
    cerr << "\nLog likelihood of hidden sequence: " << log_likelihood << endl;

    Print_Cpg(hidden_sequence);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Save_Parameters(HMM<double>& hmm)
{
    ofstream parameters_txt("new_parameters.txt");

    parameters_txt << "Initial probabilities" << endl;
    HMMVector<double> const& initial_probabilities = hmm.get_initial_probs();
    for (uint i = 0; i < initial_probabilities.get_size(); ++i)
    {
        parameters_txt << initial_probabilities(i) << endl;
    }

    parameters_txt << "Transition probabilities" << endl;
    HMMMatrix<double> const& transition_probabilities = hmm.get_trans_probs();
    for (uint i = 0; i < transition_probabilities.get_no_rows(); ++i)
    {
        for (uint j = 0; j < transition_probabilities.get_no_columns(); ++j)
        {
            parameters_txt << transition_probabilities(i, j) << "\t";
        }
        parameters_txt << endl;
    }

    parameters_txt << "Emission probabilities" << endl;
    HMMMatrix<double> const& emission_probabilities = hmm.get_emission_probs();
    for (uint i = 0; i < emission_probabilities.get_no_rows(); ++i)
    {
        for (uint j = 0; j < emission_probabilities.get_no_columns(); ++j)
        {
            parameters_txt << emission_probabilities(i, j) << "\t";
        }
        parameters_txt << endl;
    }

    parameters_txt.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc > 1 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")))
    {
        cerr << "This program searches for CpG islands in genome" << endl;
        cerr << "Usage: cat fasta file with genome and pass output to program" << endl;
        cerr << "Example: cat genome.fa | ./HMM_CpG" << endl;
        return 0;
    }

    if (argc > 1 && (!strcmp(argv[1], "--version") || !strcmp(argv[1], "-v")))
    {
        cerr << "1.0.0" << endl;
        return 0;
    }

    string mode;
    if (argc > 1 && (!strcmp(argv[1], "--train") || !strcmp(argv[1], "-t")))
    {
        mode = "train";
    }
    else
    {
        mode = "predict";
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

    Set_HMM_Parameters(initial_probabilities, transition_probabilities, emission_probabilities);

    HMM<double> hmm(initial_probabilities, transition_probabilities, emission_probabilities);

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
            return 1;
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
                    Predict(observed_sequence, hmm);
                }
                else
                {
                    Train(observed_sequence, hmm);
                }

                observed_sequence.clear();
            }
        }
    }

    //save new parameters
    if (mode == "train")
    {
        Save_Parameters(hmm);
    }

    return 0;
}
