#include "precompiled.h"
#include "HMMlib/hmm_table.hpp"
#include "HMMlib/hmm_vector.hpp"
#include "HMMlib/hmm.hpp"

//    vector<char> state_mapping(8);
//    state_mapping[0] = 'A';
//    state_mapping[1] = 'C';
//    state_mapping[2] = 'G';
//    state_mapping[3] = 'T';
//    state_mapping[4] = 'a';
//    state_mapping[5] = 'c';
//    state_mapping[6] = 'g';
//    state_mapping[7] = 't';

using hmmlib::HMM;
using hmmlib::HMMMatrix;
using hmmlib::HMMVector;

enum Mode
{
    predict,
    train
};

bool Check_Sequence(sequence& seq, uint start, uint end)
{
    uint cg_counter = 0;
    uint c_counter = 0;
    uint g_counter = 0;
    uint cpg_counter = 0;
    uint i = 0;
    uint length = end - start;

    if (length < 200)
    {
        return false;
    }

    for (i = start; i < end; ++i)
    {
        if (seq[i] == 1)
        {
            c_counter += 1;
            cg_counter += 1;
        }

        if (seq[i] == 2)
        {
            g_counter += 1;
            cg_counter += 1;
        }

        if (seq[i] == 1 && seq[i + 1] == 2)
        {
            cpg_counter += 1;
        }
    }

    if (seq[i] == 1)
    {
        c_counter += 1;
        cg_counter += 1;
    }

    if (seq[i] == 2)
    {
        g_counter += 1;
        cg_counter += 1;
    }

    if (cg_counter/(double)length < 0.5)
    {
        return false;
    }

    if (cpg_counter/(double)(c_counter*g_counter)*length < 0.6)
    {
        return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Print_Cpg(sequence& hidden_sequence, uint shift, string& chromosome_name)
{
    uint start_cpg = 0;
    bool is_open_cpg = false;
    ull i;
    for (i = 0; i < hidden_sequence.size(); ++i)
    {
        if (hidden_sequence[i] < 4)
        {
            if (!is_open_cpg)
            {
                start_cpg = i + shift + 1;
                is_open_cpg = true;
            }
        }
        else
        {
            if (is_open_cpg)
            {
                uint end_cpg = i + shift + 1;
                if (Check_Sequence(hidden_sequence, start_cpg - shift - 1, end_cpg - shift - 1))
                {
                    cout << chromosome_name << "\t" << start_cpg << "\t" << end_cpg << endl;
                }
                is_open_cpg = false;
            }
        }
    }

    if (is_open_cpg)
    {
        uint end_cpg = i + shift + 1;
        if (Check_Sequence(hidden_sequence, start_cpg, end_cpg))
        {
            cout << chromosome_name << "\t" << start_cpg << "\t" << end_cpg << endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

shared_ptr< HMM<double> > Make_Model(uint number_of_states, uint alphabet_size, ifstream& parameters)
{
    shared_ptr< HMMVector<double> > initial_probabilities(new HMMVector<double>(number_of_states));
    shared_ptr< HMMMatrix<double> > transition_probabilities(new HMMMatrix<double>(number_of_states, number_of_states));
    shared_ptr< HMMMatrix<double> > emission_probabilities(new HMMMatrix<double>(alphabet_size, number_of_states));

    Set_HMM_Parameters(initial_probabilities, transition_probabilities, emission_probabilities, parameters);

    HMM<double>* hmm_ptr = new HMM<double>(initial_probabilities, transition_probabilities, emission_probabilities);
    shared_ptr< HMM<double> > hmm_sptr(hmm_ptr);

    return hmm_sptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

shared_ptr< HMM<double> > Merge_Models(shared_ptr < HMM<double> > cpg, shared_ptr< HMM<double> > non_cpg, uint average_cpg_length, uint average_non_cpg_length)
{
    if (cpg->get_no_states() != non_cpg->get_no_states())
    {
        throw("Models states number must be same");
    }

    if (cpg->get_alphabet_size() != non_cpg->get_alphabet_size())
    {
        throw("Models alphabet size must be same");
    }

    double leave_cpg_probability = 1/(double)average_cpg_length;
    double stay_cpg_probability = 1 - leave_cpg_probability;
    double leave_non_cpg_probability = 1/(double)average_non_cpg_length;
    double stay_non_cpg_probability = 1 - leave_non_cpg_probability;

    //initial probabilities
    shared_ptr< HMMVector<double> > initial_probabilities(new HMMVector<double>(cpg->get_no_states()*2));
    uint i;
    for (i = 0; i < cpg->get_no_states(); ++i)
    {
        (*initial_probabilities)(i) = cpg->get_initial_probs()(i)/2;
    }
    for (uint j = 0; j < non_cpg->get_no_states(); ++j)
    {
        (*initial_probabilities)(j + i) = non_cpg->get_initial_probs()(j)/2;
    }

    //transition probabilities
    shared_ptr< HMMMatrix<double> > transition_probabilities(new HMMMatrix<double>(cpg->get_no_states()*2, cpg->get_no_states()*2));
    for (uint i = 0; i < transition_probabilities->get_no_rows(); ++i)
    {
        for (uint j = 0; j < transition_probabilities->get_no_columns(); ++j)
        {
            if (i < cpg->get_no_states() && j < cpg->get_no_states())
            {
                (*transition_probabilities)(i, j) = cpg->get_trans_probs()(i, j)*stay_cpg_probability;
            }
            if (i < cpg->get_no_states() && j >= cpg->get_no_states())
            {
                (*transition_probabilities)(i, j) = cpg->get_trans_probs()(i, j - cpg->get_no_states())*leave_cpg_probability;
            }
            if (i >= cpg->get_no_states() && j < cpg->get_no_states())
            {
                (*transition_probabilities)(i, j) = non_cpg->get_trans_probs()(i - cpg->get_no_states(), j)*leave_non_cpg_probability;
            }
            if (i >= cpg->get_no_states() && j >= cpg->get_no_states())
            {
                (*transition_probabilities)(i, j) =
                        non_cpg->get_trans_probs()(i - cpg->get_no_states(), j - cpg->get_no_states())*stay_non_cpg_probability;
            }
        }
    }

    //emission probabilities
    shared_ptr< HMMMatrix<double> > emission_probabilities(new HMMMatrix<double>(cpg->get_alphabet_size(), cpg->get_no_states()*2));
    for (uint i = 0; i < cpg->get_alphabet_size(); ++i)
    {
        for (uint j = 0; j < cpg->get_no_states()*2; ++j)
        {
            if (j < cpg->get_no_states())
            {
                (*emission_probabilities)(i, j) = cpg->get_emission_probs()(i, j);
            }
            else
            {
                (*emission_probabilities)(i, j) = non_cpg->get_emission_probs()(i, j - cpg->get_no_states());
            }
        }
    }

    HMM<double>* hmm = new HMM<double>(initial_probabilities, transition_probabilities, emission_probabilities);
    hmm->Save_Parameters();

    return shared_ptr< HMM<double> >(hmm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint Process_Observed(string& observed, shared_ptr< HMM<double> > hmm, Mode mode, string chromosome_name = string(""))
{
    vector<int> nucleotides_mapping(127, -1);
    nucleotides_mapping['a'] = 0;
    nucleotides_mapping['c'] = 1;
    nucleotides_mapping['g'] = 2;
    nucleotides_mapping['t'] = 3;

    uint shift = 0;
    //HMM<double> hmm = *hmm_sptr;
    sequence observed_sequence;
    observed_sequence.reserve(observed.length());
    sequence hidden_sequence;
    hidden_sequence.reserve(observed.length());

    transform(observed.begin(), observed.end(), observed.begin(), ::tolower);

    uint sum = 0;
    uint count = 0;
    for (ull i = 0; i < observed.length(); ++i)
    {
        if (observed[i] != 'a' && observed[i] != 'c' && observed[i] != 'g' && observed[i] != 't' && observed[i] != 'n')
        {
            cerr << "Unexpected character. Possible characters are a, c, g, t, n (in any register)" << endl;
            return 0;
        }

        if (observed[i] != 'n')
        {
            if (observed_sequence.size() == 0)
            {
                shift = i;
            }

            observed_sequence.push_back(nucleotides_mapping[observed[i]]);
        }
        else
        {
            if (observed_sequence.size() > 0)
            {
                ++count;
                sum += observed_sequence.size();
                if (mode == predict)
                {
                    hmm->Predict(observed_sequence, hidden_sequence);
                    Print_Cpg(hidden_sequence, shift, chromosome_name);
                }
                else
                {
                    hmm->Train(observed_sequence);
                }

                observed_sequence.clear();
                hidden_sequence.clear();
            }
        }
    }

    if (observed_sequence.size() > 0)
    {
        ++count;
        sum += observed_sequence.size();
        if (mode == predict)
        {
            hmm->Predict(observed_sequence, hidden_sequence);
            Print_Cpg(hidden_sequence, shift, chromosome_name);
        }
        else
        {
            hmm->Train(observed_sequence);
        }
    }

    return sum/count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint Train(ifstream& cpg_file, shared_ptr< HMM<double> > hmm_sptr)
{
    cerr << "Start training" << endl;
    string observed;

    string line;
    while (!cpg_file.eof())
    {
        cpg_file >> line;
        observed += line;
    }

    return Process_Observed(observed, hmm_sptr, train);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Predict(ifstream& parameters_file, string& chromosome_name)
{
    cerr << "Start prediction" << endl;

    shared_ptr< HMM<double> > hmm_sptr;
    try
    {
        hmm_sptr = Make_Model(8, 4, parameters_file);
    }
    catch (char const* msg)
    {
        cerr << msg << endl;
        return;
    }

    string observed;
    string line;
    while (!cin.eof())
    {
        cin >> line;
        observed += line;
    }

    Process_Observed(observed, hmm_sptr, predict, chromosome_name);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{    
    if (argc < 2)
    {
        cerr << "You shoud specify args. Use --help option if you don't want which.\n Your current directory should be directory where executable is stored" << endl;
        return 1;
    }

    Mode mode = predict;

    bool parameters_expectation = false;
    bool cpg_expectation = false;
    bool non_cpg_expectation = false;
    bool chromosome_name_expectation = false;
    ifstream parameters_file;
    ifstream cpg_train;
    ifstream non_cpg_train;
    ifstream cpg_probabilities;
    ifstream non_cpg_probabilities;
    bool cpg_file = false;
    bool non_cpg_file = false;
    uint average_cpg_length = 765;
    uint average_non_cpg_length = 100415;
    bool average_cpg_expectation = false;
    bool average_non_cpg_expectation = false;
    bool user_defined_average_cpg = false;
    bool user_defined_average_non_cpg = false;
    string chromosome_name = "undefined";

    //command line arguments handling
    for (int arg = 1; arg < argc; ++arg)
    {
        if (parameters_expectation)
        {
            parameters_expectation = false;
            parameters_file.open(argv[arg]);
        }

        if (cpg_expectation)
        {
            cpg_expectation = false;
            cpg_train.open(argv[arg]);

            if (!cpg_train.good())
            {
                cerr << "Bad cpg training file" << endl;
                return 1;
            }
            cpg_file = true;
        }

        if (non_cpg_expectation)
        {
            non_cpg_expectation = false;
            non_cpg_train.open(argv[arg]);

            if (!non_cpg_train.good())
            {
                cerr << "Bad non-cpg training file" << endl;
                return 1;
            }
            non_cpg_file = true;
        }

        if (average_cpg_expectation)
        {
            average_cpg_expectation = false;
            string number(argv[arg]);

            if (find_if(number.begin(), number.end(), [](char c) {return !isdigit(c);}) != number.end())
            {
                cerr << "Positive integer expected after -acl" << endl;
                return 4;
            }

            average_cpg_length = atoi(argv[arg]);
            user_defined_average_cpg = true;
        }

        if (average_non_cpg_expectation)
        {
            average_non_cpg_expectation = false;
            string number(argv[arg]);

            if (find_if(number.begin(), number.end(), [](char c) {return !isdigit(c);}) != number.end())
            {
                cerr << "Positive integer expected after -acl" << endl;
                return 4;
            }

            average_non_cpg_length = atoi(argv[arg]);
            user_defined_average_non_cpg = true;
        }

        if (chromosome_name_expectation)
        {
            chromosome_name_expectation = false;
            chromosome_name = string(argv[arg]);
        }

        if (!strcmp(argv[arg], "--help") || !strcmp(argv[arg], "-h"))
        {
            ifstream help("..\\README.md");
            string text;
            text.assign(istreambuf_iterator<char>(help), istreambuf_iterator<char>());
            cerr << text << endl;
            help.close();
            return 0;
        }

        if (!strcmp(argv[arg], "--version") || !strcmp(argv[arg], "-v"))
        {
            cerr << "1.0.0" << endl;
            return 0;
        }

        if (!strcmp(argv[arg], "--cpg"))
        {
            mode = train;
            cpg_expectation = true;
        }

        if (!strcmp(argv[arg], "--non-cpg"))
        {
            mode = train;
            non_cpg_expectation = true;
        }

        if (!strcmp(argv[arg], "--parameters") || !strcmp(argv[arg], "-p"))
        {
            parameters_expectation = true;
        }

        if (!strcmp(argv[arg], "--average-cpg-length") || !strcmp(argv[arg], "-acl"))
        {
            average_cpg_expectation = true;
        }

        if (!strcmp(argv[arg], "--average-non-cpg-length") || !strcmp(argv[arg], "-ancl"))
        {
            average_non_cpg_expectation = true;
        }

        if (!strcmp(argv[arg], "--chromosome-name") || !strcmp(argv[arg], "-cn"))
        {
            chromosome_name_expectation = true;
        }
    }

    if (!parameters_file.is_open() && mode == predict)
    {
        parameters_file.open("parameters");
    }

    if (!parameters_file.good())
    {
        cerr << "Can't read file" << endl;
        return 1;
    }

    if (mode == train)
    {
        cpg_probabilities.open("cpg_probabilities");
        if (!cpg_probabilities.good())
        {
            cerr << "Bad cpg_probabilities file" << endl;
            return 1;
        }

        non_cpg_probabilities.open("non_cpg_probabilities");
        if (!non_cpg_probabilities.good())
        {
            cerr << "Bad non_cpg_probabilities file" << endl;
            return 1;
        }

        shared_ptr< HMM<double> > cpg_hmm;
        try
        {
            cpg_hmm = Make_Model(4, 4, cpg_probabilities);
        }
        catch (char const* msg)
        {
            cerr << msg << endl;
            return 5;
        }
        cpg_probabilities.close();


        shared_ptr< HMM<double> > non_cpg_hmm;
        try
        {
            non_cpg_hmm = Make_Model(4, 4, non_cpg_probabilities);
        }
        catch (char const* msg)
        {
            cerr << msg << endl;
            return 5;
        }
        non_cpg_probabilities.close();

        if (cpg_file)
        {
            uint average_train_cpg = Train(cpg_train, cpg_hmm);
            cerr << "average training cpg length " << average_train_cpg << endl;

            cpg_train.close();

            if (!user_defined_average_cpg)
            {
                average_cpg_length = average_train_cpg;
            }
        }

        if (non_cpg_file)
        {
            uint average_train_non_cpg = Train(non_cpg_train, non_cpg_hmm);
            cerr << "average training non-cpg length " << average_train_non_cpg << endl;

            non_cpg_train.close();

            if (!user_defined_average_non_cpg)
            {
                average_non_cpg_length = average_train_non_cpg;
            }
        }

        Merge_Models(cpg_hmm, non_cpg_hmm, average_cpg_length, average_non_cpg_length);
    }
    else
    {
        Predict(parameters_file, chromosome_name);
        parameters_file.close();
    }

    return 0;
}
