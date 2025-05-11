//code of the poject of Lynn Limbach & Max Henrotin
#include <spoa/spoa.hpp>    //use chatgpt to see how to install spoa library on your computer (use linux)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdlib>
using namespace std;
using namespace spoa;  

//____________________HELPERS FOR SEQUENCES EXTRACTION___________________________
//this method extract in an array all the ACTG sequences from the .fastq file
vector<string> extract_sequences(const string& path) {
    ifstream file(path);
    if (!file.is_open()) {
        cerr << "Error: The file " << path << " was not found or could not be opened." << endl;
        exit(1);
    }

    vector<string> sequences;
    string line;
    int line_number = 0;
    while (getline(file, line)) {
        line_number++;
        if (line_number % 4 == 2) { //we assume that the fastq file is in the right format
            sequences.push_back(line);
        }
    }

    return sequences;
}

//method that filter only the sequences within a given lenght region
vector<string> keep_seq_length(const vector<string>& sequences, int minLength, int maxLength) {
    if (minLength > maxLength) {
        cerr << "minLength > maxLength" << endl;
        exit(1);
    }

    vector<string> filtered;
    for (const string& seq : sequences) {
        if (seq.length() >= minLength && seq.length() <= maxLength) {
            filtered.push_back(seq);
        }
    }
    return filtered;
}

//print a sequence list
void print_sequence_list(const vector<string>& sequences) {
    for (const string& seq : sequences) {
        cout << "- " << seq << endl << endl;
    }
}

//__________________________________METHODS FOR THE ALGORITHM____________________________
//k-medoid algorithm --> using SPOA

// function that align all the sequences from the argument
vector<string> generate_msa(const vector<string>& sequences) {
  
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps   (3 : Score de correspondance entre les bases. / -5 : Pénalité pour une ouverture de gap. / -3 : Pénalité pour une extension de gap.)
  
    spoa::Graph graph{};    //to store alignment sequences
  
    //iterate over all sequences and put them in the graph at the right place regarding the alignment_engine
    for (const auto& it : sequences) {
      auto alignment = alignment_engine->Align(it, graph);
      graph.AddAlignment(alignment, it);
    }
  
    //generate a consensus
    string consensus = graph.GenerateConsensus();
  
    //generate the multiple sequence alignment
    vector<string> msa = graph.GenerateMultipleSequenceAlignment();

    return msa;
}

//find one consensus string (maybe can be mixed together with the msa method)
string find_consensus(const vector<string>& sequences) {

    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps   (3 : Score de correspondance entre les bases. / -5 : Pénalité pour une ouverture de gap. / -3 : Pénalité pour une extension de gap.)
  
    spoa::Graph graph{};    //to store alignment sequences
  
    //iterate over all sequences and put them in the graph at the right place regarding the alignment_engine
    for (const auto& it : sequences) {
      auto alignment = alignment_engine->Align(it, graph);
      graph.AddAlignment(alignment, it);
    }
  
    //generate a consensus
    string consensus = graph.GenerateConsensus();

    return consensus;
}

//_________________________EXECUTION__________________________________

//to compile : g++ -o code code.cpp -lspoa
//OR ?? : g++ -o code code.cpp -I~/spoa/spoa/include -L~/spoa/build/lib -lspoa
//to execute : ./code fastq/J29_B_CE_IonXpress_005.fastq    //example of execution for this file



int main(int argc, char* argv[]) {
    if (argc != 2) {    //fastq/J29_B_CE_IonXpress_005.fastq  is one example of argument
        cerr << "Please provide the path of the FASTQ file." << endl;
        return 1;
    }

    //extract sequences
    string file_path = argv[1];
    vector<string> sequences = extract_sequences(file_path);
    vector<string> sequences296 = keep_seq_length(sequences, 296, 296);
    vector<string> sequences250_350 = keep_seq_length(sequences, 250, 350);
    cout << sequences.size() << "sequences in total" << endl;
    cout << sequences296.size() << "sequences of size 296" << endl;
    cout << sequences250_350.size() << "sequences of a size from 250 to 350" << endl;


    //alignment
    vector<string> msa = generate_msa(sequences296);
    //print_sequence_list(msa);


    //consensus
    cout << endl << "_______CONSENSUS_OVER_SIZE_296_SEQ__________" << endl; //for J29B correct consensus 
    string consensus296 = find_consensus(sequences296);
    cout << consensus296 << endl << endl;
    cout << endl << "_______CONSENSUS_OVER_ALL_SEQ__________" << endl;  //for J29B there is 2 errors in the principal consensus
    string consensus = find_consensus(sequences);
    cout << consensus << endl << endl;
    cout << endl << "_______CONSENSUS_OVER_ALL_SEQ_FROM_250_TO_350__________" << endl;  //for J29B there is no errors in the principal consensus
    string consensus250_350 = find_consensus(sequences250_350);
    cout << consensus250_350 << endl << endl;


    return 0;
}