//code of the poject of Lynn Limbach & Max Henrotin
#include <spoa/spoa.hpp>    //use chatgpt to see how to install spoa library on your computer (use linux)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <algorithm>
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

// Fonction pour trouver les k medoids
vector<string> find_k_medoids(const vector<string>& seq, int k) {
    std::vector<std::string> sequences = {
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };
  
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
  
    spoa::Graph graph{};
  
    for (const auto& it : sequences) {
      auto alignment = alignment_engine->Align(it, graph);
      graph.AddAlignment(alignment, it);
    }
  
    auto consensus = graph.GenerateConsensus();
  
    std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
              << consensus << std::endl;
  
    auto msa = graph.GenerateMultipleSequenceAlignment();
  
    for (const auto& it : msa) {
      std::cerr << it << std::endl;
    }
    return seq;
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

    string file_path = argv[1];
    vector<string> sequences = extract_sequences(file_path);
    vector<string> sequences296 = keep_seq_length(sequences, 296, 296);
    print_sequence_list(sequences296);
    cout << sequences.size() << endl;
    cout << sequences296.size() << endl;

    return 0;
}