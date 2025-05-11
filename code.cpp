//code of the poject of Lynn Limbach & Max Henrotin
#include <ostream>
#include <spoa/spoa.hpp>    //use chatgpt to see how to install spoa library on your computer (use linux)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cmath> 
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

//this method prints an histogram of the number of sequences of each size
void print_sequence_size_histogram(const vector<string>& sequences) {
    vector<int> sizeCount(500, 0); 

    // iterate over each sequences and incrementate the corresponding size
    for (const auto& seq : sequences) {
        int size = seq.length();
        if (size >= 1 && size <= 400) {
            sizeCount[size]++;
        }
    }

    // print the histogram
    for (int i = 1; i < 400; ++i) {
        cout << "Size " << i << " : ";
        
        int numHashes = ceil(sizeCount[i] / 10.0);
        for (int j = 0; j < numHashes; ++j) {
            cout << "#";
        }

        cout << endl;
    }
}

//__________________________________METHODS FOR THE ALGORITHM____________________________

// function that align all the sequences from the argument
vector<string> generate_msa(const vector<string>& sequences) {
  
    auto alignment_engine = spoa::AlignmentEngine::Create(
        //we use kSW because its more robust with such diverse sequences than kNW
        spoa::AlignmentType::kSW, 2, -4, -6, -1);  // linear gaps   (3 : Score de correspondance entre les bases. / -5 : Pénalité pour une ouverture de gap. / -3 : Pénalité pour une extension de gap.)
  
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
        //we use kSW because its more robust with such diverse sequences than kNW
        spoa::AlignmentType::kSW, 2, -4, -6, -1);  // linear gaps   (3 : Score de correspondance entre les bases. / -5 : Pénalité pour une ouverture de gap. / -3 : Pénalité pour une extension de gap.)
  
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


//find the k centroids:

//distance between two string
int levenshtein_distance(const std::string& s1, const std::string& s2) {
    int m = s1.size();
    int n = s2.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));

    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = std::min({
                    dp[i - 1][j] + 1,    // suppression
                    dp[i][j - 1] + 1,    // insertion
                    dp[i - 1][j - 1] + 1 // substitution
                });
            }
        }
    }

    return dp[m][n];
}
    
//k-mean clustering
vector<string> k_centroid(vector<string> sequences, int k){

    vector<string> msa_sequences = generate_msa(sequences);

    //initialize clusters
    std::vector<std::vector<std::string>> clusters(k); 
    
    //initialize centroids
    vector<string> centroids(k);
    
    for(int i=0; i<k;++i){
        centroids[i]=msa_sequences[i];   //select the two first sequences as initial clusters (so its a sort of random pick)
    }


    //10 iterations of k means
    int nbr_step = 10;
    for(int n = 0;n<nbr_step;n++){

        cout<<"K-clustering STEP : "<<n+1<<"/"<<nbr_step<<endl;

        //put each sequence in its corresponding cluster
        for (const auto& seq : msa_sequences) {
            int min_dist = std::numeric_limits<int>::max(); //virtual infinity...
            int assigned_cluster = -1;
        
            for (int i = 0; i < k; ++i) {
                int dist = levenshtein_distance(seq, centroids[i]);     //didnt found any distance method in spoa...
                if (dist < min_dist) {
                    min_dist = dist;
                    assigned_cluster = i;
                }
            }
            clusters[assigned_cluster].push_back(seq);
        }

        //new centroids
        for(int i=0; i<k;++i){
            centroids[i]=find_consensus(clusters[i]); 
        }

        for (auto& cluster : clusters) {
            cluster.clear();
        }
    }

    return centroids;
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

    //sequences size histogram
    print_sequence_size_histogram(sequences);
    
    //with the histogram we can see that we can focus essencially on sequences from size 290 to 305 (with max at 296) and that the rest might be mistakes that will compromise the data
    vector<string> sequences296 = keep_seq_length(sequences, 296, 296);
    vector<string> sequences290_305 = keep_seq_length(sequences, 250, 350);
    cout << endl;
    cout << sequences.size() << " sequences in total" << endl;
    cout << sequences296.size() << " sequences of size 296" << endl;
    cout << sequences290_305.size() << " sequences of a size from 290 to 305" << endl;


    //consensus
    cout << endl << "_______CONSENSUS_OVER_ALL_SEQ__________" << endl;  //for J29B there is one error in the principal consensus (proving what has been estimated with the histogram)
    //string consensus = find_consensus(sequences); //long to calculate
    //cout << consensus << endl << endl;
    cout << endl << "_______CONSENSUS_OVER_SIZE_296_SEQ__________" << endl; //for J29B correct consensus 
    string consensus296 = find_consensus(sequences296);
    cout << consensus296 << endl << endl;
    cout << endl << "_______CONSENSUS_OVER_ALL_SEQ_FROM_290_TO_305__________" << endl;  //for J29B correct consensus
    //string consensus290_305 = find_consensus(sequences290_305);
    //cout << consensus290_305 << endl << endl;

    //alignment (to see visually how it looks)
    vector<string> msa = generate_msa(sequences296);
    //print_sequence_list(msa);

    cout<<endl<<"2-Clustering over all sequences of lenght 296 :"<<endl;
    //vector<string> sequences296_2_centroid = k_centroid(sequences296, 2);
    //print_sequence_list(sequences296_2_centroid);
    
    cout<<endl<<"3-Clustering over all sequences of lenght from 290 to 305 :"<<endl;
    vector<string> sequences290_305_3_centroid = k_centroid(sequences290_305, 2);
    print_sequence_list(sequences290_305_3_centroid);


    return 0;
}