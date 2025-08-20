//code of the poject of Lynn Limbach & Max Henrotin
#include <cstddef>
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
#include <random>
#include <chrono>
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
        cout << seq << endl << endl;
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

string remove_gaps(string sequence){
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '-'), sequence.end());
    return sequence;
}
//to remove all the gap symbols '-' from an array of sequences
vector<string> remove_gaps(vector<string> sequences){
    for(size_t i=0; i<sequences.size();++i){
        sequences[i] = remove_gaps(sequences[i]);
    }
    return sequences;
}

//__________________________________METHODS FOR THE ALGORITHM____________________________

// function that align all the sequences from the argument
vector<string> generate_msa(const vector<string>& sequences) {
  
    auto alignment_engine = spoa::AlignmentEngine::Create(
        //we use kSW because its more robust with such diverse sequences than kNW
        spoa::AlignmentType::kSW,  2, -4, -6, -1);  // linear gaps   (3 : Score de correspondance entre les bases. / -5 : Pénalité pour une ouverture de gap. / -3 : Pénalité pour une extension de gap.)
  
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
        spoa::AlignmentType::kSW,  2, -4, -6, -1);  // linear gaps   (3 : Score de correspondance entre les bases. / -5 : Pénalité pour une ouverture de gap. / -3 : Pénalité pour une extension de gap.)
  
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



//distance between two string   --> does not seem to perform better at all and VERY SLOW !!
int levenshtein_distance(const string& s1, const string& s2) {
    int m = s1.size();
    int n = s2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));

    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min({
                    dp[i - 1][j] + 1,    // suppression
                    dp[i][j - 1] + 1,    // insertion
                    dp[i - 1][j - 1] + 1 // substitution
                });
            }
        }
    }

    return dp[m][n];
}

//compare each char one by one and increment if there is a difference
int hamming_distance(const string& s1, const string& s2) {
    if (s1.size() != s2.size()) {
        cerr << "must be same size" << std::endl;
        return -1;
    }

    int distance = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) {
            ++distance;
        }
    }
    return distance;
}


/*
//spoa distance between two string --> ??? I DONT KNOW HOW TO DO THAT
int spoa_distance(const string& s1, const string& s2){
    // Crée un aligner SPOA (alignement global ici)
    auto alignment_engine = spoa::AlignmentEngine::Create(
        //we use kSW because its more robust with such diverse sequences than kNW
        spoa::AlignmentType::kSW, 2, -4, -6, -1);

    // Crée un graphe vide
    spoa::Graph graph{};

    // Ajoute la première séquence dans le graphe
    auto alignment1 = alignment_engine->Align(s1, graph);
    graph.AddAlignment(alignment1, s1);

    // Aligne la deuxième séquence sur le graphe
    auto alignment2 = alignment_engine->Align(s2, graph);
    graph.AddAlignment(alignment2, s2);

    // Score d’alignement de la 2e séquence par rapport au graphe
    int score = alignment_engine->Score(s2, graph); //????

    //return -score;  //distance is inverse of score
}
*/
   
//--> NOT USED BECAUSE DOES NOT SEEM TO CHANGE MUCH THE OUTPUT
/*
//initialise optimaly the centroids (with kmean++) --> written by chatgpt
vector<string> kmeanspp_init_centroids(const vector<string>& sequences, int k) {
    if (sequences.empty() || k <= 0 || k > sequences.size()) {
        throw invalid_argument("Invalid input: sequences empty or k out of bounds");
    }

    vector<string> centroids;
    vector<double> distances(sequences.size(), numeric_limits<double>::max());

    random_device rd;
    mt19937 gen(rd());

    // 1. Choisir un centroïde initial au hasard
    uniform_int_distribution<> dis(0, sequences.size() - 1);
    centroids.push_back(sequences[dis(gen)]);

    while (centroids.size() < static_cast<size_t>(k)) {
        double total = 0.0;

        // 2. Calculer la distance min² de chaque séquence aux centroïdes existants
        for (size_t i = 0; i < sequences.size(); ++i) {
            for (const auto& centroid : centroids) {
                double d = levenshtein_distance(sequences[i], centroid);
                distances[i] = min(distances[i], d * d);
            }
            total += distances[i];
        }

        // 3. Roulette wheel : sélectionner le prochain centroïde avec proba ~ distance²
        uniform_real_distribution<> rand_weight(0, total);
        double r = rand_weight(gen);

        double cumulative = 0.0;
        size_t next_idx = 0;
        for (; next_idx < sequences.size(); ++next_idx) {
            cumulative += distances[next_idx];
            if (cumulative >= r)
                break;
        }

        centroids.push_back(sequences[next_idx]);
    }

    return centroids;
}
*/

//k-mean clustering
vector<string> k_centroid(vector<string> sequences, int k, int nbr_step_max){

    vector<string> msa_sequences = generate_msa(sequences);

    //initialize clusters
    vector<vector<string>> clusters(k); 
    
    //initialize centroids  --> MAYBE TRY K-MEAN++ for a better initialization and thus prevent empty clusters
    vector<string> centroids(k);
    for(int i=0; i<k;++i){
        centroids[i]=msa_sequences[i];   //select the two first sequences as initial clusters (so its a sort of random pick)
    }
    //centroids = kmeanspp_init_centroids(msa_sequences, k);


    //nbr_step iterations of k means
    for(int n = 0;n<nbr_step_max;n++){
        cout<<"K-clustering STEP : "<<n+1<<"/"<<nbr_step_max<<endl;

        //put each sequence in its corresponding cluster
        size_t nbr_seq = msa_sequences.size();
        for (size_t i = 0; i < nbr_seq; ++i) {
            //cout<<"K-clustering STEP : "<<n+1<<"/"<<nbr_step<<"      ->       sequence : "<< i+1<<"/"<<nbr_seq<<endl;
            
            const string& seq = msa_sequences[i];
            int min_dist = std::numeric_limits<int>::max(); //virtual infinity...
            int assigned_cluster = -1;
            for (int j = 0; j < k; ++j) {
                int dist = hamming_distance(seq, centroids[j]);     //didnt found any distance method in spoa...
                if (dist < min_dist) {
                    min_dist = dist;
                    assigned_cluster = j;
                }
            }
            clusters[assigned_cluster].push_back(seq);
        }

        //new centroids
        bool converged = true;
        for(int i=0; i<k;++i){
            string new_centroid = find_consensus(remove_gaps(clusters[i])); 
            
            //align the new centroid with all the sequences
            //cout<<new_centroid.size()<<endl;
            msa_sequences.push_back(new_centroid);
            msa_sequences = generate_msa(remove_gaps(msa_sequences));
            new_centroid = msa_sequences.back(); 
            msa_sequences.pop_back();
            //cout<<new_centroid.size()<<endl<<endl;

            //check convergence
            if(remove_gaps(centroids[i]) != remove_gaps(new_centroid)){
                centroids[i]=new_centroid;
                converged = false;
            }
        }
        //IF THERE IS NO CHANGE IN CENTROID THEN STOP LOOPING
        if(converged){
            cout<<"K-clustering has converged !!!"<<endl<<endl;
            for(int i = 0; i<k; ++i){
                cout<<"cluster no "<<i+1<<" : "<<clusters[i].size()<<endl;
            }
            return centroids;
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
    auto startTime = std::chrono::high_resolution_clock::now();

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
    vector<string> sequences290_305 = keep_seq_length(sequences, 290, 305);
    cout << endl;
    cout << sequences.size() << " sequences in total" << endl;
    cout << sequences296.size() << " sequences of size 296" << endl;
    cout << sequences290_305.size() << " sequences of a size from 290 to 305" << endl;


    //consensus
    cout << endl << "_______CONSENSUS_OVER_ALL_SEQ__________" << endl;  //for J29B there is one error in the principal consensus (proving what has been estimated with the histogram)
    //string consensus = find_consensus(sequences); //long to calculate
    //cout << consensus << endl << endl;
    cout << endl << "_______CONSENSUS_OVER_SIZE_296_SEQ__________" << endl; //for J29B correct consensus 
    //string consensus296 = find_consensus(sequences296);
    //cout << consensus296 << endl << endl;
    cout << endl << "_______CONSENSUS_OVER_ALL_SEQ_FROM_290_TO_305__________" << endl;  //for J29B correct consensus
    //string consensus290_305 = find_consensus(sequences290_305);
    //cout << consensus290_305 << endl << endl;

    //alignment (to see visually how it looks)
    //vector<string> msa = generate_msa(sequences296);
    //print_sequence_list(msa);

    cout<<endl<<"2-Clustering over all sequences of lenght 296 :"<<endl;
    //vector<string> sequences296_2_centroid = k_centroid(sequences296, 2, 10); //J29B : 1st variant ok, 2nd variant ok / J30B : 1st variant ok, 2nd variant 1err (with hamming distance)
    //print_sequence_list(remove_gaps(sequences296_2_centroid));

    cout<<endl<<"3-Clustering over all sequences of lenght 296 :"<<endl;
    vector<string> sequences296_3_centroid = k_centroid(sequences296, 8, 10); //J29B : 1st variant ok, 2nd variants ok, 3rd variant nope / J30B : 1st variant ok, 2nd variant 1 err, 3rd variant nope (with hamming distance)
    print_sequence_list(remove_gaps(sequences296_3_centroid));

    
    cout<<endl<<"3-Clustering over all sequences of lenght from 290 to 305 :"<<endl;
    //vector<string> sequences290_305_3_centroid = k_centroid(sequences290_305, 3, 10);   //J29B : 1st variant ok, 2nd variant 1err, 3rd variant nope / J30B : 1st variant ok, 2nd variant 1err, 3rd variant nope (with hamming_distance)
    //print_sequence_list(remove_gaps(sequences290_305_3_centroid));


    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;
    cout<<"TOTAL TIME : "<<duration.count()<<" seconds"<<endl;
    return 0;
}