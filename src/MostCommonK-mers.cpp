//============================================================================
// Name        : MostCommonK-mers.cpp
// Author      : Taylan Dogan
// Version     :
// Copyright   : The MIT License (MIT)
// Description : Given a .fastq file, finds 25 most common KMER_SIZE-length DNA subsequences.
//============================================================================

#include <iostream>
#include <algorithm>
#include <ctime>

#include <string>
#include <vector>
#include <unordered_map>

#include <fstream>
#include <sstream>

// Global variables - default settings
int N = 25;
int KMER_SIZE = 30;
int MAP_SIZE_LIMIT = 1000;
float FAIRNESS_CONST = 1.25;
bool FAIR_THRESHOLD_MODE = true;
unsigned long LINE_COUNT = 0;

using namespace std;

// K-mer struct contains a DNA-subsequence and its hash value.
struct kmer {
    string pattern;
    int hash_value;

    // Equals operator is overloaded for kmer struct.
    // It first checks hash values. If they match, proceed with string comparison.
    bool operator== (const kmer& k) const {
        // If stored hashes do not match, return false
        if(k.hash_value != this->hash_value)
            return false;
        // If they do, compare strings
        else {
            if(k.pattern.compare(this->pattern) == 0)
                return true;
            else
                return false;
        }
    }
};

namespace std {
    template <>
    struct hash<kmer> {
        std::size_t operator()(const kmer& k) const {
          return k.hash_value;
        }
    };
}


// Calculates 2 hash values, returns them as a pair.
// - actual hash is the hash used in map entry (for comparison)
// - raw hash is used to calculate next hash value faster (rolling hash)
pair<int, int> calculate_hash(string subseq, int prev_hash, char prev_subseq_first_char) {
    int actual_hash = 0, raw_hash = 0;

    // If you have the previous hash value (for rolling hash)
    if(prev_hash != -1) {
        char last_char = subseq.back();
        raw_hash = prev_hash + int(last_char) - int(prev_subseq_first_char);
    }

    // If this is the first hash you are calculating, in a given DNA sequence line
    else {
        for(char c : subseq) {
            raw_hash += int(c);
        }
    }

    actual_hash = raw_hash + int(subseq.front()) + int(subseq.back()) + int(subseq[KMER_SIZE / 2]);
    return std::make_pair(raw_hash, actual_hash);
}


// Given a DNA-sequence, extracts and stores subsequences
int process_seq(string dna_seq, unordered_map<kmer, int>& kmer_map, int most_common) {
    int most_common_count = most_common;
    int prev_hash = -1;
    char prev_first_char = '\0';
    int seq_length = dna_seq.length();

    // Given DNA sequence is shorter than 30 characters
    if (seq_length < KMER_SIZE)
        return -1;

    // Extract 30 chars-length substrings
    // Calculate their hash using rolling hash
    // Add them to the map
    for(int i = 0; i + KMER_SIZE <= seq_length ; i++) {
        string subseq = dna_seq.substr(i, KMER_SIZE);

        // Calculate hash returns (raw_hash, actual_hash)
        pair<int, int> hash_pair = calculate_hash(subseq, prev_hash, prev_first_char);
        kmer a_kmer = { subseq, hash_pair.second };

        // Find or create K-mer, increment its value, get its value, compare with the max count
        // The uncommented line below does the same,
        // However this requires 2 map searches instead of 1:
        //		kmer_map[a_kmer]++;
        //		int hit_count = kmer_map[a_kmer];

        int hit_count = ++kmer_map[a_kmer];

        // Update the count of most common k-mer
        if(hit_count > most_common_count)
            most_common_count = hit_count;

        // Update previous hash, and first char for rolling hash
        prev_hash = hash_pair.first;
        prev_first_char = subseq[0];
    }

    return most_common_count;
}


// Attempts to clean the k-mer map, using an adaptive threshold
// If cleaning is not possible, it doubles the MAP_SIZE_LIMIT
void clean_map(unordered_map<kmer, int>& kmer_map, int most_common_count) {
    /* --Explanation concerning threshold modes:
     * 1. Fair thresholding:
     * 		Compares 2 measures, fairness_factor and x_factor.
     * 		fairness_factor: ratio of LINE_COUNT / count of most common k-mer
     * 		x_factor: MAP_SIZE_LIMIT / count of a k-mer
     *
     * 		The extent of fairness can be increased by increasing the FAIRNESS_CONST.
     * 		But it is not tested thoroughly and it should be in some safe range for now.
     *
     * 		Comparing these 2 measures actually make sense, because it basically
     * 		gives a chance to the k-mers with low counts by comparing their
     * 		frequency to the frequency of the most common k-mer.
     *
     * 		However, since it is too costly to calculate x_factor for every map entry,
     * 		the threshold mode is changed after some point. (when MAP_SIZE_LIMIT surpasses 15000)
     *
     * 2. Thresholding according to the count of most common k-mer:
     * 		Calculates a threshold, which is, (most_common_count / 20) + 1
     * 		This is a much more harsh way of eliminating k-mers and may lead
     * 		to erroneous results. It is a quick and dirty solution to handle
     * 		complexity.
     *
    */

    int map_size_before = kmer_map.size();

    // Check the threshold mode and proceed accordingly
    if(FAIR_THRESHOLD_MODE) {
        int fairness_factor = (LINE_COUNT / most_common_count) * FAIRNESS_CONST;

        for(auto it = kmer_map.cbegin(); it != kmer_map.cend();) {
            int x_factor = (MAP_SIZE_LIMIT / it->second);

            if(fairness_factor < x_factor) {
                kmer_map.erase(it++);
            }

            else {
                ++it;
            }
        }
    }

    else {
        int threshold = (most_common_count / 20) + 1;

        for(auto it = kmer_map.cbegin(); it != kmer_map.cend();) {
            if(it->second <= threshold) {
                kmer_map.erase(it++);
            }

            else {
                ++it;
            }
        }
    }

    int map_size_after = kmer_map.size();

    // Increase the MAP_SIZE_LIMIT, because we could not eliminate at least 10% of the map in this iteration
    if(!((map_size_before - map_size_after) > (map_size_before / 10))) {
        MAP_SIZE_LIMIT = MAP_SIZE_LIMIT + 1000;

        if(MAP_SIZE_LIMIT >= 15000) {
            FAIR_THRESHOLD_MODE = false;
        }
        //std::cout << "Map size limit: " << MAP_SIZE_LIMIT << endl; //(DEBUG)
    }
}


void write_map_to_file(unordered_map<kmer, int>& kmer_map, string filename) {
    // Copy the map into a vector
    std::vector< pair<string, int> > subseq_vector;
    for(auto i : kmer_map) {
        subseq_vector.push_back(std::make_pair(i.first.pattern, i.second));
    }

    // Sort the vector by the hit count
    // The lambda expression used here, sorts the vector in the decreasing order
    sort(subseq_vector.begin(), subseq_vector.end(), [](const pair<string, int>& seq1, const pair<string, int>& seq2) {
            return seq1.second > seq2.second; }
    );

    // Write vector to the file
    ofstream result_file;
    result_file.open(filename);

    for(auto i : subseq_vector) {
        string vector_entry = i.first + " | " + to_string(i.second);
        result_file << vector_entry << endl;
    }

    result_file.close();
    std::cout << "The results are written to file: " << filename << endl;
}


void extract_top_N_kmers(unordered_map<kmer, int>& kmer_map, int N) {
    // Copy the map into a vector
    std::vector< pair<string, int> > subseq_vector;
    for(auto i : kmer_map) {
        subseq_vector.push_back(std::make_pair(i.first.pattern, i.second));
    }

    // Sort the vector by the hit count
    // The lambda expression used here, sorts the vector in the decreasing order
    sort(subseq_vector.begin(), subseq_vector.end(), [](const pair<string, int>& seq1, const pair<string, int>& seq2) {
            return seq1.second > seq2.second; }
    );

    std::cout << endl;
    std::cout << std::min(N, (int)subseq_vector.size()) << " most frequent " << KMER_SIZE << "-mers are:" << endl;


    for(int i = 0; i < std::min(N, (int)subseq_vector.size()); i++) {
        string vector_entry = subseq_vector.at(i).first + " | " + to_string(subseq_vector.at(i).second);
        std::cout << vector_entry << endl;
    }
}


int main(int argc, char* argv[]) {
    string input_filename = "";
    string output_filename = "output.txt";

    // Parsing arguments
    if(argc != 6) {
        std::cerr << "This program expects exactly 5 arguments. Please supply "
                << "'input file name', 'k-mer size', 'top count', 'fairness constant'"
                << " and 'output file name' in the given order." << endl;

        std::cerr << "The default values are: " << endl;
        std::cerr << "-input file name: '' " << endl;
        std::cerr << "-k-mer size: 30 " << endl;
        std::cerr << "-top count: 25 " << endl;
        std::cerr << "-fairness constant: 1.25 (should be between 1 - 2)" << endl;
        std::cerr << "-output file name: 'output.txt' " << endl;

        return 1;
    }

    else {
        std::cout << "Given arguments are: " << endl;
        for (int i = 1; i < argc; ++i) {
            std::cout << argv[i] << std::endl;
        }
        std::cout << "--------" << endl;
        std::cout << endl;

        input_filename = argv[1];
        KMER_SIZE = stoi(argv[2]);
        N = stoi(argv[3]);
        FAIRNESS_CONST = stof(argv[4]);
        output_filename = argv[5];
    }


    // Proceed with processing the given file
    clock_t begin = clock();
    unordered_map<kmer, int> kmer_map;
    int most_common_kmer_count = 0;

    // Reading .fastq file
    ifstream file(input_filename, ios::in);
    if(!file.good()) {
        std::cout << "There is a problem with the file, please check it." << endl;
    }

    else {
        string line;

        while(!file.eof()) {
            getline(file, line);
            LINE_COUNT++;

            if(LINE_COUNT % 100000 == 0)	//(DEBUG)
                cout << "Processing.. at line " << LINE_COUNT << endl; //(DEBUG)

            // If it is not a line containing DNA data, skip it.
            if(LINE_COUNT % 4 != 2) {
                continue;
            }

            // Else, proceed with processing the line
            else {
                if(!line.empty()) {
                    istringstream ss(line);
                    string dna_seq;
                    ss >> dna_seq;
                    int temp = process_seq(dna_seq, kmer_map, most_common_kmer_count);

                    if(temp == -1) {
                        std::cout << "Given DNA sequence is shorter than " << to_string(KMER_SIZE) << " chars. Aborting.." << endl;
                        break;
                    }

                    // Update most common kmer count
                    if (temp > most_common_kmer_count) {
                        most_common_kmer_count = temp;
                    }

                    //std::cout << dna_seq << endl; //(DEBUG)
                }

                // Clean map by erasing least common substrings
                if(kmer_map.size() >= MAP_SIZE_LIMIT) {
                    clean_map(kmer_map, most_common_kmer_count);
                }

                line.clear();
            }
        }

        file.close();
    }

    clean_map(kmer_map, most_common_kmer_count);
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << endl;
    std::cout << "Substring search took: " << elapsed_secs << " sec" << endl;

    write_map_to_file(kmer_map, output_filename);
    extract_top_N_kmers(kmer_map, N);

    return 0;
}
