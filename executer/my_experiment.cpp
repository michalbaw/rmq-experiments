#include <sdsl/rmq_support.hpp> // include header for range minimum queries
#include <sdsl/rmq_bitstack_fast.hpp> // include header for range minimum queries
#include <sdsl/int_vector.hpp>
#include "sdsl/memory_management.hpp"
#include <algorithm>
#include <cmath>
#include <getopt.h> 
#include <fstream>
#include <iostream>
#include <chrono>
#include <climits>
#include <random>
#include <vector>
 
#include "../rmq/includes/RMQRMM64.h"
#include "../succinct/cartesian_tree.hpp"
#include "../succinct/mapper.hpp"
#include "hardware_event.h"
#include <limits.h>

using namespace std;
using namespace std::chrono;

#define TIME_EACH_QUERY true
constexpr int INPUT_AND_QUERIES_GENERATIONS = 5;
constexpr int INPUT_ARRAY_LENGTH = 1e7;
constexpr int LOG_INPUT_ARRAY_LENGTH = 7;
constexpr int NUMBER_OF_QUERIES = 1e5;
using input_array_type = long long;
constexpr input_array_type MIN_INPUT_ARRAY_VALUE = 0;
constexpr input_array_type MAX_INPUT_ARRAY_VALUE = LLONG_MAX;

random_device rd;
mt19937 rng(rd());

vector<pair<size_t, size_t>> getRandomQueriesUniform(int number_of_queries) {
    uniform_int_distribution<size_t> distrib(0, INPUT_ARRAY_LENGTH - 1);
    vector<pair<size_t, size_t>> queries(number_of_queries);
    for (int i = 0; i < number_of_queries; i += 1) {
        auto r0 = distrib(rng);
        auto r1 = distrib(rng);
        auto [l, r] = minmax(r0, r1);
        queries[i] = {l, r + 1};
    }
    return queries;
}

int randomRange(int l, int r) {
    int x = rng();
    int offset = x % (r - l);
    if (offset < 0) {
        offset += (r - l);
    }
    return l + offset;
}

const int tens[10] = {1, (int)1e1, (int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6, (int)1e7, (int)1e8, (int)1e9};
pair<vector<pair<size_t, size_t>>, vector<size_t>> getRandomQueriesPerRange(int number_of_queries) {
    uniform_int_distribution<size_t> rangeDistrib(1, LOG_INPUT_ARRAY_LENGTH);
    uniform_int_distribution<size_t> distrib(0, INPUT_ARRAY_LENGTH - 1);

    vector<size_t> ranges(number_of_queries);
    vector<pair<size_t, size_t>> queries(number_of_queries);
    for (int i = 0; i < number_of_queries; i += 1) {
        auto range = rangeDistrib(rng);
        auto length = randomRange(tens[range - 1], tens[range]);
        auto l = randomRange(0, INPUT_ARRAY_LENGTH - length);
        ranges[i] = range;
        queries[i] = {l, l + range};
    }
    return {queries, ranges};
}

volatile long long output_variable = 0; 

vector<input_array_type> getRandomInputArray(int input_array_size) {
    uniform_int_distribution<input_array_type> distrib(MIN_INPUT_ARRAY_VALUE, MAX_INPUT_ARRAY_VALUE);
    vector<input_array_type> input_array(input_array_size);
    for (int i = 0; i < input_array_size; i += 1) {
        input_array[i] = distrib(rng);
    }
    return input_array;
}

void experimentRandomAccess(const vector<pair<size_t, size_t>>& queries, const vector<input_array_type>& input_array, bool print_benchmark_results = true, vector<size_t> ranges = {}) {
    #if TIME_EACH_QUERY
    if (ranges.size() < queries.size()) {
        std::cout << "Provide ranges";
        return;
    }
    #endif
    auto start = high_resolution_clock::now();
    for (int i = 0; i < queries.size(); ++i) {
        #if TIME_EACH_QUERY
        auto start = high_resolution_clock::now();
        #endif
        output_variable = input_array[queries[i].first];
        output_variable = input_array[queries[i].second];
        #if TIME_EACH_QUERY
        auto end = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cout << INPUT_ARRAY_LENGTH << "," << "RANDOM_ACCESS," << duration.count() << "," << ranges[i] << std::endl;
        #endif
    }
    #if !TIME_EACH_QUERY
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    if (print_benchmark_results) {
        std::cout << INPUT_ARRAY_LENGTH << "," << "RANDOM_ACCESS," << duration.count() << std::endl;
    }
    #endif
}

template<class rmq_sdsl>
void experimentSDSL(const std::string& algo_name, const vector<pair<size_t, size_t>>& queries, const vector<input_array_type>& input_array, bool print_benchmark_results = true, vector<size_t> ranges = {}) {
    #if TIME_EACH_QUERY
    if (ranges.size() < queries.size()) {
        std::cout << "Provide ranges";
        return;
    }
    #endif
    sdsl::int_vector<> int_vec(INPUT_ARRAY_LENGTH, 0, 8 * sizeof(long long)); 
    for (int i = 0; i < INPUT_ARRAY_LENGTH; ++i) { int_vec[i] = input_array[i];}
    rmq_sdsl rmq(&int_vec);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < queries.size(); ++i) {
        #if TIME_EACH_QUERY
        auto start = high_resolution_clock::now();
        #endif
        output_variable = rmq(queries[i].first, queries[i].second);
        #if TIME_EACH_QUERY
        auto end = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cout << INPUT_ARRAY_LENGTH << "," << algo_name << "," << duration.count() << "," << ranges[i] << std::endl;
        #endif
    }
    #if !TIME_EACH_QUERY
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    if (print_benchmark_results) {
        std::cout << INPUT_ARRAY_LENGTH << "," << algo_name << "," << duration.count() << std::endl;
    }
    #endif
}

void generateInputDataAndRunExperiments(bool print_benchmark_results) {
    auto [queries, ranges] = getRandomQueriesPerRange(NUMBER_OF_QUERIES);
    auto input_array = getRandomInputArray(INPUT_ARRAY_LENGTH);
    experimentRandomAccess(queries, input_array, print_benchmark_results, ranges);
    experimentSDSL<sdsl::RMQ_Fast>("UJ_LIB_SDSL", queries, input_array, print_benchmark_results, ranges);
    // experimentSDSL<sdsl::rmq_succinct_sct<>>("SDSL_SCT", queries, input_array, print_benchmark_results, ranges);
    experimentSDSL<sdsl::rmq_succinct_rec_new<true, 0, 1024,128,0>>("RMQ_SDSL_REC", queries, input_array, print_benchmark_results, ranges);
    experimentSDSL<sdsl::rmq_succinct_rec_new<true, 2048, 1024,128,0>>("RMQ_SDSL_REC_ST", queries, input_array, print_benchmark_results, ranges);
}

int main() {
    #if TIME_EACH_QUERY
    std::cout << "N,Algo,Time,Range" << std::endl; 
    #else
    std::cout << "N,Algo,Time" << std::endl; 
    #endif
    // Ignore first generations they seem to work slower.
    // Don't know why because each iteration new data is generated
    generateInputDataAndRunExperiments(false);
    for (int i = 0; i < INPUT_AND_QUERIES_GENERATIONS; ++i) {
        generateInputDataAndRunExperiments(true);
    }
}