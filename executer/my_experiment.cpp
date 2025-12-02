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

constexpr int INPUT_AND_QUERIES_GENERATIONS = 1e2;
constexpr int INPUT_ARRAY_LENGTH = 1e8;
constexpr int NUMBER_OF_QUERIES = 1e5;
using input_array_type = long long;
constexpr input_array_type MIN_INPUT_ARRAY_VALUE = 0;
constexpr input_array_type MAX_INPUT_ARRAY_VALUE = LLONG_MAX;

random_device rd;
mt19937 rng(rd());

vector<pair<size_t, size_t>> getRandomQueries(int number_of_queries) {
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

volatile long long output_variable = 0; 

vector<input_array_type> getRandomInputArray(int input_array_size) {
    uniform_int_distribution<input_array_type> distrib(MIN_INPUT_ARRAY_VALUE, MAX_INPUT_ARRAY_VALUE);
    vector<input_array_type> input_array(input_array_size);
    for (int i = 0; i < input_array_size; i += 1) {
        input_array[i] = distrib(rng);
    }
    return input_array;
}

void experimentRandomAccess(const vector<pair<size_t, size_t>>& queries, const vector<input_array_type>& input_array, bool print_benchmark_results = true) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < queries.size(); ++i) {
        output_variable = input_array[queries[i].first];
        output_variable = input_array[queries[i].second];
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (print_benchmark_results) {
        std::cout << "RANDOM_ACCESS," << duration.count() << '\n';
    }
}

void experimentUJLibSDSL(const vector<pair<size_t, size_t>>& queries, const vector<input_array_type>& input_array, bool print_benchmark_results = true) {
    sdsl::int_vector<> int_vec(INPUT_ARRAY_LENGTH, 0, 8 * sizeof(long long)); 
    for (int i = 0; i < INPUT_ARRAY_LENGTH; ++i) { int_vec[i] = input_array[i];}
    sdsl::RMQ_Fast rmq(&int_vec);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < queries.size(); ++i) {
        output_variable = rmq(queries[i].first, queries[i].second);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (print_benchmark_results) {
        std::cout << "UJ_LIB_SDSL_VECTOR," << duration.count() << '\n';
    }
}

void experimentSDSLSCT(const vector<pair<size_t, size_t>>& queries, const vector<input_array_type>& input_array, bool print_benchmark_results = true) {
    sdsl::int_vector<> int_vec(INPUT_ARRAY_LENGTH, 0, 8 * sizeof(long long)); 
    for (int i = 0; i < INPUT_ARRAY_LENGTH; ++i) { int_vec[i] = input_array[i];}
    sdsl::rmq_succinct_sct<> rmq(&int_vec);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < queries.size(); ++i) {
        output_variable = rmq(queries[i].first, queries[i].second);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (print_benchmark_results) {
        std::cout << "SDSL_SCT," << duration.count() << '\n';
    }
}

void generateInputDataAndRunExperiments(bool print_benchmark_results) {
    auto queries = getRandomQueries(NUMBER_OF_QUERIES);
    auto input_array = getRandomInputArray(INPUT_ARRAY_LENGTH);
    experimentRandomAccess(queries, input_array, print_benchmark_results);
    experimentUJLibSDSL(queries, input_array, print_benchmark_results);
    experimentSDSLSCT(queries, input_array, print_benchmark_results);
}

int main() {
    // Ignore first generations seems to work slower. 
    // Don't know why because each iteration new data is generated 
    generateInputDataAndRunExperiments(false);
    for (int i = 0; i < INPUT_AND_QUERIES_GENERATIONS; ++i) {
        generateInputDataAndRunExperiments(true);
    }
}