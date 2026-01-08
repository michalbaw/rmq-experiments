#include <algorithm>
#include <getopt.h> 
#include <random>
#include <fstream>
#include <iostream>

using ll = long long;

namespace query {

struct config {
    ll N;
    ll q; 
    ll r;
    std::string ofile;
    
    void printConfig() {
        printf("=============\nConfiguration \n=============\n");
        printf("Sequence length: %lld\n", N);
        printf("Num Queries: %lld\n", q);
        printf("Query Range Size: %lld\n", r);
        printf("Output file: %s\n", ofile.c_str());
        printf("---------------------------\n");
    }
};

void
usage(char* program)
{
    printf("%s \n", program);
    printf("where\n");
    printf("  -N    : size of the sequence\n");
    printf("  -q    : generating q queries for sequence.\n");
    printf("  -r    : range of the query interval \n");
    printf("  -f    : output file for query\n");
    exit(EXIT_FAILURE);
};

config parse_args(int argc, char* const argv[]) {
    config con;
    int op;
    while ((op = getopt(argc, argv, "n:q:r:f:")) != -1) {
        switch (op) {
            case 'n': con.N = std::stoull(std::string(optarg));
            break;
            case 'q': con.q = std::stoull(std::string(optarg));
            break;
            case 'r': con.r = std::stoull(std::string(optarg));
            break;
            case 'f': con.ofile = optarg;
            break;
            case '?':
            default:
                usage(argv[0]);
        }
    }
    return con;
}

}

int log_10(int x) {
    int a = 0;
    while (x > 1) {
        a += 1;
        x /= 10;
    }
    return a; 
}

int randomRange(int l, int r, auto rng) {
    int x = rng();
    int offset = x % (r - l);
    if (offset < 0) {
        offset += (r - l);
    }
    return l + offset;
}

const int tens[10] = {1, (int)1e1, (int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6, (int)1e7, (int)1e8, (int)1e9};

int main(int argc, char* const argv[]) {
    
    
    std::ios::sync_with_stdio(false);
    
    query::config con = query::parse_args(argc,argv);
    con.printConfig();
    
    std::random_device rd;
    std::mt19937 gen(rd());
    printf("Generating random queries...\n");
    std::ofstream os;
    os.open(con.ofile);
    os << con.q << "\n";
    
    // Use 0 for generating uniformly on the logarithm of the range
    if (con.r == 0) {
        std::uniform_int_distribution<size_t> rangeDistrib(1, log_10(con.N));
        std::uniform_int_distribution<size_t> distrib(0, con.N - 1);

        for (int i = 0; i < con.q; i += 1) {
            auto range = rangeDistrib(gen);
            auto length = randomRange(tens[range - 1], tens[range], gen);
            auto l = randomRange(0, con.N - length, gen);
            os << l << " " << l + length << '\n';
        }
    } else {
        std::uniform_int_distribution<ll> dis(0, con.N-con.r);
        for(size_t i = 0; i < con.q; ++i) {
            ll i1 = dis(gen), i2 = i1 + con.r - 1;
            os << i1 << " " << i2 << "\n";
        }
    }
    os.close();
    
    printf("Queries are written to %s\n", con.ofile.c_str());
    
    printf("Finish!\n");
    
    
    return 0;
}
