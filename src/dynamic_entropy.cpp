
/* file: sampen.cpp
 * date: 2020-3-3
 * author: phree
 *
 * description: program used compute dynamic entropy
 */

#include <iostream>
#include <time.h>
#include <algorithm>
#include <string.h>

#include "math.h"
#include "utils.h"
#include "RangeTree2.h"
#include "sampen_calculator.h"
#include "entropy_calculator.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

namespace RT = RangeTree;

struct stat {
    char filename[256];
    unsigned m;
    int r;
    unsigned ml;
} _stat;

double calculate_sampen_direct(const vector<int> &data, unsigned m, int r);
double ComputeDynamicEntropy(const vector<int> &data, unsigned m, unsigned max_level);

void phelp(char *arg0);
void ParseArgs(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    _stat.m = 0;
    _stat.r = -1;
    ParseArgs(argc, argv);

    unsigned long N;
    int *_data = readdata(_stat.filename, &N);
    vector<int> data(_data, _data + N);
    free(_data);
    for (int i = 0; i < _stat.m; i++) data.push_back(data.back());

    cout << argv[0];
    cout << "\t_stat.filename: " << _stat.filename << endl;
    cout << "\t_stat.m: " << _stat.m << endl;
    cout << "\t_stat.r: " << _stat.r << endl;
    cout << "\t_stat.ml: " << _stat.ml << endl;
    cout << "\tdata length: " << N << endl;

    // Compute sample entropy by direct method
    double result = calculate_sampen_direct(data, _stat.m, _stat.r);
    cout << "Direct: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result << endl;

    // Compute DynamicEntropy
    double dent = ComputeDynamicEntropy(data, _stat.m, _stat.ml);
    cout << "DynamicEntropy(" << _stat.m << ", " << _stat.ml << ", ";
    cout << N << ") = " << dent << endl;
    return 0;
}

double calculate_sampen_direct(
    const vector<int> &data, unsigned m, int r)
{
    sampen_calculator_d sc;
    return sc.compute_entropy(data, m, r);
}

double ComputeDynamicEntropy(const vector<int> &data, 
                             unsigned m, unsigned max_level) 
{
    DynamicEntropyCalculator dec;
    return dec.ComputeEntropy(data, m + 1, max_level);
}

void phelp(char *arg0) 
{
    char help[] = "options: \n"
                "\t-m M (default: 3) template length\n"
                "\t-r R (default: 100) tolerance\n";
    cerr << "usage: " << arg0 << "[options] INPUT_FILENAME\n";
    cerr << help;
    exit(-1);
}

void ParseArgs(int argc, char *argv[])
{
    ArgumentParser ap(argc, argv);
    string arg;
    arg = ap.getArg("-filename");
    if (arg.size()) strcpy(_stat.filename, arg.c_str());
    else throw std::invalid_argument(
        "Please specify an input file with -filename FILENAME");

    arg = ap.getArg("-m");
    if (arg.size()) _stat.m = std::stoi(arg);
    else throw std::invalid_argument("Please specify m with -m M");

    arg = ap.getArg("-ml");
    if (arg.size()) _stat.ml = std::stoi(arg);
    else throw std::invalid_argument("Please specify max level with -ml ML");

    arg = ap.getArg("-r");
    if (arg.size()) _stat.r = std::stoi(arg);
    else _stat.r = 30;
}

