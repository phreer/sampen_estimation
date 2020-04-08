/* file: sampen.cpp
 * date: 2019-11-05
 * author: phree
 *
 * description: compute sample entropy using range tree, implemented by
 * Luca Weihs
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
    unsigned sample_num;
    unsigned sample_size;
} _stat;

double calculate_sampen_direct(const vector<int> &data, unsigned m, int r);
double ComputeDynamicEntropy(const vector<int> &data, unsigned m, unsigned max_level);
double calculate_sampen_rangetree(
    const vector<int> &data, unsigned m, int r);
double calculate_sampen_kdtree(
    const vector<int> &data, unsigned m, int r);
double calculate_sampen_kdtree_grid(
    const vector<int> &data, unsigned m, int r);
double calculate_sampen_rangetree_random(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_size, unsigned sample_num);
double calculate_sampen_rangetree_hist(
    const vector<int> &data, unsigned m, int r, double sample_rate);
double calculate_sampen_nkdtree_hist(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_num, unsigned sample_size);
void phelp(char *arg0);
void ParseArgs(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    _stat.m = 0;
    _stat.r = -1;
    _stat.sample_num = 0;
    _stat.sample_size = 0;
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
    cout << "\t_stat.sample_num: " << _stat.sample_num << endl;
    cout << "\t_stat.sample_size: " << _stat.sample_size << endl;
    cout << "\tdata length: " << N << endl;

    // Compute sample entropy by direct method
    double result = calculate_sampen_direct(data, _stat.m, _stat.r);
    cout << "Direct: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result << endl;

    // Compute sample entropy by random sampling
    double sample_rate = static_cast<double>(_stat.sample_size) / N;
    double result_random = calculate_sampen_rangetree_random(
        data, _stat.m, _stat.r, _stat.sample_size, _stat.sample_num);
    cout << "Quasi-random: SampEn(" ;
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_random << endl;

    double error = result_random - result;
    cout << "Error = " << error;
    cout << ", Relative Error (Quasi-random) = " << error / result << endl;

    /* Compute sample entropy by sampling according to kd tree */
    double result_hist = calculate_sampen_nkdtree_hist(data, _stat.m, _stat.r, 
        _stat.sample_num, _stat.sample_size);
    cout << "KD Tree Sampling: SampEn(";
    cout << _stat.m << ", " << _stat.r << ", " << _stat.sample_num << ", ";
    cout << _stat.sample_size << ") = " << result_hist << endl;
    error = result_hist - result;
    cout << "Error = " << error;
    cout << ", Relative Error (kd-tree histogram) = " << error / result << endl;

    // // Compute sample entropy by kd tree
    // result = calculate_sampen_kdtree(data, _stat.m, _stat.r);
    // cout << "kd tree: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result << endl;

    // // Compute sample entropy by kd tree (grid)
    // result = calculate_sampen_kdtree_grid(data, _stat.m, _stat.r);
    // cout << "kd tree (grid): SampEn(" << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result << endl;

    // // Compute sample entropy by range tree
    // result = calculate_sampen_rangetree(data, _stat.m, _stat.r);
    // cout << "Range Tree: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result << endl;

    // Compute sample entropy by sampling according to histogram
    // double result_hist = calculate_sampen_rangetree_hist(
    //     data, _stat.m, _stat.r, 0.1);
    // cout << "Histogram: SampEn("; 
    // cout << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result_hist << endl;

    // error = (result_hist - result) / result;
    // cout << "Error (histogram) = " << error << endl;
    return 0;
}

double calculate_sampen_direct(
    const vector<int> &data, unsigned m, int r)
{
    sampen_calculator_d sc;
    return sc.compute_entropy(data, m, r);
}

double calculate_sampen_rangetree(
    const vector<int> &data, unsigned m, int r)
{
    sampen_calculator_rt sc;
    return sc.compute_entropy(data, m, r);
}

double calculate_sampen_kdtree(
    const vector<int> &data, unsigned m, int r) 
{
    sampen_calculator_kd sc;
    return sc.compute_entropy(data, m, r);
}

double calculate_sampen_nkdtree_hist(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_size, unsigned sample_num) 
{
    sampen_calculator_nkd sc(sample_size, sample_num);
    return sc.compute_entropy(data, m, r);
}

double calculate_sampen_kdtree_grid(
    const vector<int> &data, unsigned m, int r) 
{
    sampen_calculator_kdg sc;
    return sc.compute_entropy(data, m, r);
}

double calculate_sampen_rangetree_random(
    const vector<int> &data, const unsigned m, const int r, 
    const unsigned sample_size, const unsigned sample_num) 
{
    sampen_calculator_qr sc(sample_num, sample_size);
    return sc.compute_entropy(data, m, r);
}

double calculate_sampen_rangetree_hist(
    const vector<int> &data, const unsigned m, 
    const int r, double sample_rate)
{
    sampen_calculator_hg sc(sample_rate);
    return sc.compute_entropy(data, m, r);
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

    arg = ap.getArg("-r");
    if (arg.size()) _stat.r = std::stoi(arg);
    else _stat.r = 30;

    arg = ap.getArg("-sample_num");
    if (arg.size()) _stat.sample_num = std::stoi(arg);
    else throw std::invalid_argument(
        "Please specify a sample num with -sample_num SAMPLE_NUM");
    
    arg = ap.getArg("-sample_size");
    if (arg.size()) _stat.sample_size = std::stoi(arg);
    else throw std::invalid_argument(
        "Please specify a sample num with -sample_size SAMPLE_SIZE");
}