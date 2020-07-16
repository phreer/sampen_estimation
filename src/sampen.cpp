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

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

namespace RT = RangeTree;

struct stat {
    char filename[256];
    unsigned m;
    double r;
    unsigned sample_num;
    unsigned sample_size;
} _stat;

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
    for (unsigned i = 0; i < _stat.m; i++) data.push_back(data.back());
    double var = ComputeVarience(data);
    int r = static_cast<int>(_stat.r * sqrt(var));

    double sample_rate = static_cast<double>(_stat.sample_size) / N;
    
    cout << argv[0];
    cout << "+------------------------ Settings ------------------------+";
    cout << std::endl;
    cout << "\t_stat.filename: " << _stat.filename << endl;
    cout << "\t_stat.m: " << _stat.m << endl;
    cout << "\t_stat.r: " << _stat.r << endl;
    cout << "\tr_scaled: " << r << endl;
    cout << "\t_stat.sample_num: " << _stat.sample_num << endl;
    cout << "\t_stat.sample_size: " << _stat.sample_size << endl;
    cout << "\tsample rate: " << sample_rate << endl;
    cout << "\tdata length: " << N << endl;
    cout << "\tvariance: " << var << endl;
    cout << "+----------------------------------------------------------+";
    cout << std::endl;
    // Compute sample entropy by direct method
    double result = ComputeSampenDirect(data, _stat.m, r, nullptr, nullptr);
    cout << "Direct: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result << endl;

    // Compute sample entropy by quasi-random sampling
    double result_random = ComputeSampenQR2(
        data, _stat.m, r, _stat.sample_size, _stat.sample_num, 
        nullptr, nullptr);
    cout << "Quasi-random: SampEn(" ;
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_random << endl;

    double error = result_random - result;
    cout << "Error = " << error;
    cout << ", Relative Error (Quasi-random) = " << error / result << endl;

    //Compute sample entropy by random sampling
    result_random = ComputeSampenQR(
        data, _stat.m, r, _stat.sample_size, _stat.sample_num, 
        nullptr, nullptr);
    cout << "Quasi-random (Sorting): SampEn(" ;
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_random << endl;

    error = result_random - result;
    cout << "Error = " << error;
    cout << ", Relative Error (Quasi-random (Sorting)) = ";
    cout << error / result << endl;

    // Compute sample entropy using corset
    result_random = ComputeSampenCoreset(data, _stat.m, r, 
        _stat.sample_size, _stat.sample_num, nullptr, nullptr);
    cout << "Coreset: SampEn(" ;
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_random << endl;

    error = result_random - result;
    cout << "Error = " << error;
    cout << ", Relative Error (Coreset) = " << error / result << endl;

    // Compute sample entropy using unifrom distribution
    result_random = ComputeSampenUniform(data, _stat.m, r, 
        _stat.sample_size, _stat.sample_num, nullptr, nullptr);
    cout << "Uniform: SampEn(" ;
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_random << endl;

    error = result_random - result;
    cout << "Error = " << error;
    cout << ", Relative Error (Uniform) = " << error / result << endl;

    // /* Compute sample entropy by sampling according to kd tree */
    // double result_hist = calculate_sampen_nkdtree_hist(data, _stat.m, r, 
    //     _stat.sample_num, _stat.sample_size);
    // cout << "KD Tree Sampling: SampEn(";
    // cout << _stat.m << ", " << _stat.r << ", " << _stat.sample_num << ", ";
    // cout << _stat.sample_size << ") = " << result_hist << endl;
    // error = result_hist - result;
    // cout << "Error = " << error;
    // cout << ", Relative Error (kd-tree histogram) = " << error / result << endl;

    // // Compute sample entropy by kd tree
    // result = calculate_sampen_kdtree(data, _stat.m, r);
    // cout << "kd tree: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result << endl;

    // // Compute sample entropy by kd tree (grid)
    // result = calculate_sampen_kdtree_grid(data, _stat.m, r);
    // cout << "kd tree (grid): SampEn(" << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result << endl;

    // // Compute sample entropy by range tree
    // result = calculate_sampen_rangetree(data, _stat.m, r);
    // cout << "Range Tree: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result << endl;

    // Compute sample entropy by sampling according to histogram
    // double result_hist = calculate_sampen_rangetree_hist(
    //     data, _stat.m, r, 0.1);
    // cout << "Histogram: SampEn("; 
    // cout << _stat.m << ", " << _stat.r << ", ";
    // cout << N << ") = " << result_hist << endl;

    // error = (result_hist - result) / result;
    // cout << "Error (histogram) = " << error << endl;
    return 0;
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
    if (arg.size()) _stat.r = std::stod(arg);
    else _stat.r = 0.1;
    
    arg = ap.getArg("-sample_num");
    if (arg.size()) _stat.sample_num = std::stoi(arg);
    else throw std::invalid_argument(
        "Please specify a sample num with -sample_num SAMPLE_NUM");
    
    arg = ap.getArg("-sample_size");
    if (arg.size()) _stat.sample_size = std::stoi(arg);
    else throw std::invalid_argument(
        "Please specify a sample num with -sample_size SAMPLE_SIZE");
}