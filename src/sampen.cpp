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

#include "utils.h"
#include "RangeTree2.h"
#include "sampen_calculator.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;


namespace RT = RangeTree;

struct stat {
    char *filename;
    unsigned m;
    int r;
    bool t;
} _stat;

double calculate_sampen_rangetree(
    const vector<int> &data, unsigned m, int r);
double calculate_sampen_rangetree_random(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_size, unsigned sample_num);
double calculate_sampen_rangetree_hist(
    const vector<int> &data, unsigned m, int r, double sample_rate);

void phelp(char *arg0);
void pargs(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    _stat.filename = NULL;
    _stat.m = 0;
    _stat.r = -1;
    _stat.t = false;
    pargs(argc, argv);

    unsigned long N;
    int *_data = readdata(_stat.filename, &N);
    vector<int> data(_data, _data + N);
    free(_data);

    cout << argv[0];
    cout << "\t_stat.filename: " << _stat.filename << endl;
    cout << "\t_stat.m: " << _stat.m << endl;
    cout << "\t_stat.r: " << _stat.r << endl;
    cout << "\t_stat.t: " << _stat.t << endl;
    cout << "\tdata length: " << N << endl;

    // Compute sample entropy by direct method
    double result = calculate_sampen_rangetree(data, _stat.m, _stat.r);
    cout << "Range Tree: SampEn(" << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result << endl;

    // Compute sample entropy by random sampling
    unsigned sample_size = 500;
    unsigned sample_num = 50;
    double result_random = calculate_sampen_rangetree_random(
        data, _stat.m, _stat.r, sample_size, sample_num);
    cout << "Range Tree (random): SampEn(" ;
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_random << endl;

    double error = (result_random - result) / result;
    cout << "Error (random) = " << error << endl;

    // Compute sample entropy by sampling according to histogram
    double result_histogram = calculate_sampen_rangetree_hist(
        data, _stat.m, _stat.r, 0.1);
    cout << "Range Tree (histogram): SampEn("; 
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_histogram << endl;

    error = (result_histogram - result) / result;
    cout << "Error (histogram) = " << error << endl;

    return 0;
}

double calculate_sampen_rangetree(
    const vector<int> &data, const unsigned m, const int r)
{
    sampen_calculator_rt scrt;
    return scrt.compute_entropy(data, m, r);
}

double calculate_sampen_rangetree_random(
    const vector<int> &data, const unsigned m, const int r, 
    const unsigned sample_size, const unsigned sample_num) 
{
    sampen_calculator_qr scqr(sample_num, sample_size);
    return scqr.compute_entropy(data, m, r);
}

double calculate_sampen_rangetree_hist(
    const vector<int> &data, const unsigned m, 
    const int r, double sample_rate)
{
    sampen_calculator_hg schg(sample_rate);
    return schg.compute_entropy(data, m, r);
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

void pargs(int argc, char *argv[]) 
{
    if (argc < 1) phelp(argv[0]);
    int i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            if (_stat.filename) phelp(argv[0]);
            _stat.filename = argv[i];
            i++;
            break;
        }

        if (strlen(argv[i]) != 2) {
            cerr << "unrecognized option " << argv[i] << endl;
            phelp(argv[0]);
        }
        if (i + 1 >= argc) {
            cerr << "argument of option " << argv[i] << " is missing" << endl;
            phelp(argv[0]);
        }
        switch (argv[i][1]) {
        case 'm':
            if (_stat.m) {
                cerr << "template length m has been set yet\n";
                phelp(argv[0]);
            }
            _stat.m = (unsigned) atoi(argv[i+1]);
            break;
        case 'r':
            if (_stat.r != -1) {
                cerr << "tolerance r has been set yet\n";
                phelp(argv[0]);
            }
            _stat.r = (unsigned) atoi(argv[i+1]);
            break;
        default:
            cerr << "unrecognized option " << argv[i] << endl;
            phelp(argv[0]);
            break;
        }
        i = i + 2;
    }
    if (_stat.filename == NULL) {
        cerr << "INPUT_FILENAME is required\n";
        phelp(argv[0]);
    }
    if (_stat.m == 0) _stat.m = 3;
    if (_stat.r == -1) _stat.r = 100;
}
