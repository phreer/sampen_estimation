/* file: sampen_var.cpp
 * date: 2020-02-10
 * author: phree
 *
 * description: This program is used to measure the variance of the 
 *   errors of the entropies computed by sampling method.
 */

#include <iostream>
#include <vector>
#include <algorithm>

#include "sampen_calculator.h"

using std::vector;
using std::cout;
using std::endl;

struct  
{
    unsigned m = 2;
    int r = 30;
    char *filename = NULL;
    double sample_rate = -1.;
    unsigned sample_size = 0;
    unsigned sample_num = 0;
    unsigned rounds = 0;
} _status;

void parse_args(int argc, char *argv[])
{
    if (argc < 2) 
    throw std::invalid_argument("Please specify an input file");
    _status.filename = argv[argc - 1];

    ArgumentParser ap(argc, argv);
    string arg;
    arg = ap.getArg("-m");
    if (arg.size()) 
        _status.m = std::stoi(arg);
    else _status.m = 2;
    arg = ap.getArg("-r");
    if (arg.size())
        _status.r = std::stoi(arg);
    else _status.r = 30;
    arg = ap.getArg("-sr");
    if (arg.size()) 
    {
        _status.sample_rate = std::stod(arg);
    } 
    else {
        arg = ap.getArg("-sn");
        if (arg.size()) 
            _status.sample_num = std::stoi(arg);
        else throw std::invalid_argument(
            "Please specify a sample rate or a sample num");
        arg = ap.getArg("-ss");
        if (arg.size()) 
            _status.sample_size = std::stoi(arg);
        else throw std::invalid_argument(
            "Please specify a sample rate or a sample size");
    }
    arg = ap.getArg("-rounds");
    if (arg.size()) 
        _status.rounds = std::stoi(arg);
    else 
        _status.rounds = 20;
    if (_status.rounds == 1) 
        throw std::invalid_argument("rounds should be greater than 1");
}


int main(int argc, char *argv[])
{
    parse_args(argc, argv);

    unsigned long N;
    int *_data = readdata(_status.filename, &N);
    vector<int> data(_data, _data + N);
    free(_data);

    unsigned sample_size, sample_num;
    if (_status.sample_rate < 0) 
    {
        sample_size = _status.sample_size;
        sample_num = _status.sample_num;
    }
    else {
        sample_size = static_cast<unsigned>(_status.sample_rate * N);
        sample_num = static_cast<unsigned>(1. / _status.sample_rate);
    }

    cout << "+---------------------------------------------------+" << endl;
    cout << "|                 sampen_var starts                 |" << endl;
    cout << "+---------------------------------------------------+" << endl;
    cout << "m: " << _status.m << endl;
    cout << "r: " << _status.r << endl;
    cout << "rounds: " << _status.rounds << endl;
    cout << "sample size: " << sample_size << endl;
    cout << "sample num: " << sample_num << endl;
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout.precision(6);

    SampenCalculatorD sc;
    double ground_truth = sc.ComputeEntropy(
        data, _status.m, _status.r, nullptr, nullptr);
    std::cout << "SampleEntropy(" << N << ", " << _status.m << ", " << _status.r;
    std::cout << ") = " << ground_truth << std::endl;
    vector<double> results(_status.rounds);
    vector<double> errors(_status.rounds);
    for (unsigned i = 0; i < _status.rounds; i++)
    {
        SampenCalculatorQR sc(sample_num, sample_size, true);
        results[i] = sc.ComputeEntropy(
            data, _status.m, _status.r, nullptr, nullptr);
        
        errors[i] = results[i] - ground_truth;
        std::cout << "Sampen: " << results[i] << ", ";
        std::cout << "Error: " << (errors[i]) / ground_truth << std::endl;
    }
    double mean = 0.;
    std::for_each(errors.cbegin(), 
                  errors.cend(), 
                  [&](const double r) { mean += r; });
    mean /= _status.rounds * ground_truth;
    std::cout << "mean: " << mean << std::endl;
    double var = 0.;
    std::for_each(results.cbegin(), 
                  results.cend(), 
                  [&](const double r) { var += r * r; });
    var /= (_status.rounds - 1);
    std::cout << "The variance is " << var << std::endl;
    cout << "+---------------------------------------------------+" << endl;
    cout << "|                  sampen_var ends                  |" << endl;
    cout << "+---------------------------------------------------+" << endl;
}