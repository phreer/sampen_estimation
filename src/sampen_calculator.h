/* file: sampen_calculator.h
 * date: 2019-12-07
 * author: phree
 *
 * description: classes used to compute sample entropy
 */

#ifndef __SAMPEN_CALCULATOR_H__
#define __SAMPEN_CALCULATOR_H__

#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#include "random_sampler.h"

using std::vector;

double ComputeSampenAB(double A, double B, unsigned N, unsigned m)
{
    // std::cout << "A: " << A << ", B: " << B << std::endl;
    if (A > 0 && B > 0)
    {
        // std::cout << "B / A: " << B / A << std::endl;
        return -log(B / A);
    }
    else
        return -log((N - m - 1) / (N - m));
}

// base class to calculate sample entropy
class SampenCalculator
{
private:
    void _CheckDim(const vector<int> &data, unsigned m)
    {
        if (data.size() <= m)
            throw std::invalid_argument("data.size() < m");
    }
    virtual vector<long long> _ComputeAB(
        const vector<int> &data, unsigned m, int r) = 0;

public:
    virtual double ComputeEntropy(const vector<int> &data,
                                  unsigned m, int r, double *a, double *b)
    {
        _CheckDim(data, m);
        // auto start = std::chrono::system_clock::now();
        vector<long long> AB = _ComputeAB(data, m, r);
        // auto end = std::chrono::system_clock::now();
        // std::chrono::duration<double> interval = end - start;
        // std::cout << "time: " << interval.count() << "s" << std::endl;

        unsigned N = data.size();
        unsigned sample_num = AB.size() / 2;
        double result = 0;

        long long A = 0;
        long long B = 0;
        for (unsigned i = 0; i < sample_num; i++)
        {
            A += AB[i * 2];
            B += AB[i * 2 + 1];
        }
        if (a) *a = A / sample_num;
        if (b) *b = B / sample_num;
        result = ComputeSampenAB(A, B, N, m);

        return result;
    }
};

// direct method
class SampenCalculatorD : public SampenCalculator
{
private:
    virtual vector<long long> _ComputeAB(
        const vector<int> &data, unsigned m, int r) override;
};

// range tree
class SampenCalculatorRT : public SampenCalculator
{
private:
    virtual vector<long long> _ComputeAB(
        const vector<int> &data, unsigned m, int r) override;
};

// old kd tree
class SampenCalculatorKD: public SampenCalculator
{
private:
    virtual vector<long long> _ComputeAB(
        const vector<int> &data, unsigned m, int r) override;
};

// Compute sample entropy by kd tree divided according to grid
class SampenCalculatorKDG: public SampenCalculator
{
private:
    virtual vector<long long> _ComputeAB(
        const vector<int> &data, unsigned m, int r) override;
};

// Uniform distribution sampling
class SampenCalculatorUniform: public SampenCalculator
{
public:
    explicit SampenCalculatorUniform(
        unsigned sample_num_, unsigned sample_size_, bool random = false)
        : sample_num(sample_num_), sample_size(sample_size_), random(random) 
    {}
    void set_sample_num(unsigned sample_num_) { sample_num = sample_num_; }
    void set_sample_size(unsigned sample_size_)
    {
        sample_size = sample_size_;
    }

private:
    virtual vector<long long> _ComputeAB(const vector<int> &data,
                                          unsigned m, int r) override;
    unsigned sample_num;
    unsigned sample_size;
    bool random;
};

// quasi-random sampling
class SampenCalculatorQR: public SampenCalculator
{
public:
    explicit SampenCalculatorQR(unsigned sample_num_, unsigned sample_size_, 
                                bool presort = true, bool random = false) 
        : sample_num(sample_num_), sample_size(sample_size_), 
        presort(presort), random(random) 
    {}
    void set_sample_num(unsigned sample_num_) { sample_num = sample_num_; }
    void set_sample_size(unsigned sample_size_)
    {
        sample_size = sample_size_;
    }

private:
    virtual vector<long long> _ComputeAB(const vector<int> &data,
                                          unsigned m, int r) override;
    unsigned sample_num;
    unsigned sample_size;
    bool random;
    bool presort;
};


class SampenCalculatorCoreset 
{
public:
    explicit SampenCalculatorCoreset(
        unsigned sample_num_, unsigned sample_size_, bool random = false)
        : sample_num(sample_num_), sample_size(sample_size_), random(random) 
    {}
    void set_sample_num(unsigned sample_num_) { sample_num = sample_num_; }
    void set_sample_size(unsigned sample_size_)
    {
        sample_size = sample_size_;
    }
    double ComputeSampen(const vector<int> &data, unsigned m, int r, 
                         double *a, double *b);
private:
    vector<double> _ComputeAB(const vector<int> &data,
                              unsigned m, int r);
    unsigned sample_num;
    unsigned sample_size;
    bool random;
};

// Compute sample entropy using new kd tree
// Here, we require sample_size = 2^n for some non-negative integer n
class SampenCalculatorNKD : public SampenCalculator
{
public:
    explicit SampenCalculatorNKD(
        unsigned sample_num, unsigned sample_size)
        : sample_num_(sample_num), sample_size_(sample_size) 
        {
            if (!IsPowerTwo(sample_size_)) 
                throw std::invalid_argument("sample_size should be power of 2");
        }
private:
    virtual vector<long long> _ComputeAB(const vector<int> &data,
                                          unsigned m, int r) override;
    unsigned sample_num_;
    unsigned sample_size_;
};

// Compute sampen by sample using histogram
class SampenCalculatorHG : public SampenCalculator
{
private:
    double _sample_rate;
    virtual vector<long long> _ComputeAB(
        const vector<int> &data, unsigned m, int r) override;
public:
    SampenCalculatorHG(double sample_rate_)
    : _sample_rate(sample_rate_) {}
};

// base class that compute A and B with vector of points
class ABCalculatorPoint
{
public:
    virtual vector<long long> ComputeAB(
        const vector<Point> &points, int r) = 0;
};

class ABCalculatorPointD : public ABCalculatorPoint
{
public:
    virtual vector<long long> ComputeAB(
        const vector<Point> &points, int r) override;
};

class ABCalculatorPointRT : public ABCalculatorPoint
{
public:
    virtual vector<long long> ComputeAB(
        const vector<Point> &points, int r) override;
};

class ABCalculatorDirectWeighted 
{
public:
    vector<double> ComputeAB(const vector<Point> &points, 
                             const vector<double> &weights, int r);   
};


double ComputeSampenDirect(
    const vector<int> &data, unsigned m, int r, double *a, double *b);

double ComputeSampenRangetree(
    const vector<int> &data, unsigned m, int r, double *a, double *b);

double ComputeSampenKdtree(
    const vector<int> &data, unsigned m, int r, double *a, double *b);

double ComputeSampenNkdtreeHist(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_size, unsigned sample_num, double *a, double *b);

double ComputeSampenKdtreeGrid(
    const vector<int> &data, unsigned m, int r, 
    double *a, double *b);

double ComputeSampenQR(
    const vector<int> &data, const unsigned m, const int r, 
    const unsigned sample_size, const unsigned sample_num, 
    double *a, double *b);

double ComputeSampenQR2( 
    const vector<int> &data, const unsigned m, const int r, 
    const unsigned sample_size, const unsigned sample_num, 
    double *a, double *b);

double ComputeSampenCoreset(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_size, unsigned sample_num, double *a, double *b);

double ComputeSampenUniform(
    const vector<int> &data, unsigned m, int r, 
    unsigned sample_size, unsigned sample_num, double *a, double *b);

double ComputeSampenRangetreeHist(
    const vector<int> &data, const unsigned m, 
    const int r, double sample_rate, double *a, double *b);
#endif // __SAMPEN_CALCULATOR_H__
