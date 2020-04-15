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

double compute_sampen(double A, double B, unsigned N, unsigned m)
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
class sampen_calculator
{
private:
    void check_dim(const vector<int> &data, unsigned m)
    {
        if (data.size() <= m)
            throw std::invalid_argument("data.size() < m");
    }
    virtual vector<long long> _compute_AB(const vector<int> &data,
                                          unsigned m, int r) = 0;

public:
    vector<long long> compute_AB(const vector<int> &data,
                                 unsigned m, int r)
    {
        check_dim(data, m);
        return _compute_AB(data, m, r);
    }
    virtual double compute_entropy(const vector<int> &data,
                                   unsigned m, int r)
    {
        auto start = std::chrono::system_clock::now();
        vector<long long> AB = compute_AB(data, m, r);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> interval = end - start;
        std::cout << "time: " << interval.count() << "s" << std::endl;
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
        result = compute_sampen(A, B, N, m);
#ifdef DEBUG
        std::cout << "A is " << A << ", B is " << B << std::endl;
#endif
        // for (unsigned i = 0; i < sample_num; i++) 
        //     result += compute_sampen(AB[i * 2], AB[i * 2 + 1], N, m);
        // result /= sample_num;

        return result;
    }
};

// direct method
class sampen_calculator_d : public sampen_calculator
{
private:
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
};

// range tree
class sampen_calculator_rt : public sampen_calculator
{
private:
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
};

// old kd tree
class sampen_calculator_kd: public sampen_calculator
{
private:
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
};

// Compute sample entropy by kd tree divided according to grid
class sampen_calculator_kdg: public sampen_calculator
{
private:
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
};

// quasi-random sampling
class sampen_calculator_qr : public sampen_calculator
{
public:
    explicit sampen_calculator_qr(
        unsigned sample_num_, unsigned sample_size_, bool random = false)
        : sample_num(sample_num_), sample_size(sample_size_), random(random) 
    {}
    void set_sample_num(unsigned sample_num_) { sample_num = sample_num_; }
    void set_sample_size(unsigned sample_size_)
    {
        sample_size = sample_size_;
    }

private:
    virtual vector<long long> _compute_AB(const vector<int> &data,
                                          unsigned m, int r) override;
    unsigned sample_num;
    unsigned sample_size;
    bool random;
};

// quasi-random sampling
class sampen_calculator_qr2 : public sampen_calculator
{
public:
    explicit sampen_calculator_qr2(
        unsigned sample_num_, unsigned sample_size_, bool random = false)
        : sample_num(sample_num_), sample_size(sample_size_), random(random) 
    {}
    void set_sample_num(unsigned sample_num_) { sample_num = sample_num_; }
    void set_sample_size(unsigned sample_size_)
    {
        sample_size = sample_size_;
    }

private:
    virtual vector<long long> _compute_AB(const vector<int> &data,
                                          unsigned m, int r) override;
    unsigned sample_num;
    unsigned sample_size;
    bool random;
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
    double ComputeSampen(const vector<int> &data, unsigned m, int r);
private:
    vector<double> _ComputeAB(const vector<int> &data,
                                 unsigned m, int r);
    unsigned sample_num;
    unsigned sample_size;
    bool random;
};

// Compute sample entropy using new kd tree
// Here, we require sample_size = 2^n for some non-negative integer n
class sampen_calculator_nkd : public sampen_calculator
{
public:
    explicit sampen_calculator_nkd(
        unsigned sample_num, unsigned sample_size)
        : sample_num_(sample_num), sample_size_(sample_size) 
        {
            if (!IsPowerTwo(sample_size_)) 
                throw std::invalid_argument("sample_size should be power of 2");
        }
private:
    virtual vector<long long> _compute_AB(const vector<int> &data,
                                          unsigned m, int r) override;
    unsigned sample_num_;
    unsigned sample_size_;
};

// Compute sampen by sample using histogram
class sampen_calculator_hg : public sampen_calculator
{
private:
    double _sample_rate;
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
public:
    sampen_calculator_hg(double sample_rate_)
    : _sample_rate(sample_rate_) {}
};

// base class that compute A and B with vector of points
class AB_calculator_point
{
public:
    virtual vector<long long> compute_AB(
        const vector<Point> &points, int r) = 0;
};

class AB_calculator_point_d : public AB_calculator_point
{
public:
    virtual vector<long long> compute_AB(
        const vector<Point> &points, int r) override;
};

class AB_calculator_point_rt : public AB_calculator_point
{
public:
    virtual vector<long long> compute_AB(
        const vector<Point> &points, int r) override;
};

class ABCalculatorDirectWeighted 
{
public:
    vector<double> ComputeAB(const vector<Point> &points, 
                             const vector<double> &weights, int r);   
};
#endif // __SAMPEN_CALCULATOR_H__
