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

double compute_sampen(long long A, long long B, unsigned N, unsigned m)
{
    if (A * B)
        return -log(static_cast<double>(B) / A);
    else
        return -log(static_cast<double>(N - m - 1) / (N - m));
}

class sampen_calculator
{
private:
    void check_dim(const vector<int> &data, unsigned m)
    {
        if (data.size() <= m)
            throw std::logic_error("data.size() < m");
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

        long long A = AB[0];
        long long B = AB[1];
#ifdef DEBUG
        std::cout << "A is " << A << ", B is " << B << std::endl;
#endif
        unsigned N = data.size();

        return compute_sampen(A, B, N, m);
    }
};

class sampen_calculator_d : public sampen_calculator
{
private:
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
};

class sampen_calculator_rt : public sampen_calculator
{
private:
    virtual vector<long long> _compute_AB(
        const vector<int> &data, unsigned m, int r) override;
};

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

/* 
 * Compute sampen by sample using histogram
 */
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

// Classes that compute A and B with vector of points
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

#endif // __SAMPEN_CALCULATOR_H__
