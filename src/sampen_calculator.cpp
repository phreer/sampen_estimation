/* file: entropy.cpp
 * date: 2019-12-07
 * author: phree
 *
 * description: implementation of sampen
 */
#include <iostream>

#include <algorithm>
#include <chrono>
#include <string.h>
#include <math.h>
#include <utility>

#include "sampen_calculator.h"
#include "random_sampler.h"
#include "kdtree.h"
#include "utils.h"

using std::pair;

vector<long long> sampen_calculator_d::_compute_AB(
    const vector<int> &data, unsigned m, int r)
{
    AB_calculator_point_d ABc;
    vector<Point> points = get_points(data, m + 1);
    return ABc.compute_AB(points, r);
}

vector<long long> sampen_calculator_rt::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    AB_calculator_point_rt ABc;
    vector<Point> points = get_points(data, m + 1);
    return ABc.compute_AB(points, r);    
}

vector<long long> sampen_calculator_qr::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    vector<Point> points = get_points(data, m + 1);
    std::sort(points.begin(), points.end(), 
        [](const Point &p1, const Point &p2) 
        {
            for (unsigned i = 0; i < p1.dim(); i++) 
            {
                if (p1[i] > p2[i]) return true;
                if (p1[i] < p2[i]) return false;
            }
            return true;
        });
    vector<Point> sampled_points(sample_size);
    
    uniform_int_generator uig(
        0, points.size()-1, uniform_int_generator::QUASI, random);
    AB_calculator_point_d ABc;

    vector<long long> AB(2);
    vector<long long> ABs(2, 0);
    for (int i = 0; i < sample_num; i++)
    {
        for (int j = 0; j < sample_size; j++)
        {
            // generate points
            int idx = uig.get();
            // vector<int> p(data.cbegin() + idx, data.cbegin() + idx + m + 1);
            sampled_points[j] = points[idx];
        }
        AB = ABc.compute_AB(sampled_points, r);
        ABs[0] += AB[0];
        ABs[1] += AB[1];
    }
    return ABs;
}

vector<long long> sampen_calculator_qr2::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    vector<Point> points = get_points(data, m + 1);
    vector<Point> sampled_points(sample_size);
    std::sort(points.begin(), points.end(), 
        [](const Point &p1, const Point &p2) 
        {
            for (unsigned i = 0; i < p1.dim(); i++) 
            {
                if (p1[i] > p2[i]) return true;
                if (p1[i] < p2[i]) return false;
            }
            return true;
        });
    
    AB_calculator_point_d ABc;

    vector<long long> AB(2);
    vector<long long> ABs(2, 0);

    uniform_int_generator uig(
        0, points.size()-1, uniform_int_generator::QUASI, random);
    unsigned iter_num = 8;
    // double decay_rate = 0.8;
    for (unsigned i = 0; i < iter_num; i++) 
    {
        for (unsigned j = 0; j < 9182; j++)
        {
            int idx = uig.get();
            Point p = points[idx];
            for (unsigned k = 0; k < points.size(); k++)
            {
                if (points[k].within(p, m+1, r)) points[k].increaseCountByOne();
            }
        }
        std::sort(points.begin(), points.end(), 
            [] (const Point &p1, const Point &p2) {
                return (p1.count() > p2.count());
            });
    }

    for (int i = 0; i < sample_num; i++)
    {
        for (int j = 0; j < sample_size; j++)
        {
            // generate points
            int idx = uig.get();
            // vector<int> p(data.cbegin() + idx, data.cbegin() + idx + m + 1);
            sampled_points[j] = points[idx];
        }
        AB = ABc.compute_AB(sampled_points, r);
        ABs[0] += AB[0];
        ABs[1] += AB[1];
    }
    return ABs;
}

vector<long long> sampen_calculator_nkd::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    uniform_int_generator uig(
        0, data.size() - 1, uniform_int_generator::PSEUDO, random);
    AB_calculator_point_d ABc;

    vector<long long> AB(2);
    vector<long long> ABs(2 * sample_num_, 0);
    
    auto max_level = static_cast<unsigned>(log2(sample_size_));
    vector<Point> points = get_points(data, m + 1);
    NewKDTree kdtree(points, max_level);
    for (int i = 0; i < sample_num_; i++)
    {
        auto sampled_points = kdtree.Sample(sample_size_);
        AB = ABc.compute_AB(sampled_points, r);
        ABs[i * 2] = AB[0];
        ABs[i * 2 + 1] = AB[1];
    }
    return ABs;
}

vector<long long> sampen_calculator_hg::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    AB_calculator_point_d ABc;

    vector<Point> points = get_points(data, m + 1);
    int max_data = *std::max_element(data.cbegin(), data.cend());
    int min_data = *std::min_element(data.cbegin(), data.cend());

    auto start = std::chrono::system_clock::now();
    vector<vector<Point> > vec_points = sample_hist(
        points, r, max_data, min_data, _sample_rate);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> interval = end - start;
    std::cout << "Time consumed in sample_hist: ";
    std::cout << interval.count() << "s" << std::endl;
    vector<long long> ABs(2, 0);

    for (unsigned i = 0; i < vec_points.size(); i++)
    {
        vector<long long> AB(2);
        AB = ABc.compute_AB(vec_points[i], r);
        ABs[0] += AB[0];
        ABs[1] += AB[1];
    }

    return ABs;
}

vector<long long> sampen_calculator_kd::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    long long A = 0, B = 0;
    unsigned N = data.size();
	unsigned n = N - m + 1;

    vector<int> __data = vector<int>(data);

    int **datap = (int **)malloc(n * sizeof(int *));

    for (unsigned i = 0; i < n; i++)
        datap[i] = __data.data() + i;

	struct kdtree *treem = build_kdtree((const int **)datap, n - 1, m, 0, 0);
    struct kdtree *treem1 = build_kdtree((const int **)datap, n - 1, m + 1, 0, 0);

    for (unsigned i = 0; i < N - m; i++)
    {
        A += count_range_kdtree(treem, __data.data() + i, m, r);
    }

    for (unsigned i = 0; i < N - m; i++)
    {
        B += count_range_kdtree(treem1, __data.data() + i, m + 1, r);
    }

    A -= (N - m);
    B -= (N - m);

    vector<long long> result(2);
    result[0] = A;
    result[1] = B;
    return result;
}

Point ComputeMean(const vector<Point> &points)
{
    unsigned n = points.size();
    if (n == 0) throw std::invalid_argument("points.size() == 0");
    Point result(points[0]);
    for (unsigned i = 0; i < result.dim(); i++)
    {
        vector<double> data(n);
        for (unsigned j = 0; j < n; j++) 
        {
            data[j] = points[j][i];
        }
        result[i] = ComputeSum(data) / n;
    }
    return result;
}

unsigned GetInvertalIndex(const vector<double> &pmf, double x) 
{
    if (pmf.size() == 2) return 0;
    unsigned median = pmf.size() / 2;
    if (x <= pmf[median]) 
    {
        vector<double> pmf_(pmf.cbegin(), pmf.cbegin() + median + 1);
        return GetInvertalIndex(pmf_, x);
    }
    else 
    {
        vector<double> pmf_(pmf.cbegin() + median, pmf.cend());
        return GetInvertalIndex(pmf_, x) + median;
    }
}

pair<vector<vector<Point> >, vector<vector<double> > >
SampleCoreset(const vector<Point> &points, 
    unsigned sample_size, unsigned sample_num) 
{
    unsigned n = points.size();
    // Generate PMF
    Point mean = ComputeMean(points);

    vector<double> distances(n);
    for (unsigned i = 0; i < n; i++) 
        distances[i] = L1Distance(points[i], mean);
    double sum = ComputeSum(distances);
    
    vector<double> q(n);
    for (unsigned i = 0; i < n; i++) 
    {
        q[i] = 1. / 2 / n + distances[i] / 2 / sum;
    }
    vector<double> pmf(n + 1, 0);
    for (unsigned i = 1; i < n + 1; i++) 
    {
        pmf[i] = q[i - 1] + pmf[i -  1];
        // std::cout << "pmf: " << pmf[i] << std::endl;
    }    

    // Do sampling.
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    vector<vector<Point> > result(sample_num);
    vector<vector<double> > weights(sample_num);
    for (unsigned i = 0; i < sample_num; i++)
    {
        result[i] = vector<Point>(sample_size);
        weights[i] = vector<double>(sample_size);
        for (unsigned j = 0; j < sample_size; j++)
        {
            double sample = dist(e2);
            unsigned index = GetInvertalIndex(pmf, sample);
            result[i][j] = points[index];
            weights[i][j] = 1. / q[index] / n;
            // std::cout << "weights: " << weights[i][j] << std::endl;
        }
    }
    return std::make_pair(result, weights);
}

vector<double> SampenCalculatorCoreset::_ComputeAB(
    const vector<int> &data, unsigned m, int r) 
{
    vector<Point> points = get_points(data, m + 1);
    auto sampled = SampleCoreset(points, sample_size, sample_num);
    auto sampled_points = sampled.first;
    auto weights = sampled.second;
    
    vector<double> result(2, 0);
    ABCalculatorDirectWeighted ABc;
    for (unsigned i = 0; i < sample_num; i++) 
    {
        auto AB = ABc.ComputeAB(sampled_points[i], weights[i], r);
        result[0] += AB[0]; 
        result[1] += AB[1];
    }
    return result;
}

double SampenCalculatorCoreset::ComputeSampen(
    const vector<int> &data, unsigned m, int r)
{
    auto start = std::chrono::system_clock::now();
    auto AB = _ComputeAB(data, m, r);
    double sampen = compute_sampen(AB[0], AB[1], data.size(), m);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> interval = end - start;
    std::cout << "time: " << interval.count() << "s" << std::endl;
    return sampen;
}

vector<long long> sampen_calculator_kdg::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{

    long long A = 0, B = 0;
    unsigned N = data.size();

    int max_ = *std::max_element(data.cbegin(), data.cend());
    int min_ = *std::min_element(data.cbegin(), data.cend());
    double diff = static_cast<double>(max_ - min_);
    unsigned p = static_cast<unsigned>(ceil(log2(diff)));
    struct kdtree *treem = build_kdtree_grid(data.data(), N - 1, m, p);
    struct kdtree *treem1 = build_kdtree_grid(data.data(), N, m + 1, p);

    for (unsigned i = 0; i < N - m; i++)
    {
        A += count_range_kdtree(treem, data.data() + i, m, r);
    }

    for (unsigned i = 0; i < N - m; i++)
    {
        B += count_range_kdtree(treem1, data.data() + i, m + 1, r);
    }

    A -= (N - m);
    B -= (N - m);

    vector<long long> result(2);
    result[0] = A;
    result[1] = B;
    return result;
}


// compute A and B with points using direct method
// TODO: this part can be parallelized
vector<long long> AB_calculator_point_d::compute_AB(
    const vector<Point> &points, int r)
{
    vector<long long> result(2, 0);
    if (points.size() == 0) return result;

    long long A = 0;
    long long B = 0;
    unsigned m = points[0].dim() - 1;
    for (unsigned i = 0; i < points.size(); i++)
    {
        for (unsigned j = i + 1; j < points.size(); j++)
        {
            if (points[i].within(points[j], m, r)) 
            {
                A += 1;
                int diff = points[j][m] - points[i][m];
                if ((-r <= diff) && (diff <= r))
                    B += 1;
            }
        }
    }
    result[0] = A;
    result[1] = B;
    return result;
}


long long count_range_point_rt(const vector<Point> &points, 
                               const unsigned m, const int r)
{
	struct timespec tstart, tend;

    /* build tree */
    RT::RangeTree<int, int> rtree(points);
    vector<int> lower(m, 0), upper(m, 0);

    /* counting */
    long long result = 0;
    for (vector<int>::size_type i = 0; i < points.size(); i++)
    {
        for (unsigned j = 0; j < m; j++)
        {
            lower[j] = points[i][j] - r;
            upper[j] = points[i][j] + r;
        }
        result += rtree.countInRange(lower, upper);
    }
    return result;
}

vector<long long> AB_calculator_point_rt::compute_AB(
    const vector<Point> &points, int r)
{
    vector<long long> result(2, 0);
    if (points.size() == 0) return result;

    long long A = 0;
    long long B = 0;
    unsigned m = points[0].dim() - 1;
    unsigned N = points.size() + m;
    B = count_range_point_rt(points, m + 1, r);
    B -= (N - m);
    vector<Point> _points(points.size());
    std::transform(points.begin(), points.end(), _points.begin(), 
                   [](const Point &p) -> Point { return p.drop_last(); });
    A = count_range_point_rt(_points, m, r);
    A -= (N - m);
    result[0] = A;
    result[1] = B;
    return result;
}

vector<double> ABCalculatorDirectWeighted::ComputeAB(
    const vector<Point> &points, const vector<double> &weights, int r)
{
    vector<double> result(2, 0);
    if (points.size() == 0) return result;

    double A = 0;
    double B = 0;
    unsigned m = points[0].dim() - 1;
    for (unsigned i = 0; i < points.size(); i++)
    {
        for (unsigned j = i + 1; j < points.size(); j++)
        {
            if (points[i].within(points[j], m, r)) 
            {
                double w = weights[i] * weights[j];
                A += w;
                int diff = points[j][m] - points[i][m];
                if ((-r <= diff) && (diff <= r))
                    B += w;
            }
        }
    }
    result[0] = A;
    result[1] = B;
    return result;
}