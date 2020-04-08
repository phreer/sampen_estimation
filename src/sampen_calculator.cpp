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

#include "sampen_calculator.h"
#include "random_sampler.h"
#include "kdtree.h"


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
    vector<Point> sampled_points(sample_size);
    uniform_int_generator uig(
        0, data.size() - 1, uniform_int_generator::QUASI, random);
    AB_calculator_point_d ABc;

    vector<long long> AB(2);
    vector<long long> ABs(2, 0);
    for (int i = 0; i < sample_num; i++)
    {
        for (int j = 0; j < sample_size; j++)
        {
            // generate points
            int idx = uig.get();
            vector<int> p(data.cbegin() + idx, data.cbegin() + idx + m + 1);
            sampled_points[j] = RT::Point<int, int>(p, 0);
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
