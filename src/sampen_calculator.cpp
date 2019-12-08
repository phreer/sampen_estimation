/* file: entropy.cpp
 * date: 2019-12-07
 * author: phree
 *
 * description: implementation of sampen
 */

#include <algorithm>
#include "sampen_calculator.h"
#include "random_sampler.h"

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
        0, data.size() - 1, uniform_int_generator::QUASI);
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

vector<long long> sampen_calculator_hg::_compute_AB(
    const vector<int> &data, unsigned m, int r) 
{
    AB_calculator_point_d ABc;
    struct timespec tstart, tend;
    double interval = 0;

    const auto N = data.size();

    clock_gettime(CLOCK_REALTIME, &tstart);

    vector<Point> points = get_points(data, m + 1);
    int max_data = *std::max_element(data.cbegin(), data.cend());
    int min_data = *std::min_element(data.cbegin(), data.cend());
    vector<vector<Point> > vec_points = sample_hist(
        points, r, max_data, min_data, _sample_rate);

    vector<long long> ABs(2, 0);

    for (unsigned i = 0; i < vec_points.size(); i++)
    {
        vector<long long> AB(2);
        AB = ABc.compute_AB(vec_points[i], r);
        ABs[0] += AB[0];
        ABs[1] += AB[1];
    }

    clock_gettime(CLOCK_REALTIME, &tend);
    interval = (double)(tend.tv_sec - tstart.tv_sec);
    interval += (double)(tend.tv_nsec - tstart.tv_nsec) / 1e9;

    return ABs;
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
    unsigned N = points.size() + m - 1;
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
	double interval = 0, total = 0;

    /* build tree */
  	clock_gettime(CLOCK_REALTIME, &tstart);
    RT::RangeTree<int, int> rtree(points);
    vector<int> lower(m, 0), upper(m, 0);
	clock_gettime(CLOCK_REALTIME, &tend);
	interval = (double)(tend.tv_sec - tstart.tv_sec);
	interval += (double)(tend.tv_nsec - tstart.tv_nsec) / 1e9;
	total += interval;

    /* counting */
	clock_gettime(CLOCK_REALTIME, &tstart);
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
  	clock_gettime(CLOCK_REALTIME, &tend);
	interval = (double)(tend.tv_sec - tstart.tv_sec);
	interval += (double)(tend.tv_nsec - tstart.tv_nsec) / 1e9;
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