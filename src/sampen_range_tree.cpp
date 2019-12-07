/* file: main.cpp
 * date: 2019-11-05
 * author: phree
 *
 * description: compute sample entropy using range tree, implemented by
 * Luca Weihs
 */

#include <iostream>
#include <time.h>
#include <algorithm>

#include "RangeTree2.h"
#include "random_sampler.h"

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

long long count_range_points(const vector<RT::Point<int, int> > &points, 
                             const unsigned m, const int r);
long long count_range(const vector<int> &data, const unsigned m, const int r);
long long count_range_random(const vector<int> &data, 
                             const unsigned m, 
                             const int r, 
                             const unsigned sample_size);

double calculate_sampen_rangetree(
    const vector<int> &data, const unsigned m, const int r);
double calculate_sampen_rangetree_random(
    const vector<int> &data, const unsigned m, const int r, 
    const unsigned sample_size, const unsigned sample_num);
double calculate_sampen_rangetree_hist(
    const vector<int> &data, const unsigned m, const int r);

int *readdata(char *filenm, unsigned long *filelen);
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
        data, _stat.m, _stat.r);
    cout << "Range Tree (histogram): SampEn("; 
    cout << _stat.m << ", " << _stat.r << ", ";
    cout << N << ") = " << result_histogram << endl;

    error = (result_histogram - result) / result;
    cout << "Error (histogram) = " << error << endl;

    return 0;
}

long long count_range_points(const vector<RT::Point<int, int> > &points, 
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
	if (_stat.t) 
        cout << "time consumed in building tree: " << interval << "s\n";
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
	if (_stat.t) 
        cout << "time consumed in couting tree: " << interval << "s\n";
    return result;
}

long long count_range(const vector<int> &data, const unsigned m, const int r)
{
    if (data.size() <= 0)
        std::logic_error("data length is zero");
    if (m <= 0)
        std::logic_error("m is zero");
    if (data.size() < m + 1)
        std::logic_error("data length is too small (less than m + 1");

    // generate points
    vector<RT::Point<int, int> > points(data.size() - m + 1);
    for (vector<int>::size_type i = 0; i < data.size() - m + 1; i++)
    {
        vector<int> p(data.cbegin() + i, data.cbegin() + i + m);
        points[i] = RT::Point<int, int>(p, 0);
    }
    return count_range_points(points, m, r);
}

long long count_range_random(const vector<int> &data, 
                             const unsigned m, const int r, 
                             const unsigned sample_size)
{
    if (data.size() <= 0)
        std::logic_error("data length is zero");
    if (m <= 0)
        std::logic_error("m is zero");
    if (data.size() < m + 1)
        std::logic_error("data length is too small (less than m + 1");

    // generate points

    vector<RT::Point<int, int> > sampled_points(sample_size);
    uniform_int_generator uig(
        0, data.size() - 1, uniform_int_generator::QUASI);
    for (int j = 0; j < sample_size; j++)
    {
        int idx = uig.get();
        vector<int> p(data.cbegin() + idx, data.cbegin() + idx + m);
        sampled_points[j] = RT::Point<int, int>(p, 0);
    }
    return count_range_points(sampled_points, m, r);
}

vector<long long> _calculate_AB(const vector<int> &data, 
                                const unsigned m, const int r)
{
    struct timespec tstart, tend;
	double interval = 0;

    const auto N = data.size();

  	clock_gettime(CLOCK_REALTIME, &tstart);
    long long B = count_range(data, m + 1, r);
    B -= (N - m);
    vector<int> _data(data.cbegin(), data.cend() - 1);
    long long A = count_range(_data, m, r);
    A -= (N - m);
	clock_gettime(CLOCK_REALTIME, &tend);
	interval = (double)(tend.tv_sec - tstart.tv_sec);
	interval += (double)(tend.tv_nsec - tstart.tv_nsec) / 1e9;
    if (_stat.t) 
	    cout << "time consumed by range tree method: " << interval << "s\n";
    vector<long long> result(2);
    result[0] = A;
    result[1] = B;
    return result;
}

vector<long long> _calculate_AB_random(const vector<int> &data, 
                                       const unsigned m, 
                                       const int r, 
                                       const unsigned sample_size, 
                                       const unsigned sample_num)
{
    struct timespec tstart, tend;
	double interval = 0;

  	clock_gettime(CLOCK_REALTIME, &tstart);
    long long B = count_range_random(data, m + 1, r, sample_size);
    B -= (sample_size - m);
    vector<int> _data(data.cbegin(), data.cend() - 1);
    long long A = count_range_random(_data, m, r, sample_size);
    A -= (sample_size - m);

	clock_gettime(CLOCK_REALTIME, &tend);
	interval = (double)(tend.tv_sec - tstart.tv_sec);
	interval += (double)(tend.tv_nsec - tstart.tv_nsec) / 1e9;
	if (_stat.t) 
    cout << "time consumed by range tree method: " << interval << "s\n";

    vector<long long> result(2);
    result[0] = A;
    result[1] = B;
    return result;
}

vector<long long> _calculate_AB_hist(const vector<int> &data,
                                     const unsigned m,
                                     const int r,
                                     const double sample_rate)
{
    struct timespec tstart, tend;
    double interval = 0;

    const auto N = data.size();

    clock_gettime(CLOCK_REALTIME, &tstart);

    vector<Point> points = get_points(data, m + 1);
    int max_data = *std::max_element(data.cbegin(), data.cend());
    int min_data = *std::min_element(data.cbegin(), data.cend());
    vector<vector<Point> > vec_points = sample_hist(
        points, r, max_data, min_data, sample_rate);

#ifdef DEBUG
    for (unsigned i = 0; i < vec_points.size(); i++)
    {
        for (unsigned j = 0; j < vec_points[i].size(); j++)
        {
            vec_points[i][j].print();
            if (vec_points[i][j].dim() != vec_points[i][0].dim())
            throw std::logic_error("dimensions does not coincide.");
        }
    }
#endif
    vector<long long> As(vec_points.size());
    vector<long long> Bs(vec_points.size());

    for (unsigned i = 0; i < vec_points.size(); i++)
    {
        Bs[i] = count_range_points(vec_points[i], m+1, r);
        // WARN: Should I substract m?
        Bs[i] -= (vec_points[i].size() - m);
    }

    // Drop last dimension 
    for (unsigned i = 0; i < vec_points.size(); i++)
    {
        std::transform(vec_points[i].begin(), 
                       vec_points[i].end(), 
                       vec_points[i].begin(), 
                       [](const Point &p) -> Point 
                           { return p.drop_last(); });
    }

    for (unsigned i = 0; i < vec_points.size(); i++)
    {
        As[i] = count_range_points(vec_points[i], m, r);
        Bs[i] -= (vec_points[i].size() - m);
    }

    
    long long B = 0;
    std::for_each(Bs.begin(), Bs.end(), [&](long long n) { B += n; });
    long long A = 0;
    std::for_each(As.begin(), As.end(), [&](long long n) { A += n; });

    clock_gettime(CLOCK_REALTIME, &tend);
    interval = (double)(tend.tv_sec - tstart.tv_sec);
    interval += (double)(tend.tv_nsec - tstart.tv_nsec) / 1e9;
    if (_stat.t)
    {
        cout << "time consumed by range tree method (hist): "; 
        cout << interval << "s\n";
    }

    vector<long long> result(2);
    result[0] = A;
    result[1] = B;
    return result;
}

double calculate_sampen_rangetree(
    const vector<int> &data, const unsigned m, const int r)
{
    int N = data.size();
    if (N <= 0)
        throw std::logic_error("data length is zero");
    if (m <= 0)
        throw std::logic_error("m is zero");
    if (N < m + 1)
        throw std::logic_error("data length is too small (less than m + 1");

    vector<long long> AB = _calculate_AB(data, m, r);
    long long A = AB[0];
    long long B = AB[1];
    if (A * B) return -log(static_cast<double>(B) / A);
    else return log(static_cast<double>(N -m - 1) / (N - m));
}


double calculate_sampen_rangetree_random(
    const vector<int> &data, const unsigned m, const int r, 
    const unsigned sample_size, const unsigned sample_num) 
{
    int N = data.size();
    if (N <= 0)
        std::logic_error("data length is zero");
    if (m <= 0)
        std::logic_error("m is zero");
    if (N < m + 1)
        std::logic_error("data length is too small (less than m + 1");

    double A = 0., B = 0.;
    for (int i = 0; i < sample_num; i++)
    {
        vector<long long> AB = _calculate_AB_random(
            data, m, r, sample_size, sample_num);
        A += static_cast<double>(AB[0]);
        B += static_cast<double>(AB[1]);
    }
    A /= sample_num;
    B /= sample_num;

    if (A * B) return -log(B / A);
    else return log((N -m - 1) / (N - m));
}

double calculate_sampen_rangetree_hist(
    const vector<int> &data, const unsigned m, const int r)
{
    int N = data.size();
    if (N <= 0)
        throw std::logic_error("data length is zero");
    if (m <= 0)
        throw std::logic_error("m is zero");
    if (N < m + 1)
        throw std::logic_error("data length is too small (less than m + 1");

    double sample_rate = 0.1;
    vector<long long> AB = _calculate_AB_hist(data, m, r, sample_rate);
    long long A = AB[0];
    long long B = AB[1];
    if (A * B) return -log(static_cast<double>(B) / A);
    else return log(static_cast<double>(N -m - 1) / (N - m));
}

int *readdata(char *filenm, unsigned long *filelen)
{
	FILE *ifile;
	unsigned long maxdat = 0L, npts = 0L, i;
	int *data = NULL, y;
	int maxval = INT_MIN, minval = INT_MAX;
	if (strcmp(filenm, "-") == 0)
	{
		filenm = "standard input";
		ifile = stdin;
	}
	else if ((ifile = fopen(filenm, "rt")) == NULL)
	{
		cerr << "could not open file " << filenm << endl;;
		exit(1);
	}
	while (fscanf(ifile, "%d", &y) == 1)
	{
		if (++npts >= maxdat)
		{
			int *s;

			maxdat += 4096; /* allow the input buffer to grow (the
				   increment is arbitrary) */
			if ((s = (int *) realloc(data, maxdat * sizeof(int))) == NULL)
			{
				cerr << "insufficient memory, truncating input at row ";
                cerr << maxdat << endl;
				break;
			}
			data = s;
		}
		data[npts - 1] = y;
		if (y > maxval) maxval = y;
		if (y < minval) minval = y;
	}

	fclose(ifile);

	if (npts < 1)
	{
        cerr << "file " << filenm << " contains no data\n";
		exit(-1);
	}

	*filelen = npts;
	for (i = 0; i < *filelen; i++) data[i] -= minval;
	return data;
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
