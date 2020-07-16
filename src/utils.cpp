#include <iostream>
#include <algorithm>

#include <string.h>
#include <limits.h>
#include "utils.h"

using std::cerr;
using std::endl;

vector<Point> GetPoints(const vector<int> &data, unsigned m)
{
    vector<Point> result(data.size() - m + 1);
    for (unsigned i = 0; i < result.size(); i++)
    {
        Point p(vector<int>(data.begin() + i, data.begin() + i + m), 0);
        result[i] = p;
    }
    return result;
}

bool IsPowerTwo(unsigned n)
{
	if (n == 1) return true;
	else if (n % 2 == 0) return IsPowerTwo(n / 2);
	else return false;
}

double ComputeVarience(const vector<int> &data) 
{
	vector<double> data_(data.cbegin(), data.cend());
	double avg = ComputeSum(data_) / data.size();
	std::for_each(data_.begin(), data_.end(), [avg] (double &x) 
		{ x -= avg; x *= x; });
	double sum = ComputeSum(data_); 
	sum /= data.size();
	return sum;
}

double ComputeSum(const vector<double> &data)
{
	unsigned n0 = 1 << 10;
	unsigned p = data.size() / n0;
	double sum = 0.;
	sum = std::accumulate(data.cbegin() + p * n0, data.cend(), 0.);
	if (p == 0) return sum;
	else 
	{
		vector<double> temp_sum(p, 0);
		for (unsigned i = 0; i < p; i++) 
		{
			temp_sum[i] = std::accumulate(
				data.cbegin() + i * n0, data.cbegin() + (i + 1) * n0, 0.);
		}
		temp_sum.push_back(sum);
		return ComputeSum(temp_sum);
	}
}

double L1Distance(const Point &p1, const Point &p2) 
{
	double max_diff = 0; 
	for (unsigned i = 0; i < p1.dim(); i++) 
	{
		double diff = abs(p1[i] - p2[i]);
		if (diff > max_diff) max_diff = diff;
	}
	return max_diff;
}

double EclideanDistance(const Point &p1, const Point &p2)
{
	double sum = 0.;
	for (unsigned i = 0; i < p1.dim(); i++)
	{
		double diff = p1[i] - p2[i];
		sum += diff * diff;
	}
	return sqrt(sum);
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

string ArgumentParser::getArg(const string &arg) 
{
	auto iter = std::find(arg_list.cbegin(), arg_list.cend(), arg);
	if (iter != arg_list.cend() && ++iter != arg_list.cend()) 
	{
		return *iter;
	}
	else return string("");
}

bool ArgumentParser::isOption(const string &opt) 
{
	return find(arg_list.cbegin(), arg_list.cend(), opt) != arg_list.cend();
}
