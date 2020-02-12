#include <iostream>
#include <algorithm>

#include <string.h>
#include <limits.h>
#include "utils.h"

using std::cerr;
using std::endl;

vector<Point> get_points(const vector<int> &data, unsigned m)
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
	double avg = std::accumulate(data.cbegin(), data.cend(), 0) / data.size();
	double sum = 0; 
	std::for_each(data.cbegin(), data.cend(), [&sum, avg] (const int x) {
		sum += (x - avg) * (x - avg);
	});
	sum /= data.size();
	return sum;
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
