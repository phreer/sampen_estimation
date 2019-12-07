#include <iostream>
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