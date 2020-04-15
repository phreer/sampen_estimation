/* file: utils.h
 * date: 2019-12-07
 * author: phree
 *
 * description: utility functions and definitions
 */
#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <string>

#include "RangeTree2.h"

using std::vector;
using std::string;

namespace RT = RangeTree;

typedef RT::Point<int, int> Point;

int *readdata(char *filenm, unsigned long *filelen);

vector<Point> get_points(const vector<int> &data, unsigned m);

bool IsPowerTwo(unsigned n);
double ComputeVarience(const vector<int> &data);
double ComputeSum(const vector<double> &data);
double EclideanDistance(const Point &p1, const Point &p2);
double L1Distance(const Point &p1, const Point &p2);

class ArgumentParser
{
public:
    ArgumentParser(int argc, char *argv[]) : arg_list(argv, argv + argc) {}
    string getArg(const string &arg);
    bool isOption(const string &opt);
private:
    vector<string> arg_list;
};

#endif // __UTILS_H__


