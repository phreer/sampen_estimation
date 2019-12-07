/* file: utils.h
 * date: 2019-12-07
 * author: phree
 *
 * description: utility functions and definitions
 */
#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include "RangeTree2.h"

using std::vector;
namespace RT = RangeTree;

typedef RT::Point<int, int> Point;

int *readdata(char *filenm, unsigned long *filelen);

vector<Point> get_points(const vector<int> &data, unsigned m);

#endif // __UTILS_H__


