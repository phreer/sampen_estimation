/* file: entropy_calculator.h
 * date: 2020-03-03
 * author: phree
 *
 * description: This file contains lasses to calculate entropy. 
 */


#ifndef __ENTROPY_CALCULATOR_H__
#define __ENTROPY_CALCULATOR_H__

#include <vector>

using std::vector;

/* This class used to compute DynamicEntropy, which is based on blocks 
 * genreated by kd-tree. 
 */ 
class DynamicEntropyCalculator 
{
public:
    double ComputeEntropy(const vector<int> &data, unsigned m, unsigned max_level);
};

#endif // __ENTROPY_CALCULATOR_H__

