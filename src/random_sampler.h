/* file: random_sampler.h
 * date: 2019-11-25
 * author: phree
 *
 * description: class used to generate random sequence
 */

#ifndef __RANDOM_SAMPLER_H__
#define __RANDOM_SAMPLER_H__

#include <vector>
#include <random>
#include <gsl/gsl_qrng.h>
#include "RangeTree2.h"

using std::vector;

namespace RT = RangeTree;

typedef RT::Point<int, int> Point;

class uniform_int_generator
{
public:
    enum random_type {PSEUDO, QUASI, SHUFFLE};
    uniform_int_generator(int _rangel, int _ranger, random_type _rtype): 
            rangel(_rangel), ranger(_ranger), rtype(_rtype)
    {
        init_state();
    }
    int get();
private:
    int rangel, ranger;
    const uniform_int_generator::random_type rtype;
    std::uniform_int_distribution<int> uid;
    std::default_random_engine eng;
    gsl_qrng *qrng;
    int sample;

    void init_state();
};

/*
 * Convert a sequence of data to vector of Points
 */
vector<Point> get_points(const vector<int> &data, unsigned m);

/* 
 * Sample a vector of Points by sampling according to histogram
 * 
 * @param vec: the vector of Points
 * @param r: the width of grid of the histogram
 * @param max_data: the maximum of the original data
 * @param min_data: the minimum of the original data
 * @param sample_rate: the sampling rate
 * @return a vector of vectors of Points
 */
vector<vector<Point> > sample_hist(const vector<Point> &vec, int r, 
                                  int max_data, int min_data, 
                                  double sample_rate);

#endif // __RANDOM_SAMPLER_H__
