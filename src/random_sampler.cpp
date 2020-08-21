#include <iostream>
#include <chrono>
#include <algorithm>
#include <list>

#include <math.h>
#include <stdlib.h>
#include <random>

#include <gsl/gsl_qrng.h>

#include "random_sampler.h"
#include "tensor.h"

using std::vector;
using std::list;

unsigned powu(unsigned x, unsigned u) 
{
    unsigned result = x;
    for (unsigned i = 1; i < u; i++)
    {
        result *= x;
    }
    return result;
}

// random generator function
int randomfunc(int j)
{
    return rand() % j;
}

vector<unsigned> random_permutation(unsigned n)
{
    vector<unsigned> result(n);
    for (unsigned i = 0; i < n; i++)
    result[i] = i;
    std::random_shuffle(result.begin(), result.end(), randomfunc);
    return result;
}

void uniform_int_generator::init_state() 
{
    if (rangel > ranger) 
    throw std::invalid_argument("rangel is larger than ranger.");

    if (real_random) 
    {
        std::random_device rd;
        eng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }
    else 
    {
        eng.seed(0); 
    }
    if (rtype == PSEUDO) 
    {
        uid = std::uniform_int_distribution<int>(rangel, ranger);
    }
    else if (rtype == QUASI) 
    {
        qrng = gsl_qrng_alloc(gsl_qrng_sobol, 1);
    }
    else if (rtype == SHUFFLE)
    {
        throw std::runtime_error("SHUFFULE not implemented.");
    }
}


int uniform_int_generator::get() 
{
    switch (rtype)
    {
    case PSEUDO:
        sample = uid(eng);
        break;
    case QUASI:
        double v;
        gsl_qrng_get(qrng, &v);
        sample = static_cast<int>(v * (ranger - rangel) + rangel);
        break;
    case SHUFFLE:
        throw std::runtime_error("SHUFFULE not implemented.");
        break;
    }
    return sample;
}


vector<vector<Point> > sample_hist(const vector<Point> &vec, int r, 
                                   int max_data, int min_data, 
                                   double sample_rate)
{
    // Construct histogram
    unsigned num_grid = (max_data - min_data) / r + 1;
    unsigned dim = vec[0].dim();

    vector<unsigned> size(dim, num_grid);

    tensor<list<const Point *> > hist(size);
    for (unsigned i = 0; i < vec.size(); i++)
    {
        vector<unsigned> idx(dim);
        for (unsigned j = 0; j < dim; j++)
        {
            idx[j] = (vec[i][j] - min_data) / r;
        }
        hist[idx].push_front(vec.data() + i);
    }

    unsigned num_sample = static_cast<unsigned>(1 / sample_rate) + 1;
    unsigned size_sample = static_cast<unsigned>(vec.size() * sample_rate) * 3;
    vector<vector<Point> > results(num_sample, vector<Point>(size_sample));
    vector<unsigned> counts(num_sample, 0);

    // Do sampling!
    for (unsigned i = 0; i < powu(num_grid, size.size()); i++)
    {
        unsigned len = hist[i].size();
        vector<unsigned> perm = random_permutation(
            (len / size_sample + 1) * size_sample);

        for (unsigned j = 0; j < len; j++)
        {
            unsigned idx = perm[j] % num_sample;
            results[idx][counts[idx]] = *hist[i].front();
#ifdef DEBUG
            std::cout << "counts[" << idx << "] = " << counts[idx] << std::endl;
            hist[i].front()->print();
#endif
            hist[i].pop_front();
            counts[idx]++;
        }
    }

    for (unsigned i = 0; i < results.size(); i++)
    {
        results[i].resize(counts[i]);

#ifdef DEBUG
        std::cout << counts[i] << std::endl;
        std::for_each(results[i].begin(), results[i].end(), 
                      [](const Point &p) { p.print(); });
#endif
    }
    return results;
}
