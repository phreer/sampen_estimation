/* file: tensor.h
 * date: 2019-11-25
 * author: phree
 *
 * description: class tensor
 */

#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <vector>

template <typename T>
class tensor
{
public:

    explicit tensor(vector<unsigned> _size): size(_size) {
        unsigned len = 1;
        for (unsigned i = 0; i < size.size(); i++)
        len *= size[i];
        data = vector<T>(len);
    }

    T& operator[] (const vector<unsigned> &idx)
    {
        return data[get_idx(idx)];
    }

    const T& operator[] (const vector<unsigned> &idx) const 
    {
        return data[get_idx(idx)];
    }

    T& operator[] (unsigned idx)
    {
        return data[idx];
    }

    const T& operator[] (unsigned idx) const 
    {
        return data[idx];
    }

private:
    std::vector<T> data;
    std::vector<unsigned> size;
    unsigned get_idx(const vector<unsigned> &idx) 
    {
//#ifdef DEBUG
//        if (idx.size() != size.size())
//        throw std::logic_error("idx.size() != size.size()");
//#endif
        unsigned uidx = 0;
        for (unsigned i = 0; i < idx.size(); i++)
        {
            uidx *= size[i];
            uidx += idx[i];
        }
        return uidx;
    }
};

#endif // __TENSOR_H__
