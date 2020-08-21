/* file: kdtree.h
 * date: 2019-12-10
 * author: phree
 *
 * description: data structure and functions for kd tree. Besides, a new 
 *   kd tree written on 2020-2-10 using cpp is provided, which stores points 
 *   in its nodes.
 */
#ifndef __KDTREE_H__
#define __KDTREE_H__

#include <vector>
#include <memory>

#include "utils.h"

using std::vector;
using std::shared_ptr;

class KDTreeNode 
{
public:
    KDTreeNode() = default;
    explicit KDTreeNode(const vector<Point *> &point_ptrs, vector<int> ranges, 
        bool is_leaf, unsigned level, unsigned dim, unsigned max_level) : 
        point_ptrs_(point_ptrs), ranges_(ranges), is_leaf_(is_leaf), 
        level_(level)
    {
        if (!is_leaf) BuildKDTreeNode_(dim, max_level);
    }
    const unsigned count() const { return point_ptrs_.size(); }
    const unsigned level() const { return level_; }
    const bool is_leaf() const { return is_leaf_; }
    const vector<Point *> point_ptrs() const { return point_ptrs_; }
    vector<shared_ptr<const KDTreeNode> > GetLeafNodePtrs() const 
    {
        vector<shared_ptr<const KDTreeNode> > node_ptrs; 
        GetLeafNodePtrs_(node_ptrs); 
        return node_ptrs;
    }
    double get_volume() const 
    {
        double result = 1; 
        for (unsigned i = 0; 2 * i < ranges_.size(); i++)
        {
            result *= (ranges_[2 * i] - ranges_[2 * i + 1]);
        }
        return result;
    }
private:
    void BuildKDTreeNode_(unsigned dim, unsigned max_level);
    void GetLeafNodePtrs_(
        vector<shared_ptr<const KDTreeNode> > &node_ptrs) const;
    // Pointers to points in this node
    vector<Point *> point_ptrs_;
    vector<int> ranges_;
    bool is_leaf_;
    unsigned level_;
    shared_ptr<KDTreeNode> lc_;
    shared_ptr<KDTreeNode> rc_;
};

class NewKDTree 
{
public:
    NewKDTree(const vector<Point> &points, unsigned max_level) :
        points_(points), node_ptrs_(static_cast<unsigned>(pow(2, max_level)))
    {
        const unsigned N = points.size();
        if (N == 0) 
            throw std::invalid_argument("N == 0");
        if (!IsPowerTwo(N))
            throw std::invalid_argument("N should be a power of 2");
        if (max_level == 0) 
            throw std::invalid_argument("max_level == 0");
        if (pow(2, max_level) > N)
            throw std::invalid_argument("2 ^ max_depth > N");
        dim_ = points_[0].dim();
        BuildKDTree_(max_level);
    }
    vector<Point> Sample(unsigned sample_size);
    vector<shared_ptr<const KDTreeNode> > get_node_ptrs() const 
    {
        return node_ptrs_;
    }
private:
    void BuildKDTree_(unsigned max_level);
    vector<Point> points_;
    // The leaf node pointers
    vector<shared_ptr<const KDTreeNode> > node_ptrs_;
    // All data points
    // Dimemsion of the kd-tree, namely k in the kd-tree
    unsigned dim_;
    KDTreeNode root_;
};

struct kdtree {
    struct kdtree *lc;
    struct kdtree *rc;
    int *range;
    unsigned nump; /* the number of points in the (sub)tree. */
    unsigned dim;
    unsigned level;
};



/*
 * Create a kd tree divided by data points
 * Arguments:
 * int **data: an array of pointers to int *
 * unsigned long n: the number of the pointers to int *
 * unsigned m: the length of template
 * unsigned dim: the current discrimiant
 * unsigned level: level of the kd tree node being created
 */
struct kdtree *build_kdtree(const int **data, unsigned long n,
                             unsigned m, unsigned dim,
                             unsigned level);
 
struct kdtree *build_kdtree_grid(const int *data, unsigned long N,
                                 unsigned m, unsigned p);

long long count_range_kdtree(struct kdtree *tree, const int *point, 
                             unsigned m, int r);

/* 
 * Create a kd tree node given range, m, ...
 *
 */
struct kdtree *_create_kdtree_node(
    const int *range, unsigned m, unsigned dim,
    unsigned level, unsigned long nump);
#endif // __KDTREE_H__
 