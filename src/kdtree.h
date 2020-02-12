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
    vector<shared_ptr<KDTreeNode> > GetLeafNodePtrs() const 
    {
        vector<shared_ptr<KDTreeNode> > node_ptrs; 
        GetLeafNodePtrs_(node_ptrs); 
        return node_ptrs;
    }
private:
    void BuildKDTreeNode_(unsigned dim, unsigned max_level);
    void GetLeafNodePtrs_(vector<shared_ptr<KDTreeNode> > &node_ptrs) const;
    vector<Point *> point_ptrs_;
    bool is_leaf_;
    shared_ptr<KDTreeNode> lc_;
    shared_ptr<KDTreeNode> rc_;
    vector<int> ranges_;
    unsigned level_;
};

class NewKDTree 
{
public:
    NewKDTree(const vector<Point> &points, unsigned dim, int max_level) :
        points_(points), dim_(dim), node_ptrs_(static_cast<unsigned>(pow(
            2, max_level)))
    {
        const unsigned N = points.size();
        if (pow(2, max_level) > N)
            throw std::logic_error("2 ^ max_depth > N");
        BuildKDTree_(max_level);
    }
    vector<Point> Sample(unsigned sample_size);
private:
    void BuildKDTree_(unsigned max_level);
    KDTreeNode root_;
    // The number of leaf node pointers.
    vector<shared_ptr<KDTreeNode> > node_ptrs_;
    vector<Point> points_;
    unsigned dim_;
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
 