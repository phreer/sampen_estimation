/* file: kdtree.h
 * date: 2019-12-10
 * author: phree
 *
 * description: data structure and functions for kd tree
 */
#ifndef __KDTREE_H__
#define __KDTREE_H__

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
 