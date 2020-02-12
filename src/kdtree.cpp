#include <iostream>

#include <string.h>
#include "utils.h"
#include "kdtree.h"
#include "random_sampler.h"

void KDTreeNode::BuildKDTreeNode_(unsigned dim, unsigned max_level)
{
    vector<Point *> point_ptrs(point_ptrs_);
    unsigned k = level_ % dim;
    unsigned level = level_ + 1;
    bool is_leaf = (level == max_level);
    std::sort(point_ptrs.begin(), point_ptrs.end(), 
        [k] (const Point *p1, const Point *p2) 
        {
            return (*p1)[k] < (*p2)[k];
        });
    unsigned mid = point_ptrs.size() / 2;
    
    vector<int> ranges(ranges_);
    ranges[2 * k] = (*point_ptrs[0])[k];
    ranges[2 * k + 1] = (*point_ptrs[mid - 1])[k];
    lc_ = std::make_shared<KDTreeNode>(
        vector<Point *>(point_ptrs.cbegin(), point_ptrs.cbegin() + mid), 
        ranges, is_leaf, level, dim, max_level);
    
    ranges[2 * k] = (*point_ptrs[mid])[k];
    ranges[2 * k + 1] = (*point_ptrs[point_ptrs.size() - 1])[k];
    rc_ = std::make_shared<KDTreeNode>(
        vector<Point *>(point_ptrs.cbegin() + mid, point_ptrs.cend()), 
        ranges, is_leaf, level, dim, max_level);
}

void KDTreeNode::GetLeafNodePtrs_(
    vector<shared_ptr<KDTreeNode> > &node_ptrs) const
{
    if (lc_->is_leaf()) node_ptrs.push_back(lc_);
    else lc_->GetLeafNodePtrs_(node_ptrs);
    if (rc_->is_leaf()) node_ptrs.push_back(rc_);
    else rc_->GetLeafNodePtrs_(node_ptrs);
}

void NewKDTree::BuildKDTree_(unsigned max_level) 
{
    vector<int> ranges(2 * dim_);
    for (int i = 0; i < dim_; i++) 
    {
        ranges[2 * i] = std::min_element(points_.cbegin(), points_.cend(), 
            [i] (const Point &p1, const Point &p2) { 
                return p1[i] < p2[i]; 
            })->operator[](i);
        ranges[2 * i + 1] = std::max_element(points_.cbegin(), points_.cend(), 
            [i] (const Point &p1, const Point &p2) { 
                return p1[i] < p2[i]; 
            })->operator[](i);
    }
    vector<Point *> point_ptrs(points_.size());
    for (int i = 0; i < point_ptrs.size(); i++) 
    point_ptrs[i] = &points_[i];
    bool is_leaf = (max_level == 0);
    root_ = KDTreeNode(point_ptrs, ranges, is_leaf, 0, dim_, max_level);
    node_ptrs_ = root_.GetLeafNodePtrs();
}

vector<Point> NewKDTree::Sample(unsigned sample_size) 
{
    if (sample_size != node_ptrs_.size()) 
    throw std::logic_error("sample_size != node_ptrs_.size()");
    vector<Point> result(sample_size);
    uniform_int_generator uig(0, node_ptrs_[0]->count() - 1, 
        uniform_int_generator::PSEUDO, true);
    for (int i = 0; i < sample_size; i++) 
    {
        result[i] = *(node_ptrs_[i]->point_ptrs()[uig.get()]);
    }
    return result;
}

/* Some utility functions */
int _compar_int(const void *p1, const void *p2)
{
    return *(int *)p1 - *(int *)p2;
}

int _compar_intp(const void *p1, const void *p2)
{
    return **(int **)p1 - **(int **)p2;
}

void shift(int **data, unsigned long np, int dim)
{
    unsigned long i;
    for (i = 0; i < np; i++)
        data[i] += dim;
}

void* _get_max_min(void **numsp, unsigned long N, 
                   int (*compar)(const void *p1, const void *p2)) 
{
    if (N == 0) return NULL;
    void **result = (void **)malloc(2 * sizeof(void *));

    result[0] = result[1] = *numsp;
    unsigned long i = 0;
    for (i = 0; i < N / 2; i++)
    {   
        unsigned long idx = 2 * i;
        unsigned long maxi = compar(numsp + idx, numsp + idx + 1) > 0 ? 
            idx : idx + 1;
        unsigned long mini = 4 * i + 1 - maxi;
        if (compar(numsp + maxi, result) > 0) 
            result[0] = *(numsp + maxi);
        if (compar(result + 1, numsp + mini) > 0) 
            result[1] = *(numsp + mini);
    }
    if (N % 2) 
    {
        unsigned long idx = N - 1;
        if (compar(numsp + idx, result) > 0) result[0] = *(numsp + idx);
        if (compar(result + 1, numsp + idx) > 0) result[1] = *(numsp + idx);
    }
    return result;
}


/*
 * create a kd tree divided by data points
 * Arguments:
 *     int **data: an array of pointers to int *
 *     unsigned long n: the number of the pointers to int *
 *     unsigned m: the length of template
 *     unsigned dim: the current discrimiant
 *     unsigned level: level of the kd tree node being created
 * Returns:
 *     a point ter to the built kd tree
 */
struct kdtree *build_kdtree(const int **data, unsigned long n,
                           unsigned m, unsigned dim,
                           unsigned level)
{
    if (n <= 0)
        return NULL;
    unsigned i;
    int stop = 1;

    int **_data = static_cast<int **>(malloc(n * sizeof(int *)));
    memcpy(_data, data, n * sizeof(int *));

    /* sort by coordinate dim */
    shift(_data, n, dim);
    qsort(_data, n, sizeof(int *), _compar_intp);
    shift(_data, n, -dim);

    int *range = (int *)malloc(2 * m * sizeof(int));
    range[2 * dim] = _data[0][dim];
    range[2 * dim + 1] = _data[n - 1][dim];

    for (i = 0; i < m; i++)
    {
        if (i != dim)
        {
            auto max_min = static_cast<const int ** const>(
                _get_max_min((void **)_data, n, _compar_intp));
            range[2 * i] = *max_min[1];
            range[2 * i + 1] = *max_min[0];
            free(max_min);
        }
        if (range[2 * i] != range[2 * i + 1])
            stop = 0;
        shift(_data, n, 1);
    }
    shift(_data, n, -m);

    struct kdtree *result = _create_kdtree_node(range, m, dim, level, n);

    if (!stop)
    {
        result->lc = build_kdtree(
            (const int **)_data, n / 2, m, (dim + 1) % m, level + 1);
        result->rc = build_kdtree(
            (const int **)_data + n / 2, (n + 1) / 2, m, (dim + 1) % m, level + 1);
    }
    
    free(range);
    free(_data);
    return result;
}


struct kdtree *_create_kdtree_node(
    const int *range, unsigned m, unsigned dim,
    unsigned level, unsigned long nump)
{
    struct kdtree *result = (struct kdtree *)malloc(sizeof(struct kdtree));
    result->range = (int *)malloc(2 * m * sizeof(int));
    memcpy(result->range, range, 2 * m * sizeof(int));
    result->lc = result->rc = NULL;
    result->nump = nump;
    result->dim = dim;
    result->level = level;
    return result;
}


/*
 * Insert a data point of length m to 
 */
void insert_point(struct kdtree *tree, const int *point,
                  unsigned m, unsigned p)
{
    unsigned i = 0, j = 0;
    struct kdtree *curr = tree;
    curr->nump++;

    for (i = 0; i < p; i++)
    {
        for (j = 0; j < m; j++)
        {
            int med = curr->range[2 * j] + curr->range[2 * j + 1];
            med >>= 1;
            if (point[j] <= med)
            {
                if (!curr->lc)
                {
                    curr->lc = _create_kdtree_node(
                        curr->range, m, (j + 1) % m, i * m + j + 1, 0);
                    curr->lc->range[2 * j + 1] = med;
                }
                curr = curr->lc;
            }
            else
            {
                if (!curr->rc)
                {
                    curr->rc = _create_kdtree_node(
                        curr->range, m, (j + 1) % m, i * m + j + 1, 0);
                    curr->rc->range[2 * j] = med + 1;
                }
                curr = curr->rc;
            }
            curr->nump++;
        }
    }
}

struct kdtree *build_kdtree_grid(const int *data, unsigned long N,
                                 unsigned m, unsigned p)
{
    unsigned long i = 0;
    if (N < m)
        return NULL;
    struct kdtree *result = (struct kdtree *)malloc(sizeof(struct kdtree));
    result->nump = 0;
    result->dim = 0;
    result->lc = result->rc = 0;
    result->level = 0;
    result->range = (int *)malloc(2 * m * sizeof(int));
    for (i = 0; i < m; i++)
    {
        result->range[2 * i] = 0;
        result->range[2 * i + 1] = (1 << p) - 1;
    }

    for (i = 0; i < N - m + 1; i++)
    {
        insert_point(result, data + i, m, p);
    }
    return result;
}

long long count_range_kdtree(struct kdtree *tree, const int *point,
                             unsigned m, int r)
{
    /* case 0, [point - r, point + r] does NOT intersect the range of tree
     * case 1, the range of tree is within [point - r, point + r] 
     * case 2, the range of tree intersects [point - r, point + r] and 
     * is NOT contained in it
     */
    if (!tree) return 0;
    enum CASE
    {
        NOT_INTER,
        WITHIN,
        INTER
    };
    enum CASE _case = WITHIN;
    unsigned i;
    for (i = 0; i < m; i++)
    {
        if (tree->range[2 * i] > point[i] + r ||
            tree->range[2 * i + 1] < point[i] - r)
        {
            _case = NOT_INTER;
            break;
        }
        if (tree->range[2 * i] < point[i] - r ||
            tree->range[2 * i + 1] > point[i] + r)
        {
            _case = INTER;
        }
    }
    switch (_case)
    {
    case NOT_INTER:
        return 0;
    case WITHIN:
        return tree->nump;
    case INTER:
    {
        long long count = 0;
        if (tree->lc) count += count_range_kdtree(tree->lc, point, m, r);
        if (tree->rc) count += count_range_kdtree(tree->rc, point, m, r);
        return count;
    }
    }
}