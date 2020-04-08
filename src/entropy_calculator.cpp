#include "entropy_calculator.h"
#include "kdtree.h"

double DynamicEntropyCalculator::ComputeEntropy(const vector<int> &data, 
                                                unsigned m, 
                                                unsigned max_level) 
{
    if (log2(data.size()) <= max_level - 1) 
        throw std::invalid_argument("log2(data.size()) <= max_level");
    
    vector<Point> points = get_points(data, m);
    NewKDTree tree(points, max_level);
    vector<shared_ptr<const KDTreeNode> > node_ptrs = tree.get_node_ptrs();

    double ent = 0; 
    for (unsigned i = 0; i < node_ptrs.size(); i++)
    {
        double p = node_ptrs[i]->count() / node_ptrs[i]->get_volume();
        std::cout << "volume: " << node_ptrs[i]->get_volume();
        std::cout << ", p: " << p << std::endl;
        ent -= p * log(p);
    }
    return ent;
}