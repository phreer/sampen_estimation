#include <iostream>
#include <vector>

#include "utils.h"

using namespace std;
int main()
{
    vector<double> data(53453222, 1);
    double sum = ComputeSum(data);
    cout << sum << endl;
    return 0;
}