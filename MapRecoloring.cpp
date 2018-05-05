#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class MapRecoloring {
public:
    MapRecoloring() = default;
    vector<int> recolor(int H, vector<int> regions, vector<int> oldColors) {
        // number of regions = max element in regions + 1
        int reg = *max_element(regions.begin(), regions.end()) + 1;
        vector<int> ret(reg);
        for (int i = 0; i < reg; ++i) {
            ret[i] = i;
        }
        return ret;
    }
};
