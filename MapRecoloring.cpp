// C++11
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class MapRecoloring {
public:
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
// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
    MapRecoloring mr;
    int H, S, R;
    cin >> H >> S;
    vector<int> regions(S);
    getVector(regions);
    cin >> R;
    vector<int> oldColors(R);
    getVector(oldColors);

    vector<int> ret = mr.recolor(H, regions, oldColors);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}
