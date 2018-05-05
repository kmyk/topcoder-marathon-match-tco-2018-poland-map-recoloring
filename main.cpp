#include "MapRecoloring.cpp"

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
    MapRecoloring *mr = new MapRecoloring();
    int H, S, R;
    cin >> H >> S;
    vector<int> regions(S);
    getVector(regions);
    cin >> R;
    vector<int> oldColors(R);
    getVector(oldColors);

    vector<int> ret = mr->recolor(H, regions, oldColors);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}
