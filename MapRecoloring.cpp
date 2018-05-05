#include <bits/stdc++.h>
#define REP(i, n) for (int i = 0; (i) < int(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < int(n); ++ (i))
#define REP_R(i, n) for (int i = int(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = int(n) - 1; (i) >= int(m); -- (i))
#define ALL(x) begin(x), end(x)

using ll = long long;
using namespace std;
template <class T> using reversed_priority_queue = priority_queue<T, vector<T>, greater<T> >;
template <class T> inline void chmax(T & a, T const & b) { a = max(a, b); }
template <class T> inline void chmin(T & a, T const & b) { a = min(a, b); }
template <typename T> ostream & operator << (ostream & out, vector<T> const & xs) { REP (i, int(xs.size()) - 1) out << xs[i] << ' '; if (not xs.empty()) out << xs.back(); return out; }

constexpr int MAX_H = 200;
constexpr int MAX_W = 200;
constexpr int MAX_R = 4000;
constexpr int MAX_C = 5;

vector<vector<int> > construct_graph(int H, int W, int R, vector<int> const & regions) {
    vector<vector<int> > g(R);
    auto func = [&](int z, int nz) {
        if (regions[z] != regions[nz]) {
            g[regions[z]].push_back(regions[nz]);
            g[regions[nz]].push_back(regions[z]);
        }
    };
    REP (y, H) REP (x, W) {
        int z = y * W + x;
        if (x + 1 < W) func(z, z + 1);
        if (y + 1 < H) func(z, z + W);
    }
    REP (i, R) {
        sort(ALL(g[i]));
        g[i].erase(unique(ALL(g[i])), g[i].end());
    }
    return g;
}

vector<array<int, MAX_C> > count_old_colors(int HW, int R, vector<int> const & regions, vector<int> const & old_colors) {
    vector<array<int, MAX_C> > cnt(R);
    REP (z, HW) {
        cnt[regions[z]][old_colors[z]] += 1;
    }
    return cnt;
}

int mex_destruct(vector<int> & xs) {
    int y = 0;
    sort(xs.rbegin(), xs.rend());
    while (not xs.empty() and y >= xs.back()) {
        if (y == xs.back()) ++ y;
        xs.pop_back();
    }
    return y;
}

vector<int> color_greedy(int R, vector<vector<int> > const & g) {
    vector<int> color(R, -1);
    vector<int> used;
    REP (i, R) {
        for (int j : g[i]) if (color[j] != -1) {
            used.push_back(color[j]);
        }
        color[i] = mex_destruct(used);
        used.clear();
    }
    return color;
}

int calculate_score_delta(int R, vector<int> const & paint, vector<array<int, MAX_C> > const & old_color_count) {
    int delta = 0;
    REP (i, R) {
        if (paint[i] < MAX_C) {
            delta += old_color_count[i][paint[i]];
        }
    }
    return delta;
}

vector<int> permute_paint(int R, int C, int k, vector<int> const & paint, vector<array<int, MAX_C> > const & old_color_count) {
    // prepare delta[i][j]; it is the delta of score if color i becomes color j
    vector<array<int, MAX_C> > delta(k);
    REP (i, R) {
        REP (j, C) {
            delta[paint[i]][j] += old_color_count[i][j];
        }
    }
    // try all permutations
    vector<int> sigma(k);  // inverted
    int highscore = INT_MIN;
    vector<int> tau;
    auto func = [&](int chosen) {
        // permute C colors
        REP (i, k) if (chosen & (1 << i)) {
            tau.push_back(i);
        }
        do {
            int score = 0;
            REP (i, C) {
                score += delta[tau[i]][i];  // since tau is inverted here
            }
            if (highscore < score) {
                if ((int)tau.size() < k) {
                    REP (i, k) if (not (chosen & (1 << i))) {
                        tau.push_back(i);  // put remaining parts
                    }
                }
                highscore = score;
                sigma = tau;
            }
        } while (next_permutation(tau.begin(), tau.begin() + C));
        tau.clear();
    };
    // choose C colors from k colors
    for (int x = (1 << C) - 1; x < (1 << k); ) {  // enumerate x \subseteq k s.t. |x| = C
        func(x);
        // update x
        int t = x | (x - 1);
        x = (t + 1) | (((~ t & - ~ t) - 1) >> (__builtin_ctz(x) + 1));
    }
    // apply tau = sigma^{-1}
    tau.resize(k, -1);
    REP (i, k) {
        tau[sigma[i]] = i;
    }
    vector<int> npaint(R);
    REP (i, R) {
        npaint[i] = tau[paint[i]];
    }
    return npaint;
}

vector<int> apply_permutation(vector<int> const & sigma, vector<int> xs) {
    for (int & x : xs) {
        x = sigma[x];
    }
    return xs;
}

vector<int> solve(int H, int W, int R, int C, vector<int> const & regions, vector<int> const & old_colors) {
    cerr << "H = " << H << endl;
    cerr << "W = " << W << endl;
    cerr << "R = " << R << endl;
    cerr << "C = " << C << endl;

    vector<vector<int> > g = construct_graph(H, W, R, regions);
    vector<int> paint = color_greedy(R, g);
    int k = *max_element(ALL(paint)) + 1;
    cerr << "the number of color = " << k << endl;

    vector<array<int, MAX_C> > old_color_count = count_old_colors(H * W, R, regions, old_colors);
    cerr << "the sum of delta = " << calculate_score_delta(R, paint, old_color_count) << "  (before permutation)" << endl;
    paint = permute_paint(R, C, k, paint, old_color_count);
    int delta = calculate_score_delta(R, paint, old_color_count);
    ll score = 100000ll * k - H * W + delta;
    cerr << "the sum of delta = " << delta << endl;
    cerr << "the raw score = " << score << endl;
    return paint;
}


class MapRecoloring {
public:
    vector<int> recolor(int H, vector<int> regions, vector<int> oldColors) {
        int W = regions.size() / H;
        int R = *max_element(regions.begin(), regions.end()) + 1;
        int C = *max_element(oldColors.begin(), oldColors.end()) + 1;
        return solve(H, W, R, C, regions, oldColors);
    }
};
