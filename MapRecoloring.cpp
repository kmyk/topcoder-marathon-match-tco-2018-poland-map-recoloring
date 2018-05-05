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

class xor_shift_128 {
public:
    typedef uint32_t result_type;
    xor_shift_128(uint32_t seed) {
        set_seed(seed);
    }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator() () {
        uint32_t t = (a ^ (a << 11));
        a = b; b = c; c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }
    static constexpr uint32_t max() { return numeric_limits<result_type>::max(); }
    static constexpr uint32_t min() { return numeric_limits<result_type>::min(); }
private:
    uint32_t a, b, c, d;
};

constexpr double ticks_per_sec = 2800000000;
constexpr double ticks_per_sec_inv = 1.0 / ticks_per_sec;
inline double rdtsc() { // in seconds
    uint32_t lo, hi;
    asm volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return (((uint64_t)hi << 32) | lo) * ticks_per_sec_inv;
}
constexpr int TLE = 10; // sec

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

template <class RandomEngine>
vector<int> color_greedy(int R, vector<vector<int> > const & g, RandomEngine & gen) {
    vector<int> color(R, -1);
    vector<int> used;
    vector<int> order(R);
    iota(ALL(order), 0);
    shuffle(ALL(order), gen);
    for (int i : order) {
        for (int j : g[i]) if (color[j] != -1) {
            used.push_back(color[j]);
        }
        color[i] = mex_destruct(used);
        used.clear();
    }
    return color;
}

inline int get_k(vector<int> const & paint) {
    return *max_element(ALL(paint)) + 1;
}

int calculate_score_skipped(int R, vector<int> const & paint, vector<array<int, MAX_C> > const & old_color_count) {
    int skipped = 0;
    REP (i, R) {
        if (paint[i] < MAX_C) {
            skipped += old_color_count[i][paint[i]];
        }
    }
    return skipped;
}

void prepare_color_regin_lookup(int R, int k, vector<int> const & paint, vector<vector<int> > & lookup) {
    lookup.resize(k);
    REP (c, k) lookup[c].clear();
    REP (i, R) lookup[paint[i]].push_back(i);
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

/**
 * @param r is an index of a region
 * @param c is a forbidden color; ignored if -1
 * @note return -1 if there are no paintable colors
 */
template <class RandomEngine>
int get_random_paintable_color(int r, int c, int k, vector<int> const & paint, vector<vector<int> > const & g, RandomEngine & gen) {
    vector<bool> used(k);
    if (c != -1) used[c] = true;
    for (int j : g[r]) {
        used[paint[j]] = true;
    }
    int cnt = count(ALL(used), 0);
    if (cnt == 0) return -1;
    int nth = uniform_int_distribution<int>(0, cnt - 1)(gen);
    int nc = 0;
    while (true) {
        while (used[nc]) ++ nc;
        if (not nth) break;
        -- nth;
        ++ nc;
    }
    return nc;
}

void remove_an_unused_color(int unused_c, int R, int & k, vector<int> & paint, vector<vector<int> > & lookup) {
    REP (r, R) {
        assert (paint[r] != unused_c);
        if (paint[r] > unused_c) paint[r] -= 1;
    }
    REP3 (c, unused_c, k - 1) {
        lookup[c].swap(lookup[c + 1]);
    }
    lookup.pop_back();
    k -= 1;
}

void remove_unused_colors(int R, int & k, vector<int> & paint, vector<vector<int> > & lookup) {
    for (int c = 0; c < k; ++ c) {  // k may be modified in the loop
        if (lookup[c].empty()) {
            remove_an_unused_color(c, R, k, paint, lookup);
        }
    }
}

vector<int> solve(int H, int W, int R, int C, vector<int> const & regions, vector<int> const & old_colors) {
    double clock_begin = rdtsc();

    // debug print
#ifdef LOCAL
    ll seed = -1;
    if (getenv("SEED")) {
        seed = atoll(getenv("SEED"));
    }
    cerr << "seed = " << seed << endl;
#endif
    cerr << "H = " << H << endl;
    cerr << "W = " << W << endl;
    cerr << "R = " << R << endl;
    cerr << "C = " << C << endl;

    xor_shift_128 gen(42);
    vector<vector<int> > g = construct_graph(H, W, R, regions);
    vector<array<int, MAX_C> > old_color_count = count_old_colors(H * W, R, regions, old_colors);

    vector<int> paint = color_greedy(R, g, gen);
    int k = get_k(paint);
    int skipped = calculate_score_skipped(R, paint, old_color_count);
    vector<vector<int> > lookup; prepare_color_regin_lookup(R, k, paint, lookup);

    vector<int> answer;
    int answer_k = INT_MAX;
    int answer_skipped = INT_MIN;
    auto update_answer = [&]() {
        if (make_pair(- answer_k, answer_skipped) < make_pair(- k, skipped)) {
            answer = paint;
            answer_k = k;
            answer_skipped = skipped;
            cerr << "k = " << answer_k << ", skipped = " << answer_skipped << endl;
        }
    };
    update_answer();

    ll iteration = 0;
    for (; iteration % 100 != 0 or rdtsc() - clock_begin < 0.95 * TLE; ++ iteration) {
        int c = min_element(ALL(lookup), [&](vector<int> const & a, vector<int> const & b) { return a.size() < b.size(); }) - lookup.begin();

        if (iteration % 2 == 0) {
            // shuffle non-target colors
            vector<int> order(R);
            iota(ALL(order), 0);
            shuffle(ALL(order), gen);
            for (int r : order) {
                int nc = get_random_paintable_color(r, c, k, paint, g, gen);
                if (nc != -1) {
                    paint[r] = nc;
                }
            }
            // update
            paint = permute_paint(R, C, k, paint, old_color_count);
            skipped = calculate_score_skipped(R, paint, old_color_count);
            prepare_color_regin_lookup(R, k, paint, lookup);

        } else {
            // remove target color
            for (int i = 0; i < (int)lookup[c].size(); ++ i) {  // lookup[c].size() may be modified in the loop
                int r = lookup[c][i];
                int nc = get_random_paintable_color(r, c, k, paint, g, gen);
                if (nc != -1) {
                    paint[r] = nc;
                    swap(lookup[c][i], lookup[c].back());
                    lookup[c].pop_back();
                    lookup[nc].push_back(r);
                    -- i;
                }
            }
        }
        remove_unused_colors(R, k, paint, lookup);
        update_answer();
    }

    // debug print
    ll score = 100000ll * answer_k + H * W - answer_skipped;
    cerr << "the number of color = " << k << endl;
    cerr << "the sum of skipped = " << answer_skipped << endl;
    cerr << "the raw score = " << score << endl;
#ifdef LOCAL
    if (seed != -1) {
        cerr << "{\"seed\":" << seed
             << ",\"H\":" << H
             << ",\"W\":" << W
             << ",\"R\":" << R
             << ",\"C\":" << C
             << ",\"k\":" << answer_k  // the number of color, smaller is better
             << ",\"skipped\":" << answer_skipped  // the number of cells which skipped to paint, larger is better
             << ",\"score\":" << score  // smaller is better
             << ",\"iteration\":" << iteration
             << ",\"time\":" << rdtsc() - clock_begin
             << "}" << endl;
    }
#endif

    return answer;
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
