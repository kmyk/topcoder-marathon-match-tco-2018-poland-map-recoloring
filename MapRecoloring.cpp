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

template <class RandomEngine>
vector<int> color_greedy(int R, vector<vector<int> > const & g, RandomEngine & gen) {
    vector<int> color(R, -1);
    vector<int> order(R);
    iota(ALL(order), 0);
    shuffle(ALL(order), gen);
    for (int i : order) {
        int used = 0;
        for (int j : g[i]) if (color[j] != -1) {
            used |= 1 << color[j];
        }
        color[i] = __builtin_ctz(~ used);
    }
    return color;
}

inline int get_C(vector<int> const & paint) {
    return *max_element(ALL(paint)) + 1;
}

int calculate_Pc(int R, vector<int> const & paint, vector<array<int, MAX_C> > const & old_color_count) {
    int Pc = 0;
    REP (i, R) {
        if (paint[i] < MAX_C) {
            Pc += old_color_count[i][paint[i]];
        }
    }
    return Pc;
}

pair<int, int> get_fewest_color(int C, vector<int> const & paint) {
    vector<int> cnt(C);
    for (int c : paint) {
        cnt[c] += 1;
    }
    int c = min_element(ALL(cnt)) - cnt.begin();
    return make_pair(c, cnt[c]);
}

vector<int> permute_paint(int R, int C0, int C, vector<int> const & paint, vector<array<int, MAX_C> > const & old_color_count) {
    // prepare delta[i][j]; it is the delta of score if color i becomes color j
    vector<array<int, MAX_C> > delta(C);
    REP (i, R) {
        REP (j, C0) {
            delta[paint[i]][j] += old_color_count[i][j];
        }
    }
    // try all permutations
    vector<int> sigma(C);  // inverted
    int highscore = INT_MIN;
    vector<int> tau;
    auto func = [&](int chosen) {
        // permute C0 colors
        REP (i, C) if (chosen & (1 << i)) {
            tau.push_back(i);
        }
        do {
            int score = 0;
            REP (i, C0) {
                score += delta[tau[i]][i];  // since tau is inverted here
            }
            if (highscore < score) {
                if ((int)tau.size() < C) {
                    REP (i, C) if (not (chosen & (1 << i))) {
                        tau.push_back(i);  // put remaining parts
                    }
                }
                highscore = score;
                sigma = tau;
            }
        } while (next_permutation(tau.begin(), tau.begin() + C0));
        tau.clear();
    };
    // choose C0 colors from C colors
    for (int x = (1 << C0) - 1; x < (1 << C); ) {  // enumerate x \subseteq C s.t. |x| = C0
        func(x);
        // update x
        int t = x | (x - 1);
        x = (t + 1) | (((~ t & - ~ t) - 1) >> (__builtin_ctz(x) + 1));
    }
    // apply tau = sigma^{-1}
    tau.resize(C, -1);
    REP (i, C) {
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
int get_random_paintable_color(int r, int c, int C, vector<int> const & paint, vector<vector<int> > const & g, RandomEngine & gen) {
    vector<bool> used(C);
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

void remove_unused_colors(int R, int & C, vector<int> & paint) {
    while (true) {
        int c, cnt; tie(c, cnt) = get_fewest_color(C, paint);
        if (cnt) break;
        REP (r, R) {
            assert (paint[r] != c);
            if (paint[r] > c) paint[r] -= 1;
        }
        C -= 1;
    }
}

vector<int> solve(int H, int W, int R, int C0, vector<int> const & regions, vector<int> const & old_colors) {
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
    cerr << "C0 = " << C0 << endl;

    xor_shift_128 gen(42);
    vector<vector<int> > g = construct_graph(H, W, R, regions);
    vector<array<int, MAX_C> > old_color_count = count_old_colors(H * W, R, regions, old_colors);

    vector<int> paint = color_greedy(R, g, gen);
    int C = get_C(paint);
    int Pc = calculate_Pc(R, paint, old_color_count);

    vector<int> answer;
    int answer_C = INT_MAX;
    int answer_Pc = INT_MIN;
    auto update_answer = [&]() {
        if (make_pair(- answer_C, answer_Pc) < make_pair(- C, Pc)) {
            answer = paint;
            answer_C = C;
            answer_Pc = Pc;
            cerr << "C = " << answer_C << ", Pc = " << answer_Pc << endl;
        }
    };
    update_answer();

    ll iteration = 0;
    vector<int> order(R);
    iota(ALL(order), 0);
    for (; iteration % 100 != 0 or rdtsc() - clock_begin < 0.95 * TLE; ++ iteration) {
        // shuffle non-target colors
        int c = get_fewest_color(C, paint).first;
        shuffle(ALL(order), gen);
        for (int r : order) {
            int nc = get_random_paintable_color(r, c, C, paint, g, gen);
            if (nc != -1) {
                paint[r] = nc;
            }
        }

        // update
        paint = permute_paint(R, C0, C, paint, old_color_count);
        Pc = calculate_Pc(R, paint, old_color_count);
        remove_unused_colors(R, C, paint);
        update_answer();
    }

    // debug print
    ll score = 100000ll * answer_C + H * W - answer_Pc;
    cerr << "color = " << C << endl;
    cerr << "P^c = " << answer_Pc << endl;
    cerr << "raw score = " << score << endl;
#ifdef LOCAL
    if (seed != -1) {
        cerr << "{\"seed\":" << seed
             << ",\"H\":" << H
             << ",\"W\":" << W
             << ",\"R\":" << R
             << ",\"C0\":" << C0
             << ",\"C\":" << answer_C  // the number of color, smaller is better
             << ",\"P\":" << H * W - answer_Pc  // the number of cells painted, smaller is better
             << ",\"Pc\":" << answer_Pc
             << ",\"rawScore\":" << score
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
        int C0 = *max_element(oldColors.begin(), oldColors.end()) + 1;
        return solve(H, W, R, C0, regions, oldColors);
    }
};
