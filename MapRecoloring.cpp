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

/**
 * @note C < 32 is required
 */
uint32_t get_unpaintablity_array(int r, vector<int> const & paint, vector<vector<int> > const & g) {
    uint32_t used = 0;
    for (int r1 : g[r]) {
        used |= 1u << paint[r1];
    }
    return used;
}

/**
 * @note returned values are in [0, 32), or -1 if x = 0
 */
template <class RandomEngine>
int get_random_bit(uint32_t x, RandomEngine & gen) {
    int cnt = __builtin_popcount(x);
    if (cnt == 0) return -1;
    int i = uniform_int_distribution<int>(0, cnt - 1)(gen);
    while (true) {
        int lsb = x & - x;
        if (i -- == 0) return __builtin_ctz(lsb);
        x &= ~ lsb;
    }
}

template <class RandomEngine>
vector<int> color_greedy(int R, vector<vector<int> > const & g, RandomEngine & gen) {
    vector<int> paint(R, -1);
    vector<int> order(R);
    iota(ALL(order), 0);
    shuffle(ALL(order), gen);
    for (int r : order) {
        uint32_t used = 0;
        for (int r1 : g[r]) if (paint[r1] != -1) {
            used |= 1u << paint[r1];
        }
        paint[r] = __builtin_ctz(~ used);
    }
    return paint;
}

inline int get_C(vector<int> const & paint) {
    return *max_element(ALL(paint)) + 1;
}

int calculate_P(int HW, int R, vector<int> const & paint, vector<array<int, MAX_C> > const & old_color_count) {
    int P = HW;
    REP (i, R) {
        if (paint[i] < MAX_C) {
            P -= old_color_count[i][paint[i]];
        }
    }
    return P;
}

vector<int> get_color_frequency(int C, vector<int> const & paint) {
    vector<int> cnt(C);
    for (int c : paint) {
        cnt[c] += 1;
    }
    return cnt;
}

pair<int, int> get_fewest_color(int C, vector<int> const & paint) {
    vector<int> cnt = get_color_frequency(C, paint);
    int c = min_element(ALL(cnt)) - cnt.begin();
    return make_pair(c, cnt[c]);
}

vector<int> get_primary_color(int R, vector<array<int, MAX_C> > const & old_color_count) {
    vector<int> primary_color(R);
    REP (r, R) {
        primary_color[r] = min_element(ALL(old_color_count[r])) - old_color_count[r].begin();
    }
    return primary_color;
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

int remove_unused_colors(int R, int C, vector<int> & paint) {
    while (true) {
        int c, cnt; tie(c, cnt) = get_fewest_color(C, paint);
        if (cnt) break;
        REP (r, R) {
            assert (paint[r] != c);
            if (paint[r] > c) paint[r] -= 1;
        }
        C -= 1;
    }
    return C;
}

vector<int> list_target_regions(int R, int C, vector<int> const & paint) {
    vector<int> lookup;
    REP (r, R) {
        if (paint[r] == C - 1) {
            lookup.push_back(r);
        }
    }
    return lookup;
}

constexpr int R_LIMIT = 0;
constexpr int C_BASE = 7;

double get_score(int R, int C0, int C, int P, vector<int> const & freq) {
    if (C > C_BASE) {
        return freq.back() * 0.01;
    } else if (C == C_BASE and R < R_LIMIT) {
        return freq.back() * 0.1;
    } else {
        return P * 0.05;
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
// xor_shift_128 gen((random_device()()));
    vector<vector<int> > g = construct_graph(H, W, R, regions);
    vector<array<int, MAX_C> > old_color_count = count_old_colors(H * W, R, regions, old_colors);
    vector<int> primary_color = get_primary_color(R, old_color_count);

    vector<int> answer;
    int answer_C = INT_MAX;
    int answer_P = INT_MIN;
    auto update_answer = [&](vector<int> const & paint, int C, int P) {
        if (make_pair(C, P) < make_pair(answer_C, answer_P)) {
            answer = paint;
            answer_C = C;
            answer_P = P;
            vector<int> freq = get_color_frequency(C, paint);
#ifdef LOCAL
            cerr << "C = " << answer_C << ", P = " << answer_P << ", freq = (" << freq << ")" << endl;
#endif
        }
    };

    vector<int> paint = color_greedy(R, g, gen);
    int C = get_C(paint);
    paint = permute_paint(R, C0, C, paint, old_color_count);
    int P = calculate_P(H * W, R, paint, old_color_count);
    update_answer(paint, C, P);
    vector<int> freq = get_color_frequency(C, paint);
    vector<int> lookup = list_target_regions(R, C, paint);
    double score = get_score(R, C0, C, P, freq);

    ll iteration = 0;
    double t = rdtsc() - clock_begin;
    vector<int> order(R);
    iota(ALL(order), 0);
    for (; t < 0.95 * TLE; ++ iteration) {
        if (iteration % 101 == 0) t = rdtsc() - clock_begin;
        double temperature = 1 - t / TLE;

// if(iteration % 100000 == 0)
// cerr << "C = " << answer_C << ", P = " << answer_P << ", freq = (" << freq << ")" << endl;
        // modify a cell
        int r, prev_paint_r;
        int lookup_index = -1;
        if (C > C_BASE or (C == C_BASE and R < R_LIMIT)) {
            lookup_index = uniform_int_distribution<int>(0, lookup.size() - 1)(gen);
            r = lookup[lookup_index];
        } else {
            r = uniform_int_distribution<int>(0, R - 1)(gen);
        }
        while (true) {
            prev_paint_r = paint[r];
            uint32_t used = get_unpaintablity_array(r, paint, g);
            if ((C < C_BASE or (C == C_BASE and R >= R_LIMIT))
                    and paint[r] != primary_color[r]
                    and not (used & (1u << primary_color[r]))
                    and bernoulli_distribution(0.5)(gen)) {
                paint[r] = primary_color[r];
                break;
            } else {
                if (C > C_BASE) used |= (1u << (C - 1));
                if (C == C_BASE and R < R_LIMIT and bernoulli_distribution(0.9)(gen)) used |= (1u << (C - 1));
                used |= (1u << paint[r]);
                used ^= (1u << C) - 1;
                if (used) {
                    paint[r] = get_random_bit(used, gen);
                    break;
                }
            }
            r = g[r][uniform_int_distribution<int>(0, g[r].size() - 1)(gen)];
        }
        assert (prev_paint_r != paint[r]);

        // update
        int next_C = remove_unused_colors(R, C, paint);  // destructive
        assert (next_C <= C);
        if (next_C < C) {
            C = next_C;
            paint = permute_paint(R, C0, C, paint, old_color_count);
            P = calculate_P(H * W, R, paint, old_color_count);
            update_answer(paint, C, P);
            score = get_score(R, C0, C, P, freq);
            freq = get_color_frequency(C, paint);
            lookup = list_target_regions(R, C, paint);
        } else {
            int next_P = P;
            if (prev_paint_r < MAX_C) next_P += old_color_count[r][prev_paint_r];
            if (paint[r]     < MAX_C) next_P -= old_color_count[r][paint[r]];
            update_answer(paint, next_C, next_P);
            freq[prev_paint_r] -= 1;
            freq[paint[r]] += 1;
            double next_score = get_score(R, C0, C, next_P, freq);
            double delta = score - next_score;
            if (delta >= 0 or bernoulli_distribution(exp(delta / temperature))(gen)) {
                P = next_P;
                score = next_score;
                if (C >= 7 or R >= R_LIMIT) {
                    if (prev_paint_r == C - 1) {
                        if (lookup_index == -1 or lookup[lookup_index] != r) {
                            lookup_index = find(ALL(lookup), r) - lookup.begin();
                        }
                        swap(lookup[lookup_index], lookup.back());
                        lookup.pop_back();
                    }
                    if (paint[r] == C - 1) {
                        lookup.push_back(r);
                    }
                }
            } else {
                freq[prev_paint_r] += 1;
                freq[paint[r]] -= 1;
                paint[r] = prev_paint_r;
            }
        }
    }

    // debug print
    ll raw_score = 100000ll * answer_C + answer_P;
    cerr << "color = " << C << endl;
    cerr << "P = " << answer_P << endl;
    cerr << "raw score = " << raw_score << endl;
#ifdef LOCAL
    if (seed != -1) {
        cerr << "{\"seed\":" << seed
             << ",\"H\":" << H
             << ",\"W\":" << W
             << ",\"R\":" << R
             << ",\"C0\":" << C0
             << ",\"C\":" << answer_C  // the number of color, smaller is better
             << ",\"P\":" << answer_P  // the number of cells painted, smaller is better
             << ",\"rawScore\":" << raw_score
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
