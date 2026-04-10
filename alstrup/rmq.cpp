#ifndef ALSTRUP_RMQ
#define ALSTRUP_RMQ

#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
using namespace std;
#define fwd(i, a, n) for (int i = (a); i < (n); i++)
#define rep(i, n) fwd(i, 0, n)
#define sz(X) int(ssize(X))
#define eb emplace_back
using pii = pair<int, int>; using vi = vector<int>;
using ll = long long; using ld = long double;

template <uint32_t block_bit_len, class T>
struct SparseTable {
    vector<T> data;
    vector<uint8_t> indices;
    vector<vector<uint8_t>> st;
    uint32_t num_blocks;
    uint32_t max_level;
    uint32_t min_idx = 0;

    SparseTable() = default;

    inline uint32_t log2_floor(uint32_t i) const {
        return i ? 31 - __builtin_clz(i) : 0;
    }

    inline uint32_t read_bits(const vector<uint8_t>& mem, size_t bit_offset, uint32_t k) const {
        if (k == 0) return 0;
        size_t byte_idx = bit_offset / 8;
        uint32_t bit_shift = bit_offset % 8;
        uint64_t block;
        std::memcpy(&block, &mem[byte_idx], sizeof(uint64_t));
        return (block >> bit_shift) & ((1ULL << k) - 1);
    }

    inline void write_bits(vector<uint8_t>& mem, size_t bit_offset, uint32_t k, uint32_t val) {
        if (k == 0) return;
        size_t byte_idx = bit_offset / 8;
        uint32_t bit_shift = bit_offset % 8;
        uint64_t block;
        std::memcpy(&block, &mem[byte_idx], sizeof(uint64_t));
        uint64_t mask = ((1ULL << k) - 1);
        block &= ~(mask << bit_shift);
        block |= (static_cast<uint64_t>(val) & mask) << bit_shift;
        std::memcpy(&mem[byte_idx], &block, sizeof(uint64_t));
    }

    size_t getSize() const {
        size_t total_bytes = sizeof(*this);
        
        total_bytes += data.capacity() * sizeof(T);
        total_bytes += indices.capacity() * sizeof(uint8_t);
        
        total_bytes += st.capacity() * sizeof(vector<uint8_t>);
        
        for (const auto& row : st) {
            total_bytes += row.capacity() * sizeof(uint8_t);
        }
        
        return total_bytes;
    }

    SparseTable(vector<T> data_, vector<uint8_t> indices_)
        : data(std::move(data_)), indices(std::move(indices_)) {
        num_blocks = sz(data);
        if (num_blocks == 0) return;

        max_level = log2_floor(num_blocks);
        st.resize(max_level + 1);

        for (uint32_t k = 1; k <= max_level; ++k) {
            size_t num_entries = num_blocks - (1 << k) + 1;
            size_t bits_needed = num_entries * k;
            size_t bytes_needed = (bits_needed + 7) / 8;
            st[k].resize(bytes_needed + 8, 0);

            uint32_t half_range = 1 << (k - 1);
            
            for (uint32_t i = 0; i < num_entries; ++i) {
                uint32_t offset_left = (k == 1) ? 0 : read_bits(st[k - 1], i * (k - 1), k - 1);
                uint32_t offset_right = (k == 1) ? 0 : read_bits(st[k - 1], (i + half_range) * (k - 1), k - 1);
                
                uint32_t abs_left = i + offset_left;
                uint32_t abs_right = i + half_range + offset_right;

                uint32_t winning_offset = (data[abs_left] <= data[abs_right]) ? offset_left : (half_range + offset_right); 
                write_bits(st[k], i * k, k, winning_offset);
            }
        }
        for (uint32_t idx = 0; idx < data.size(); ++idx) {
            if (data[idx] < data[min_idx]) min_idx = idx;
        }
    }

    uint64_t get(uint32_t s, uint32_t e) const {
        if (s >= e) return 0;
        if (s <= min_idx && min_idx < e) {return min_idx + indices[min_idx];}
        uint32_t len = e - s;
        uint32_t k = log2_floor(len);
        uint32_t abs_left = s;
        uint32_t abs_right = e - (1 << k);

        if (k > 0) {
            abs_left += read_bits(st[k], s * k, k);
            abs_right += read_bits(st[k], (e - (1 << k)) * k, k);
        }

        uint32_t argmin_block;
        if (abs_right == abs_left) {
            argmin_block = abs_right;
        } else if (abs_right < s + len) {
            argmin_block = abs_left;
        } else if (e - len < abs_left) {
            argmin_block = abs_right;
        } else {
            argmin_block = (data[abs_left] <= data[abs_right]) ? abs_left : abs_right;
        }
        return (static_cast<uint64_t>(argmin_block) << block_bit_len) + indices[argmin_block];
    }
};

template<class T>
struct RMQF {
    static constexpr int B_BITS = 5;
    static constexpr int B = 1 << B_BITS; // 32, not larger!
    
    SparseTable<B_BITS, T> s;
    vector<uint32_t> m;
    vector<T> a, c;

    RMQF(vector<T> A = {}) : m(sz(A)), a(A), c(sz(A)) {
        if (A.empty()) return;
        
        vector<T> b(sz(a) / B + 1);
        vector<uint8_t> b_idx(sz(a) / B + 1);
        
        uint32_t mi = 0;
        rep(i, sz(a)) {
            int blk = i / B;
            if (i % B == 0 || a[i] < b[blk]) {
                b[blk] = a[i];
                b_idx[blk] = i % B;
            }
            
            mi <<= 1;
            while (mi && a[i] < a[i - __builtin_ctz(mi)])
                mi ^= (1u << __builtin_ctz(mi));
            m[i] = mi ^= 1; 
            c[i] = a[i - __lg(m[i])];
        }
        
        s = SparseTable<B_BITS, T>(b, b_idx);
    }

    T get(int l, int r) {
        if (r - l + 1 < B)
            return a[r - __lg(m[r] & ((1u << (r - l + 1)) - 1))];
            
        T k = min(c[r], c[l + B - 1]);
        int st_l = (l + B - 1) / B;
        int st_r = r / B - 1;
        
        if (st_l <= st_r) {
            k = min(k, a[s.get(st_l, st_r + 1)]);
        }
        return k; 
    }

    size_t getSize() const {
        size_t total_bytes = sizeof(*this);
        
        total_bytes += m.capacity() * sizeof(uint32_t);
        total_bytes += a.capacity() * sizeof(T);
        total_bytes += c.capacity() * sizeof(T);
        
        total_bytes += (s.getSize() - sizeof(s));
        
        return total_bytes;
    }
};

template<class T>
struct SparseTableRMQ {
	vector<vector<T> > s;
	SparseTableRMQ(vector<T> a = {}) : s(1, a) {
		if (!sz(a)) return;
		rep(d, __lg(sz(a))) {
			s.eb(sz(a) - (1 << d) * 2 + 1);
			rep(j, sz(s[d + 1]))
				s[d + 1][j] = min(s[d][j], s[d][j + (1 << d)]);
		}
	}
	T get(int l, int r) {
		int d = __lg(r - l + 1);
		return min(s[d][l], s[d][r - (1 << d) + 1]);
	}
};

template<class T>
struct RMQ_Alstrup {
	static constexpr int B = 32; // not larger!
	SparseTableRMQ<pair<T, uint32_t>> s;
    using pair_type = pair<T, uint32_t>;
    vector<pair_type> val_mask;
	vector<T> a;

	RMQ_Alstrup(vector<T> A = {}) : a(A), val_mask(sz(A)) {
		int nb = (sz(a) + B - 1) / B;
        // block minimums
        vector<pair<T, uint32_t>> b(nb);
        // monotone queue
		uint32_t mi = 0;
		rep(i, sz(a)) {
            // update block minimum and block minimum index
            b[i / B] = (i % B ? min(b[i / B], pair_type{a[i], i}) : pair_type{a[i], i}); 
            // move queue 1 left
			mi <<= 1;
            // pop all values larger then the current one from the queue
			while (mi && a[i] < a[i - __builtin_ctz(mi)])
				mi ^= (1u << __builtin_ctz(mi));
			// add current value to queue
            val_mask[i].second = mi ^= 1;
            // calculate smallest value in the current 32-block
            val_mask[i].first = a[i - __lg(mi)];
		}
        // recursibly construct ST
		s = SparseTableRMQ<pair<T, uint32_t>>(b);
	}
	uint32_t get(int l, int r) {
		// if the interval is small we can use montone queues
        auto [right_val, right_mask] = val_mask[r];
        if (r - l + 1 < B) {
			return r - __lg(right_mask & ((1u << (r - l + 1)) - 1));
        }

        auto [left_val, left_mask] = val_mask[l + B - 1];
        uint32_t right_idx  = r - __lg(right_mask);
        uint32_t left_idx   = l + B - 1 - __lg(left_mask);

        uint32_t k_idx = right_idx;
        auto k_val = right_val;
        if (left_val <= right_val) {
            k_idx = left_idx;
            k_val = left_val;
        }
        l = (l + B - 1) / B, r = r / B - 1;
        if (l <= r) {
            // RMQ works on pairs {value, index}. It's proven to be faster and easier to implement that way.
            auto inner_idx = s.get(l, r).second;
            if (a[inner_idx] < k_val || (a[inner_idx] == k_val && inner_idx <= k_idx)) {
                k_idx = inner_idx;
            }
        }
		return k_idx; }
};

// template<typename T>
// vector<pair<T, size_t>> indexed_vec(const vector<T>& v) {
//     vector<pair<T, size_t>> ret(v.size());
//     for (int i = 0; i < v.size(); ++i) {
//         ret[i] = {v[i], i};
//     }
//     return ret;
// }
// template<class T>
// struct RMQ_Alstrup {
// 	RMQ_Alstrup_Return_Value<pair<T, size_t>> rmq_value;
//  	RMQ_Alstrup(vector<T> A = {}) : rmq_value{indexed_vec(A)} {
// 	}
// 	size_t get(int l, int r) {
//         return rmq_value.get(l, r).second;
//     }
// };

#endif