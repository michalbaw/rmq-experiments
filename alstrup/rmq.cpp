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
    using indexed_pair = pair<T, uint32_t>;
    using val_idx_mask_t = tuple<T, uint32_t, uint32_t>;
    SparseTableRMQ<indexed_pair> s; // sparse table
    vector<val_idx_mask_t> val_idx_mask;
	vector<indexed_pair> a;
	RMQ_Alstrup(vector<T> A = {}) : a(sz(A)), val_idx_mask(sz(A)) {
		rep(i, sz(a)) {
            a[i] = indexed_pair{A[i], i};
        }
        // preprocess block minimum values
		vector<indexed_pair> b(sz(a) / B + 1);
		// monotone queue
		uint32_t mi = 0;
		rep(i, sz(a)) {
			// calculate block minimum value
			b[i / B] = (i % B ? min(b[i / B], a[i]) : a[i]);
			// shift queue one left
			mi <<= 1;
			// need to pop all values larger then the current one
			while (mi && a[i] < a[i - __builtin_ctz(mi)])
				mi ^= (1u << __builtin_ctz(mi));
			// add the rightmost element to the queue
            mi ^= 1;
            auto [val, idx] = a[i - __lg(mi)];
			val_idx_mask[i] = val_idx_mask_t{val, idx, mi};
		}
		// recursivly construct RMQ on blocks
		s = SparseTableRMQ(b);
	}
	uint32_t get(int l, int r) {
		// if we are in the range of monotone queues, use them
		auto [rval, ridx, rmask] = val_idx_mask[r];
        if (r - l + 1 < B) { 
			return r - __lg(rmask & ((1u << (r - l + 1)) - 1));
		}
		// Take the min value to be min of right end 32-block and left end 32-block
		auto [lval, lidx, lmask] = val_idx_mask[l + B - 1];
        indexed_pair r_val_idx = indexed_pair{rval, ridx};
        indexed_pair l_val_idx = indexed_pair{lval, lidx};
        auto k = min(l_val_idx, r_val_idx);
		// compute block indicies to query recusively
		l = (l + B - 1) / B, r = r / B - 1;
		// only query if necessary
		if (l <= r) k = min(k, s.get(l, r));
		return k.second;
    }
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