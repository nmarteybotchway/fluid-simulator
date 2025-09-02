#pragma once
#include <vector>
#include <cassert>
#include <algorithm>

struct Grid {
    int size;                  // grid dimension (square grid)
    std::vector<float> data;   // underlying storage

    Grid(int n = 0, float init = 0.0f)
        : size(n), data(n * n, init) {}

    inline float& operator()(int x, int y) {
        assert(x >= 0 && x < size && y >= 0 && y < size && "Grid index out of bounds");
        return data[x + y * size];
    }

    inline const float& operator()(int x, int y) const {
        assert(x >= 0 && x < size && y >= 0 && y < size && "Grid index out of bounds");
        return data[x + y * size];
    }

    size_t flatSize() const { return data.size(); }

    void fill(float value) {
        std::fill(data.begin(), data.end(), value);
    }
};
