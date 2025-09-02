#pragma once
#include <vector>
#include <cassert>
#include <algorithm>

/**
 * @struct Grid
 * @brief Encapsulates a 2D square grid of floats with convenient access and initialization.
 *
 * Provides bounds-checked access via operator(), fill operation, and flat size retrieval.
 */
struct Grid {
    int size;               ///< Dimension of the square grid (size x size)
    std::vector<float> data;///< Underlying 1D storage for grid values

    /**
     * @brief Constructs a square grid with optional initial value.
     * @param n Number of cells along one side
     * @param init Initial value for all cells (default 0.0f)
     */
    explicit Grid(const int n = 0, const float init = 0.0f)
        : size(n), data(n * n, init) {}

    /**
     * @brief Access a cell for reading or writing.
     * @param x X-coordinate
     * @param y Y-coordinate
     * @return Reference to the value at (x, y)
     */
    inline float& operator()(int x, int y) {
        assert(x >= 0 && x < size && y >= 0 && y < size && "Grid index out of bounds");
        return data[x + y * size];
    }


    /**
     * @brief Const access to a cell for reading.
     * @param x X-coordinate
     * @param y Y-coordinate
     * @return Const reference to the value at (x, y)
     */
    inline const float& operator()(int x, int y) const {
        assert(x >= 0 && x < size && y >= 0 && y < size && "Grid index out of bounds");
        return data[x + y * size];
    }

    /**
     * @brief Returns the total number of elements in the underlying vector.
     * @return Flat size (size*size)
     */
    [[nodiscard]] size_t flatSize() const { return data.size(); }

    /**
     * @brief Fills the grid with a constant value.
     * @param value Value to assign to every cell
     */
    void fill(const float value) {
        std::fill(data.begin(), data.end(), value);
    }
};
