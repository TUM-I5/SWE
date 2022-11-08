/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema
 * @author Sebastian Rettenberger
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#pragma once

#include <cstdio>

#include "Float1D.hpp"

namespace Tools {

  /**
   * class Float2D is a very basic helper class to deal with 2D float arrays:
   * indices represent columns (1st index, "horizontal"/x-coordinate) and
   * rows (2nd index, "vertical"/y-coordinate) of a 2D grid;
   * values are sequentially ordered in memory using "column major" order.
   * Besides constructor/deconstructor, the class provides overloading of
   * the []-operator, such that elements can be accessed as a[i][j].
   */
  template <class T>
  class Float2D {
  private:
    int rows_;
    int cols_;

    T* data_;

    bool allocateMemory_;

  public:
    /**
     * Constructor:
     * takes size of the 2D array as parameters and creates a respective Float2D object;
     * allocates memory for the array, but does not initialise value.
     * @param cols number of columns (i.e., elements in horizontal direction)
     * @param rows rumber of rows (i.e., elements in vertical directions)
     */
    Float2D(int cols, int rows, bool allocateMemory = true):
      rows_(rows),
      cols_(cols),
      data_(nullptr),
      allocateMemory_(allocateMemory) {

      if (allocateMemory_) {
        data_ = new T[rows * cols];
      }
    }

    /**
     * Constructor:
     * takes size of the 2D array as parameters and creates a respective Float2D object;
     * this constructor does not allocate memory for the array, but uses the allocated memory
     * provided via the respective variable #data
     * @param cols number of columns (i.e., elements in horizontal direction)
     * @param rows rumber of rows (i.e., elements in vertical directions)
     * @param data pointer to a suitably allocated region of memory to be used for thew array elements
     */
    Float2D(int cols, int rows, T* data):
      rows_(rows),
      cols_(cols),
      data_(data),
      allocateMemory_(false) {}

    /**
     * Constructor:
     * takes size of the 2D array as parameters and creates a respective Float2D object;
     * this constructor does not allocate memory for the array, but uses the allocated memory
     * provided via the respective variable data
     * @param cols number of columns (i.e., elements in horizontal direction)
     * @param rows rumber of rows (i.e., elements in vertical directions)
     * @param data pointer to a suitably allocated region of memory to be used for thew array elements
     */
    Float2D(Float2D<T>& data, bool shallowCopy):
      rows_(data.rows_),
      cols_(data.cols_),
      allocateMemory_(!shallowCopy) {

      if (shallowCopy) {
        data_           = data.data_;
        allocateMemory_ = false;
      } else {
        data_ = new T[rows_ * cols_];
        for (int i = 0; i < rows_ * cols_; i++) {
          data_[i] = data.data_[i];
        }
        allocateMemory_ = true;
      }
    }

    ~Float2D() {
      if (allocateMemory_) {
        delete[] data_;
      }
    }

    T* operator[](int i) { return (data_ + (rows_ * i)); }

    const T* operator[](int i) const { return (data_ + (rows_ * i)); }

    T* getData() { return data_; }

    int getRows() const { return rows_; }

    int getCols() const { return cols_; }

    Float1D<T> getColProxy(int i) { return Float1D<T>(data_ + (rows_ * i), rows_); }

    Float1D<T> getRowProxy(int j) { return Float1D<T>(data_ + j, cols_, rows_); }

    static void toString(const Float2D<T>& toPrint) {
      for (int row = 0; row < toPrint.rows_; row++) {
        for (int col = 0; col < toPrint.cols_; col++) {
          printf("%f ", toPrint[col][row]);
        }
        printf("\n");
      }
    }
  };

} // namespace Tools
