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

namespace Tools {

  /**
   * class Float1D is a proxy class that can represent, for example,
   * a column or row vector of a Float2D array, where row (sub-)arrays
   * are stored with a respective stride.
   * Besides constructor/deconstructor, the class provides overloading of
   * the []-operator, such that elements can be accessed as v[i]
   * (independent of the stride).
   * The class will never allocate separate memory for the vectors,
   * but point to the interior data structure of Float2D (or other "host"
   * data structures).
   */
  template <class T>
  class Float1D {
  private:
    int rows_;
    int stride_;
    T*  data_;

  public:
    Float1D(T* data, int rows, int stride = 1):
      rows_(rows),
      stride_(stride),
      data_(data) {}

    ~Float1D() = default;

    T& operator[](int i) { return data_[i * stride_]; }

    const T& operator[](int i) const { return data_[i * stride_]; }

    T* getData() { return data_; }

    const T* getData() const { return data_; }

    int getSize() const { return rows_; }
  };

} // namespace Tools
