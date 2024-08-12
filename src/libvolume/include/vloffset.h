/***************************************************************************
                          vloffset.h  -  description
                             -------------------
    begin                : Mon Apr 12 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef _vlOffset_h
#define _vlOffset_h

#include <stdio.h>

#include <vector>


#include "vltriple.h"

/**
 * This class defines an offset. Offset is the relative location of one voxel
 * with respect to another.
 *
 * The offset can be accessed using either x(), y(), z() and t() functions or
 * using the overloaded [] operator. The offset is optimized for 3D & 4D access.
 *
 */
class vlOffset : public vlTriple<int16> {

public:
  /// Constructor in which the number of dimensions is specified
  vlOffset(const uint16 dim=3)
    : vlTriple<int16>(0,0,0),
      m_dummy(0), m_dim(dim)
  {
    if(dim > 3) {
      m_offset.resize(dim-3);
    }
  };

  /// Constructor in which a triple is given
  vlOffset(const vlTriple<int16> & triple)
    : vlTriple<int16>(triple),
      m_dummy(0), m_dim(3)
  {

  };

  /// Special constructor for 3D
  vlOffset(const int16 x, const int16 y, const int16 z)
    : vlTriple<int16>(x,y,z),
      m_dummy(0), m_dim(3)
  {

  };

  /// Special constructor for 4D
  vlOffset(const int16 x, const int16 y, const int16 z, const int16 t)
    : vlTriple<int16>(x,y,z),
      m_dummy(0), m_dim(4)
  {
    m_offset.resize(1);
    m_offset[0] = t;
  };

  ~vlOffset(){};

  /**
   * Returns the t component of the offset. Also available via offset[3]
   */
  int16 t() const
  {
    return (m_offset[0]);
  };

  /// Returns the offset at the given dimension
  int16 operator[](const uint16 dim) const
  {
    switch(dim) {
    case 0:
      return (m_x);
      break;
    case 1:
      return (m_y);
      break;
    case 2:
      return (m_z);
      break;
    default:
      if(dim < m_dim)
        return (m_offset[dim-3]);
      fprintf(stderr,"WARNING : Accessing offset dimension which does not exist.\n");
      return 0;
      break;
    }
  };

  /// Returns the offset at the given dimension
  int16 & operator[](const uint16 dim)
  {
    switch(dim) {
    case 0:
      return (m_x);
      break;
    case 1:
      return (m_y);
      break;
    case 2:
      return (m_z);
      break;
    default:
      if(dim < m_dim)
        return (m_offset[dim-3]);
      fprintf(stderr,"WARNING : Accessing offset dimension which does not exist.\n");
      return m_dummy;
      break;
    }
  };


  /// Returns the dimensionality of the offset
  uint16 dim() const { return m_dim; };

  /// overriding << to enable writing offset to a stream
  friend std::ostream & operator << (std::ostream & os, vlOffset const * const o) {
    return (os << *o);
  }

  /// overriding << to enable writing offset to a stream
  friend std::ostream & operator << (std::ostream & os, vlOffset const & o) {
    uint16 i=o.m_offset.size(), j(1);
    os << "(";
    os << o.m_x;
    if(o.m_dim > 1)
      os << "," << o.m_y;
    if(o.m_dim > 2)
      os << "," << o.m_z;
    if(o.m_dim > 3) {
      os << "," << o.m_offset[0];
      while(i!=j)
        os << "," << o.m_offset[j++];
    }
    os << ")";
    return os;
  }

private:
  /// store the offsets greater than the first 3 dims
  std::vector<int16> m_offset;

  /**
   * Some dummy which will be returned when offset of dimension greater
   * than the real is accessed.
   */
  int16 m_dummy;

  /// Dimensions of the offset
  uint16 m_dim;
};

#endif // _vlOffset_h
