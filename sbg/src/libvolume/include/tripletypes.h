/***************************************************************************
                          tripletypes.h  -  description
                             -------------------
    begin                : Thu Apr 1 2004
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


#ifndef _tripleTypes_h       // prevent multiple declarations
#define _tripleTypes_h

#include "vltriple.h"
#include <ostream>


typedef vlTriple<int16> vlPoint3i;
typedef vlTriple<uint16> vlPoint3ui;
typedef vlTriple<float> vlPoint3f;
typedef vlTriple<double> vlPoint3d;



/**
 * This datatype defines the dimensions of a volume. Each dimension is
 * stored as 16 bit unsigned int. Thus, the range for each dimension is
 * 0 - 2^16-1, i.e. 0 - 65535.
 */
class vlDim : public vlTriple<uint16> {
public:
  /// Default constructor
  vlDim(uint16 xdim=0, uint16 ydim=0, uint16 zdim=0)
    : vlTriple<uint16>(xdim, ydim, zdim) { };

  /// overriding << to enable writing triple to a stream
  friend std::ostream & operator << (std::ostream & os, vlDim const * const dim) {
    return (os << *dim);
  }

  /// overriding << to enable writing triple to a stream
  friend std::ostream & operator << (std::ostream & os, vlDim const & dim) {
    return os<<dim.m_x<<"x"<<dim.m_y<<"x"<<dim.m_z;
  }

};



/**
 * This datatype defines the stepping distance in a volume. Each step distance
 * is the distance between consecutive voxels in voxel units assuming :
 * - The data is stored linearly
 * - The data is on a rectilinear grid
 */
class vlStep : public vlTriple<uint32> {
public:
  /// Default constructor
  vlStep(uint32 xdim=0, uint32 ydim=0, uint32 zdim=0)
    : vlTriple<uint32>(xdim, ydim, zdim) { };

  /// Constructing the step when dimensions are given
  vlStep(const vlDim & dim)
    : vlTriple<uint32>(1, dim.x(), dim.x()*dim.y()) { };

  /// overriding << to enable writing triple to a stream
  friend std::ostream & operator << (std::ostream & os, vlStep const * const step) {
    return (os << *step);
  }

  /// overriding << to enable writing triple to a stream
  friend std::ostream & operator << (std::ostream & os, vlStep const & step) {
    return os<<step.m_x<<"x"<<step.m_y<<"x"<<step.m_z;
  }
};


/**
 * This data structure will store the unit distance between voxels.
 */
class vlUnit : public vlTriple<float> {

public:
  /**
   * The type of unit that is stored in this object.
   */
  enum UnitType {
    Meter,
    MilliMeter,
    MicroMeter,
    NoUnitType
  };

	/// Default constructor
	vlUnit(float x=1.0, float y=1.0, float z=1.0, UnitType type=NoUnitType)
		: vlTriple<float>(x, y, z) { };

	/// Returns of the type of the unit
	UnitType type() const { return (m_type); };

protected:
	/// The type of unit represented - meter, micro-meter etc.
	UnitType m_type;

};

#endif // _tripleTypes_h
