/***************************************************************************
                          vlcolor.h  -  description
                             -------------------
    begin                : Wed Apr 14 2004
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


#ifndef _vlColor_h       // prevent multiple declarations
#define _vlColor_h

#include <stdio.h>

#include "vltriple.h"

typedef enum _ColorType{
  	RGB,
		HSV,
		NoColorType
  } ColorType;

/**
 * Class to store color with 3 components (R,G and B)
 */
template <class T>
class vlColor3 : public vlTriple<T> {
public:
  /**
   * Enum to define what type of color it is : RGB, HSV etc.
   */



  /// Default constructor
  vlColor3<T>(T x=0.0, T y=0.0, T z=0.0, ColorType type=NoColorType)
		: vlTriple<T>(x, y, z),
			m_type(type) {  };

  /// Default destructor
  ~vlColor3<T>() {};

  /// Returns the type of the color : RGB, CMY etc.
	ColorType type() { return (m_type); };

protected:
	/// The type of color stored - RGB, CMY etc.
	ColorType m_type;
};

typedef vlColor3<float> vlColor3f;
typedef vlColor3<double> vlColor3d;
typedef vlColor3<int8> vlColor3b;
typedef vlColor3<uint8> vlColor3ub;
typedef vlColor3<int16> vlColor3s;
typedef vlColor3<uint16> vlColor3us;
typedef vlColor3<int32> vlColor3i;
typedef vlColor3<uint32> vlColor3ui;
typedef vlColor3<int64> vlColor3l;
typedef vlColor3<uint64> vlColor3ul;




#endif // _vlColor_h
