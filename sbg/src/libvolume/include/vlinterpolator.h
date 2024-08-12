/***************************************************************************
                          vlinterpolator.h  -  description
                             -------------------
    begin                : Tue Apr 13 2004
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


 #ifndef _vlInterpolator_h
#define _vlInterpolator_h

#include <string>
#include "vlenums.h"
#include "vlvoliterbase.h"


/**
 * This class defines an interpolator. An interpolator is a special class that
 * does interpolation. Interpolation is needed to find value at a non-grid
 * location. This interpolation is used by iterators to get value at non-grid
 * location when requested by user. Every interpolation has a generic implementation
 * which will work with any volume data layout. Each layout can optionally re-implement
 * an interpolation to get speed-up.
 *
 * This is a super base class of an interpolator object which is independent of the
 * datatype and the layout type. This class can be used to pass around interpolator
 * object pointers. To implement an interpolator, subclass from vlInterpolator.
 *
 * @author Sarang Lakare <sarang#users.sourceforge.net>
 * @see vlInterpolator
 */
class vlInterpolatorSuperBase
{
public:
  /// Gives the interpolation type of this interpolation
  virtual vlInterpolationType type() = 0;

  /// Gives the name of the interpolation
  virtual std::string name() = 0;

  /// Returns the layout for which this interpolator is implemented
  virtual vlLayoutType layout() = 0;
};


/**
 * Template of superclass vlInterpolatorSuperBase with type specified.
 * @see vlInterpolatorSuperBase
 */
template <typename DataType>
class vlInterpolatorBase : public vlInterpolatorSuperBase
{

};


/**
 * Template of class vlInterpolatorBase with layout specified.
 * @see vlInterpolatorSuperBase
 */
template <typename DataType, vlLayoutType Layout>
class vlInterpolator : public vlInterpolatorBase<DataType>
{
public:
  /**
   * Returns the value at a given position.
   *
   * @param iter: Iterator over the volume (at the end the iterator is set in the position?).
   * @param position: Position to check.
   * @param check: If check is true, check for errors/bounds etc.
   */
  virtual DataType getValueAt(vlVolIterConst<DataType, Layout> & iter,
                              const vlPoint3f & position, bool check=true) = 0;

  /**
   * Returns the value at a given offset of the position of the iterator.
   *
   * @param iter: Iterator over the volume.
   * @param offset: offset to apply.
   * @param check: If check is true, check for errors/bounds etc.
   */
  virtual DataType getValueAtOffset(vlVolIterConst<DataType, Layout> & iter,
                                    const vlPoint3f & offset, bool check=true) = 0;

};



#endif // _vlInterpolator_
