/***************************************************************************
                          vlvoxelop_centraldiff.h  -  description
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


 #ifndef _vlVoxelOpCentralDiff_h
#define _vlVoxelOpCentralDiff_h

#include "vlvoxeloperator.h"


/**
 * This class performs central difference gradient calculation over a iterator of a Volume.
 * The return value is a double.
 *
 * @todo Specialize this class for different datatypes so that
 * the casting to double can be removed.
 */
template <typename DataType, vlLayoutType Layout>
class vlVoxelOpCentralDiff : public vlVoxelOperator<DataType, Layout>
{
public:
  vlVoxelOpType type() { return vlVoxelOp::CentralDiff; };

  std::string name() { return ("CentralDiff"); };

  /// Returns the layout for which this voxel operator is implemented
  vlLayoutType layout() { return Layout; };



  /**
   * Gives the value of type "DataType". This type of value is not computed hence
   * this function always returns false. In short, do not use this for this operator.
   * ???
   */
  bool getValue(vlVolIterConst<DataType, Layout> & iter, DataType & value)
  { return (false); };

  /**
   * Performs voxel operation at the voxel pointed to by the iterator and returns
   * value at the given position. The value returned is of type double. So make sure
   * you do value.dl() to get the value.
   */
  bool getValue(vlVolIterConst<DataType, Layout> & iter, vlVoxelOpValue & value)
  {
    double dx = (double)(iter.getRelativeX(1)-iter.getRelativeX(-1));
    double dy = (double)(iter.getRelativeY(1)-iter.getRelativeY(-1));
    double dz = (double)(iter.getRelativeZ(1)-iter.getRelativeZ(-1));

    value.setDl(sqrt(dx*dx+dy*dy+dz*dz));

    return (true);
  };



};

#endif // _vlVoxelOpCentralDiff_h
