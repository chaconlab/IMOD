/***************************************************************************
                          vlvoxeloperator.h  -  description
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



 #ifndef _vlVoxelOp_h
#define _vlVoxelOp_h


#include "vlvoliterbase.h"



/**
 * This class defines a voxel operator. A voxel operator is a special class
 * that performs certain operation at a given voxel. An example operation could
 * be computing gradient at given voxel. These voxel operators are be used by
 * iterators to perform voxel operations when requested by the user.
 *
 * Every voxel operator will have a generic implementation. Every layout can optionally
 * override the generic implementation and implement the voxel operation for itself.
 *
 * This is a super base class of a voxel operator object. To implement a voxel
 * operator, inherit your class from vlVoxelOp. This class is useful for passing
 * voxel operator objects around in a layout and datatype independent manner.
 *
 * @author Sarang Lakare <sarang#users.sourceforge.net>
 * @see vlVoxelOperator
 */
class vlVoxelOpSuperBase
{
public:
  /// Gives the operation type of this voxel operation
  virtual vlVoxelOpType type() = 0;

  /// Gives the name of the voxel operation
  virtual std::string name() = 0;

  /// Returns the layout for which this voxel operator is implemented
  virtual vlLayoutType layout() = 0;


};

/**
* Basic Class for Operator Objects
*/
template <typename DataType>
class vlVoxelOpBase : public vlVoxelOpSuperBase
{

};

/**
* Class Template for Operator Classesr. All the Operator object must
* impelements this class
*/
template <typename DataType, vlLayoutType Layout>
class vlVoxelOperator : public vlVoxelOpBase<DataType>
{
public:
  /**
   * Performs voxel operation at the voxel pointed to by the iterator and returns
   * value at the given position.
   * Used when the result is of the same type that the volume data
   */
  virtual bool getValue(vlVolIterConst<DataType, Layout> & iter, DataType & value) = 0;

  /**
   * Performs voxel operation at the voxel pointed to by the iterator and returns
   * value at the given position.
   * Used when the result is not of the same type that the volume data
   */
  virtual bool getValue(vlVolIterConst<DataType, Layout> & iter, vlVoxelOpValue & value) = 0;

};

#endif // _vlVoxelOp_h
