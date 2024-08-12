/***************************************************************************
                          vlvoldatalayoutbase.h  -  description
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

 #ifndef _vlVolDataLayout_h
#define _vlVolDataLayout_h

#include <typeinfo>


#include "vlvoldata.h"
#include "vlenums.h"

/**
 * This class will store the 3D data of vlVolume in different layouts.
 * This class is a base class. For every data layout, define a class which will
 * subclass this one and implement all the virtual methods.
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <class DataType>
class vlVolDataLayoutBase : public vlVolData {

public:
  /// Default constructor
  vlVolDataLayoutBase(const vlDim & dim, const vlDataType dataType,
                      const vlUnit & units, const vlLayoutType layout)
    : vlVolData(dim, dataType, units, layout) { };

  /// Virtual destructor
  virtual ~vlVolDataLayoutBase() { };

  /// Returns the type name of the datatype - compilar dependant
  const char * typeInfo() const { return (typeid(DataType).name()); };

  /// Returns true if the data is valid
  virtual bool isValid() const = 0;

  /// Clears the volume data with the given data
  virtual bool clear(const uint8 data) = 0;

  /// Get the voxel value at the given position
  virtual DataType getVoxel(const vlPoint3ui & position) const = 0;

  /// Get the voxel value at the given position
  virtual DataType getVoxel(const vlPoint3f & position) const = 0;

  /// Get the voxel value at the given position in the buffer
   virtual DataType getVoxelSpecial(const int position) const = 0;

  /// Set the voxel value at the given position to 'voxel'
  virtual bool setVoxelSpecial(const vlPoint3ui & position, const DataType voxel) = 0;

  /// Set the voxel value at the given position to 'voxel'
  virtual bool setVoxelSpecial(const int & position, const DataType voxel) = 0;

  /// Set the voxel value at the given position to 'voxel'
  virtual bool setVoxel(const vlPoint3ui & position, const DataType voxel) = 0;
  /// Increase the value of the voxels around the given position with interpolations
  virtual bool incVoxel(const vlPoint3f & position, const DataType voxel) = 0;

  /**
   * Returns the pointer to the voxel at the given position. Do not use
   * this unless it is absolutely necessary. Use iterators instead.
   *
   * @param position        location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  virtual DataType * getVoxelPtr(const vlPoint3ui & position) const = 0;

  /**
   * Returns the pointer to the voxel at the given position. This returns a void pointer,
   * so make sure you cast it to the correct type. To avoid type-conflicts, "use iterators".
   *
   * @param position        location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  void * getVoxelVoidPtr(const vlPoint3ui & position) const
  { return ((void*)(getVoxelPtr(position))); };

};

/**
 * Class Template for the Data Container Class.
 * This class inherits from the Virtual Class vlVolDataLayoutBase.
 * It is the base class for all classes thar stores Data Buffers of any kind
 * and in any type of Layout
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 * @see vlVolDataLayoutBase
 */
template <typename DataType, vlLayoutType Layout>
class vlVolDataLayout : public vlVolDataLayoutBase<DataType> {
public:
  /// Default constructor
  vlVolDataLayout(const vlDim & dim, const vlDataType dataType,
                  const vlUnit & units)
    : vlVolData(dim, dataType, units, Layout) { };

  /// Virtual destructor
  virtual ~vlVolDataLayout() { };

  /// Returns true if the data is valid
  bool isValid() const;

  /// Clears the volume data with the given data
  bool clear(const uint8 data);

  /// Get the voxel value at the given position
  DataType getVoxel(const vlPoint3ui & position) const;

  /// Get the voxel value at the given position
  DataType getVoxel(const vlPoint3f & position) const;

  /// Set the voxel value at the given position to 'voxel'
  bool setVoxel(const vlPoint3ui & position, const DataType voxel);

  /// Increase the value of the voxels around the given position with interpolations
  bool incVoxel(const vlPoint3f & position, const DataType voxel) ;

  /**
   * Returns the pointer to the voxel at the given position. Do not use
   * this unless it is absolutely necessary. Use iterators instead.
   *
   * @param position        location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  DataType * getVoxelPtr(const vlPoint3ui & position) const;
};

#endif // _vlVolDataLayoutBase_h

