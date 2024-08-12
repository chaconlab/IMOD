/***************************************************************************
                          vlvoldata.h  -  description
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

#ifndef _vlVolData_h
#define _vlVolData_h


#include "tripletypes.h"
#include "vlenums.h"



/**
 * This class will store the 3D data of a vlVolume in different layouts.
 * This class is a base class. vlVolDataLayout inherits this class
 * and defines additional API for each data layout. To add a new data
 * layout, subclass vlVolDataLayout.
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 * @see vlVolDataLayout
 */
class vlVolData {
// Member functions
public:
  /// default destructor
  virtual ~vlVolData();

  /// Returns the dimensions of the 3D data
  virtual vlDim dim() const;

  /// Returns the stepping in voxel distance along x, y and z axis
  virtual vlStep stepping() const;

  /// Returns the actual distance between the voxels along x, y and z
  virtual vlUnit units() const;

  /// Returns the number of bits that form a voxel
  virtual uint16 bitsPerVoxel() const;

  /// Returns the number of bytes that form a voxel
  virtual uint16 bytesPerVoxel() const;

  /// Returns the total number of voxels in the volume = XDim x YDim x ZDim
  virtual uint64 voxelCount() const;

  /// Returns the type of the data stored
  virtual vlDataType dataType() const;

  /// Returns the type name of the stored data
  virtual const char * typeInfo() const = 0;

  /// Returns the layout in which the data is stored
  virtual vlLayoutType layout() const;

  /// Returns the state of the dirty flag
  virtual bool isDirty() const;

  /// Set/Reset the dirty flag
  virtual void setDirty(bool dirty);

  /// Returns true if the data is valid
  virtual bool isValid() const = 0;

  /// Clears the volume data with the given data
  virtual bool clear(const uint8 data = 0x00) = 0;

  /**
   * Returns the pointer to the voxel at the given position. This returns a void pointer,
   * so make sure you cast it to the correct type. To avoid type-conflicts, "use iterators".
   *
   * @param position        location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  virtual void * getVoxelVoidPtr(const vlPoint3ui & position) const = 0;

protected:
  /// Constructor - only for derived classes
  vlVolData(const vlDim & dim,
            const vlDataType dataType,
            const vlUnit & units,
            const vlLayoutType layout);

  /// dimension for the data along x, y and z
  vlDim m_dim;

  /// step distance for the data along x, y and z
  vlStep m_step;

  /// Unit size of each voxel
  vlUnit m_units;

  /// data type of the stored data
  vlDataType m_dataType;

  /// The layout in which the data is stored
  vlLayoutType m_layout;

  /// The number of bits and bytes each voxel consumes
  uint16 m_bitsPerVoxel, m_bytesPerVoxel;

  /// Dirty flag - true indicates volume has changed
  bool m_dirty;
};

#endif // _vlVolData_h
