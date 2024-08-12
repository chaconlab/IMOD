/***************************************************************************
                          vlvoldata.cpp  -  description
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

/**
 * This class will store the data of vlVolume in different layouts.
 * This class is a virtual base class.
 *
 * @author Sarang Lakare
 */

#include "vlvoldata.h"


/**
 * default constructor
 *
 * @param dim             dimensions of the data to be created.
 * @param dataType        type of data to be stored.
 */
vlVolData::vlVolData(const vlDim & dim, const vlDataType dataType,
                    const vlUnit & units, const vlLayoutType layout)
                    : m_dim(dim),
                    m_step(dim),
                    m_units(units),
                    m_dataType(dataType),
                    m_layout(layout),
                    m_dirty(false)
{
}

/**
 * default destructor.
 */
vlVolData::~vlVolData()
{

}

/// Returns the dimensions of the 3D data
vlDim vlVolData::dim() const
{
  return (m_dim);
}

/// Returns the stepping in voxel distance along x, y and z axis
vlStep vlVolData::stepping() const
{
  return (m_step);
}

/// Returns the actual distance between the voxels along x, y and z
vlUnit vlVolData::units() const
{
  return (m_units);
}

/// Returns the number of bits that form a voxel
uint16 vlVolData::bitsPerVoxel() const
{
  return (m_bitsPerVoxel);
}

/// Returns the number of bytes that form a voxel
uint16 vlVolData::bytesPerVoxel() const
{
  return (m_bytesPerVoxel);
}

/// Returns the total number of voxels in the volume = XDim x YDim x ZDim
uint64 vlVolData::voxelCount() const
{
  return (m_dim.x()*m_dim.y()*m_dim.z());
}

/// Returns the type of the data stored
vlDataType vlVolData::dataType() const
{
  return (m_dataType);
}

/// Returns the layout in which the data is stored
vlLayoutType vlVolData::layout() const
{
  return (m_layout);
}

/// Returns the state of the dirty flag
bool vlVolData::isDirty() const
{
  return (m_dirty);
}

/// Set/Reset the dirty flag
void vlVolData::setDirty(bool dirty)
{
  m_dirty = dirty;
}
