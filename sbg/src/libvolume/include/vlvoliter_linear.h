/***************************************************************************
                          vlvoliter_linear.h  -  description
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


 #ifndef _vlVolIterLinear_h
#define _vlVolIterLinear_h

#include "vlvoliter_linear_const.h"

/**
 * This class defines an iterator for "linear" data volumes. It is derived
 * from vlVolIterator.
 * Establece la implementacion de los metodos de modificacion de datos del
 * volumen para un layout lineal. Los metodos de no modificacion
 * los hereda de vlVolIterConst
 *
 * @todo Make & use local copies of m_step, m_dimLimit etc
 * @todo Replace m_rData with m_pData and check performance
 * @todo Add robust error checks - easy if moved to m_pData
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <typename DataType>
class vlVolIter<DataType, vlLayout::Linear>
	: public vlVolIterConst<DataType, vlLayout::Linear>,
		public vlVolIterBase<DataType> {
public:

	/// Constructor when volume pointer is provided
	vlVolIter<DataType, vlLayout::Linear>(vlVolume * vol)
		: vlVolIterConst<DataType, vlLayout::Linear>(vol)
	{
	};

	/// Constructor when volume data pointer is provided
	vlVolIter<DataType, vlLayout::Linear>(vlVolData * data)
		: vlVolIterConst<DataType, vlLayout::Linear>(data)
	{
	};

	/// Constructor when volume data in linear layout pointer is provided
	vlVolIter<DataType, vlLayout::Linear>(vlVolDataLayout<DataType, vlLayout::Linear> * data)
		: vlVolIterConst<DataType, vlLayout::Linear>(data)
	{
	};

  /// default destructor
  virtual ~vlVolIter<DataType, vlLayout::Linear>() { };

  /// Set the value of voxel at current position
  bool set(const DataType & value);

  /// set value at a voxel offset from the current one
  bool setRelative(const vlOffset & offset, const DataType & value);

	/// Same as getRelative but with NBC = No Boundary Check, thus faster.
	bool setRelativeNBC(const vlOffset & offset, const DataType & value);

	/// Get the voxel along X axis
	bool setRelativeX(const int32 offset, const DataType & value);

	/// Get the voxel along Y axis
	bool setRelativeY(const int32 offset, const DataType & value);

	/// Get the voxel along Z axis
	bool setRelativeZ(const int32 offset, const DataType & value);

  bool setNeighbor(const DataType & value);

  /**
   * set value at a voxel deltaOffset from the current position.
   * This function is available only to Linear layout.
   */
  bool setRelative(const int32 deltaOffset, const DataType & value);
};

// inline implementation of the API

/**
 * Set the value of voxel at current position
 *
 * @param DataType
 * @param value
 * @return
 */
template <typename DataType>
inline bool vlVolIter<DataType, vlLayout::Linear>::set(const DataType & value)
{
  *(this->m_pCurrVoxel) = value;
  return true;
}


template <typename DataType>
inline bool vlVolIter<DataType, vlLayout::Linear>::setRelative(const int32 deltaOffset, const DataType & value)
{
  // check if going outside the volume
  if(0 <= (this->m_pCurrVoxel-this->m_pFirstVoxel)+deltaOffset < this->m_voxelCount) {
    *(this->m_pCurrVoxel+deltaOffset) = value;
    return (true);
  }

//  cout << "WARNING : Writing outside the volume!" << endl;
  return (false);
}


template <typename DataType>
inline bool vlVolIter<DataType, vlLayout::Linear>::setRelative(const vlOffset & offset, const DataType & value)
{
 int32 deltaOffset = offset.x()*this->m_step.x() + offset.y()*this->m_step.y() + offset.z()*this->m_step.z();
  // check if going outside the volume
  if((int)((this->m_pCurrVoxel-this->m_pFirstVoxel)+deltaOffset) <= (int)(this->m_voxelCount)) {
    *(this->m_pCurrVoxel+deltaOffset) = value;
    return (true);
  }

//  cout << "WARNING : Writing outside the volume!" << endl;
  return (false);
}

template <typename DataType>
bool vlVolIter<DataType, vlLayout::Linear>::setRelativeNBC(const vlOffset & offset, const DataType & value)
{
 int32 deltaOffset = offset.x()*this->m_step.x() + offset.y()*this->m_step.y() + offset.z()*this->m_step.z();
  // check if going outside the volume
  if(((int)(this->m_pCurrVoxel-this->m_pFirstVoxel)+deltaOffset) <= (int)(this->m_voxelCount)) {
    *(this->m_pCurrVoxel+deltaOffset) = value;
    return (true);
  }

//  cout << "WARNING : Writing outside the volume!" << endl;
  return (false);
}

template <typename DataType>
bool vlVolIter<DataType, vlLayout::Linear>::setRelativeX(const int32 offset, const DataType & value)
{
 int32 deltaOffset = offset*this->m_step.x();
  // check if going outside the volume
  if(int((this->m_pCurrVoxel-this->m_pFirstVoxel)+deltaOffset) <= (int)(this->m_voxelCount)) {
    *(this->m_pCurrVoxel+deltaOffset) = value;
    return (true);
  }

//  cout << "WARNING : Writing outside the volume!" << endl;
  return (false);
}

template <typename DataType>
bool vlVolIter<DataType, vlLayout::Linear>::setRelativeY(const int32 offset, const DataType & value)
{
 int32 deltaOffset = offset*this->m_step.y();
  // check if going outside the volume
  if((int)((this->m_pCurrVoxel-this->m_pFirstVoxel)+deltaOffset) <= (int)(this->m_voxelCount)) {
    *(this->m_pCurrVoxel+deltaOffset) = value;
    return (true);
  }

//  cout << "WARNING : Writing outside the volume!" << endl;
  return (false);
}

template <typename DataType>
bool vlVolIter<DataType, vlLayout::Linear>::setRelativeZ(const int32 offset, const DataType & value)
{
 int32 deltaOffset = offset*this->m_step.z();
  // check if going outside the volume
  if((int)((this->m_pCurrVoxel-this->m_pFirstVoxel)+deltaOffset) <= (int)(this->m_voxelCount)) {
    *(this->m_pCurrVoxel+deltaOffset) = value;
    return (true);
  }

//  cout << "WARNING : Writing outside the volume!" << endl;
  return (false);
}

template <class DataType>
bool vlVolIter<DataType, vlLayout::Linear>::setNeighbor(const DataType & value)
{
  if(this->m_neighborhood.getOffsets().empty()) {
    return (false);
  }

  return (setRelative(*(this->m_currNeighbor), value));
}

#endif // _vlVolIterLinear_h
