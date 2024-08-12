/***************************************************************************
                          vlvoldata_linear.h  -  description
                             -------------------
    begin                : Fri Apr 2 2004
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

#ifndef _vlVolDataLinear_h
#define _vlVolDataLinear_h

#include <stdlib.h>
#include <typeinfo>
#include <string.h>

#include "vlvoldatalayoutbase.h"
#include "vlmacros.h"



template <typename A, vlLayoutType Layout >
class vlVolIter;

template <typename A, vlLayoutType Layout >
class vlVolIterConst;

/**
 * This class implements a linear storage layout for the 3D volumetric data.
 * The data voxels are stored as a linear array. The data is stored in the
 * X-Y-Z order.
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <typename DataType>
class vlVolDataLayout<DataType, vlLayout::Linear>
  : public vlVolDataLayoutBase<DataType> {
public:


  // declare the Iter class to be friend
  friend class vlVolIterConst<DataType, vlLayout::Linear>;
  friend class vlVolIter<DataType, vlLayout::Linear>;

  /// Default constructor
  vlVolDataLayout<DataType, vlLayout::Linear>(const vlDim & dim, const vlUnit & units = vlUnit(1.0, 1.0, 1.0));

  /// Default destructor
  virtual ~vlVolDataLayout<DataType, vlLayout::Linear>();

  /// Get the voxel value at the given position
  DataType getVoxel(const vlPoint3ui & position) const;

  /// Get the voxel value at the given position
  DataType getVoxel(const vlPoint3f & position) const;

  /// SPECIAL.Get the voxel value at the given position in the buffer in fast way
  DataType getVoxelSpecial(const int position) const;

  /// SPECIAL.Set the voxel at given position to the given value in a fast way
  bool setVoxelSpecial(vlPoint3ui const & position, DataType const voxel);

  /// SPECIAL.Set the voxel at given position to the given value in a fast way
  bool setVoxelSpecial(int const & position, DataType const voxel);

  /// Set the voxel at given position to the given value
  bool setVoxel(vlPoint3ui const & position, DataType const voxel);

  /// Increase the value of the voxels around the given position with interpolations
  bool incVoxel(const vlPoint3f & position, const DataType voxel);

  /// Returns true if the data is valid
  bool isValid() const { return (m_pData != 0L); };

  /// Clears the volume data with the given data
  bool clear(const uint8 data = 0x00);

  /**
   * Returns the pointer to the voxel at the given position. This returns a void pointer,
   * so make sure you cast it to the correct type. To avoid type-conflicts, "use iterators".
   *
   * @param position        location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  void * getVoxelVoidPtr(const vlPoint3ui & position) const;

  /**
   * Returns the pointer to the voxel at the given position. Do not use
   * this unless it is absolutely necessary. Use iterators instead.
   *
   * @param position        location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  DataType * getVoxelPtr(const vlPoint3ui & position) const;

protected:
  /// Return a const pointer to the stored data
  DataType * dataPtr() const { return m_pData; };


private:
  /// pointer to the linear data
  DataType * m_pData;

  /// this is dimension-1 in each direction. for speeding up bounds checking - save one subtraction
  vlDim m_dimLimit;


};


/**
 * default constructor.
 *
 * @param DataType        the type of data to be stored.
 * @param dim             dimensions of the data.
 * @param dataType        type of the data to be stored.
 */
template <class DataType>
vlVolDataLayout<DataType, vlLayout::Linear>::vlVolDataLayout(const vlDim & dim, const vlUnit & units)
  : vlVolDataLayoutBase<DataType>(dim, UnsignedInt8, units, vlLayout::Linear)
{
  // allocate memory for data storage
  m_pData = (DataType *) calloc(dim.x()*dim.y()*dim.z(), sizeof(DataType));

  if(!m_pData) {
	  fprintf(stderr,"vlVolDataLinear::vlVolDataLinear() : Unable to allocate memory for 3D data : [ %d %d %d ]\n",dim.x(),dim.y(),dim.z());
    this->m_dim.x(0);this->m_dim.y(0);this->m_dim.z(0);
    this->m_bytesPerVoxel = 0;
    this->m_bitsPerVoxel = 0;
  } else {
    DataType temp;
    getVariableDataType(temp, this->m_dataType);
    this->m_bytesPerVoxel = sizeof(DataType);
    this->m_bitsPerVoxel = this->m_bytesPerVoxel*8;
    this->m_dimLimit = vlDim(dim.x()-1, dim.y()-1, dim.z()-1);
  }

}

/**
 * default destructor
 *
 * @param DataType        the type of data to be stored.
 */
template <class DataType>
vlVolDataLayout<DataType, vlLayout::Linear>::~vlVolDataLayout()
{
  free(m_pData);
}


/**
 * get the voxel at the given position
 *
 * @param DataType
 * @param position        position of the voxel (x,y,z)
 * @return voxel value at the given position
 */
template <class DataType>
inline DataType vlVolDataLayout<DataType, vlLayout::Linear>::getVoxel(const vlPoint3ui & position) const
{
  if (position.x() > m_dimLimit.x() || position.y() > m_dimLimit.y() || position.z() > m_dimLimit.z()) {
    // given position is outside the dimensions of the data
    return 0;
    // THROW EXCEPTION or WARN
  }

  // return the voxel data
  return *(this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x()));
}



/**
 * SPECIAL. get the voxel at int position
 *
 * @param DataType
 * @param position        position in the buffer
 * @return voxel value at the given position
 */
template <class DataType>
inline DataType vlVolDataLayout<DataType, vlLayout::Linear>::getVoxelSpecial(const int position) const
{
  // return the voxel data
  return *(m_pData + position);
}

/**
 * SPECIAL.sets the value of the voxel at the given position
 *
 * @param DataType
 * @param position        position of the voxel
 * @param voxel           voxel value
 * @return always return true. It does not check if the position is valid
 */
template <class DataType>
inline bool vlVolDataLayout<DataType, vlLayout::Linear>::setVoxelSpecial(const vlPoint3ui & position, const DataType voxel)
{
  // set the voxel data
  *(this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x())) = voxel;
 this->m_dirty = true;
  // return success
  return true;
}

/**
 * SPECIAL.sets the value of the voxel at the given position
 *
 * @param DataType
 * @param position        Secuencial position in the buffer
 * @param voxel           voxel value
 * @return always return true. It does not check if the position is valid
 */
template <class DataType>
inline bool vlVolDataLayout<DataType, vlLayout::Linear>::setVoxelSpecial(const int &position, const DataType voxel)
{
// set the voxel data
  *(m_pData + position) = voxel;
  this->m_dirty = true;
  // return success
  return true;
}



/**
 * get the voxel at the given position
 * IMPLEMENTED ONLY FOR TYPE FLOAT
 *
 * @param DataType
 * @param position        position of the voxel (x,y,z)
 * @return voxel value at the given position
 */
template <class DataType>
inline DataType vlVolDataLayout<DataType, vlLayout::Linear>::getVoxel(const vlPoint3f & position) const
{
  if (position.x() > m_dimLimit.x() || position.y() > m_dimLimit.y() || position.z() > m_dimLimit.z()) {
    // given position is outside the dimensions of the data
    return 0;
    // THROW EXCEPTION or WARN
  }

  fprintf(stderr,"ERROR : Voxel at a float position not implemented yet.\n");

  // return the voxel data
  return (0);
}

template <> inline float vlVolDataLayout<float, vlLayout::Linear>::getVoxel(const vlPoint3f & position) const
{
  if (position.x() >= m_dimLimit.x() || position.y() >= m_dimLimit.y() || position.z() >= m_dimLimit.z()) {
    // given position is outside the dimensions of the data
    return 0;
    // THROW EXCEPTION or WARN
  }

  int x0,x1,y0,y1,z0,z1;
  float a,b,c;
  float f111,f112,f121,f122,f211,f212,f221,f222;
  float aux;

  x0=floor(position.x());
  x1=x0+1;
  if(x1>=m_dimLimit.x())
   x1--;

  y0=floor(position.y());
  y1=y0+1;
  if(y1>=m_dimLimit.y())
   y1--;

  z0=floor(position.z());
  z1=z0+1;
  if(z1>=m_dimLimit.z())
   z1--;


  a=position.x()-x0; b=position.y()-y0; c=position.z()-z0;

  f111=*(m_pData + (z0*m_step.z()+y0*m_step.y()+x0));
  f112=*(m_pData + (z1*m_step.z()+y0*m_step.y()+x0));
  f121=*(m_pData + (z0*m_step.z()+y1*m_step.y()+x0));
  f122=*(m_pData + (z1*m_step.z()+y1*m_step.y()+x0));
  f211=*(m_pData + (z0*m_step.z()+y0*m_step.y()+x1));
  f212=*(m_pData + (z1*m_step.z()+y0*m_step.y()+x1));
  f221=*(m_pData + (z0*m_step.z()+y1*m_step.y()+x1));
  f222=*(m_pData + (z1*m_step.z()+y1*m_step.y()+x1));

  aux=((1-a)*((1-b)*((1-c)*f111+c*f112)+b*((1-c)*f121+c*f122))+
           a *((1-b)*((1-c)*f211+c*f212)+b*((1-c)*f221+c*f222)) );

  // return the voxel data
  return (aux);
}

/**
 * sets the value of the voxel at the given position
 *
 * @param DataType
 * @param position        position of the voxel
 * @param voxel           voxel value
 * @return true if the voxel is set, false otherwise
 */
template <class DataType>
inline bool vlVolDataLayout<DataType, vlLayout::Linear>::setVoxel(const vlPoint3ui & position, const DataType voxel)
{
  if (position.x() > m_dimLimit.x() || position.y() > m_dimLimit.y() || position.z() > m_dimLimit.z()) {
    // given position is outside the dimensions of the data
    return false;
    // THROW EXCEPTION or WARN
  }

  // set the voxel data
  *(this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x())) = voxel;

  this->m_dirty = true;

  // return success
  return true;
}

/**
 * Increases the value of the voxel around the given position with interpolations
 * ONLY IMPLEMENTED FOR FLOAT TYPE
 *
 * @param DataType
 * @param position        position between the voxels
 * @param voxel           voxel value
 * @return true if the voxel is set, false otherwise
 */
 template <class DataType>
 inline bool vlVolDataLayout<DataType, vlLayout::Linear>::incVoxel(const vlPoint3f & position, DataType voxel)
 {
   return true;
 }

/**
* Increase the values of the voxels around the given position.
* It interpolates the value between the adjacent voxels
*
* @param position: Real Position for the interpolation
* @param voxel: Value to be interpolated
*/
template <> inline bool vlVolDataLayout<float, vlLayout::Linear>::incVoxel(const vlPoint3f & position, float voxel)
{
  if (position.x() > m_dimLimit.x() || position.y() > m_dimLimit.y() || position.z() > m_dimLimit.z()
      || position.x() < 0 || position.y() < 0 || position.z() < 0) {
    // given position is outside the dimensions of the data
    return false;
    // THROW EXCEPTION or WARN
  }
  int x0,x1,y0,y1,z0,z1;
  float a,b,c,ab,ab1,a1b,a1b1;
  x0=(int)floor(position.x());
  x1=x0+1;
  y0=(int)floor(position.y());
  y1=y0+1;
  z0=(int)floor(position.z());
  z1=z0+1;


  a = x1 - position.x();
  b = y1 - position.y();
  c = z1 - position.z();
  ab = a * b;
  ab1 = a * ( 1 - b );
  a1b = ( 1 - a ) * b;
  a1b1 = ( 1 - a ) * ( 1 - b );
  a = ( 1 - c );

  // set the voxel data
  if(z0>0)
  {
    if(y0>0)
    {
      if(x0>0)
        *(m_pData + (z0*m_step.z()+y0*m_step.y()+x0)) += ab*c*voxel;
      if(x1<m_dimLimit.x())
        *(m_pData + (z0*m_step.z()+y0*m_step.y()+x1)) += a1b*c*voxel;
    }
    if(y1<m_dimLimit.y())
    {
      if(x0>0)
        *(m_pData + (z0*m_step.z()+y1*m_step.y()+x0)) += ab1*c*voxel;
      if(x1<m_dimLimit.x())
        *(m_pData + (z0*m_step.z()+y1*m_step.y()+x1)) += a1b1*c*voxel;
    }
  }

  if(z1<m_dimLimit.z())
  {
    if(y0>0)
    {
      if(x0>0)
        *(m_pData + (z1*m_step.z()+y0*m_step.y()+x0)) += ab*a*voxel;
      if(x1<m_dimLimit.x())
        *(m_pData + (z1*m_step.z()+y0*m_step.y()+x1)) += a1b*a*voxel;
    }
    if(y1<m_dimLimit.y())
    {
      if(x0>0)
        *(m_pData + (z1*m_step.z()+y1*m_step.y()+x0)) += ab1*a*voxel;
      if(x1<m_dimLimit.x())
        *(m_pData + (z1*m_step.z()+y1*m_step.y()+x1)) += a1b1*a*voxel;
    }
  }

  m_dirty = true;

  // return success
  return true;
}



template <class DataType>
inline bool vlVolDataLayout<DataType, vlLayout::Linear>::clear(const uint8 data)
{
  if(memset(this->m_pData, data, sizeof(DataType)*this->m_step.z()*this->m_dim.z())) {
    this->m_dirty = true;
    return (true);
  }else
    return (false);
}

template <class DataType>
inline void * vlVolDataLayout<DataType, vlLayout::Linear>::getVoxelVoidPtr(const vlPoint3ui & position) const
{
  if (position.x() > m_dimLimit.x() || position.y() > m_dimLimit.y() || position.z() > m_dimLimit.z()) {
    // given position is outside the dimensions of the data
    return 0;
    // THROW EXCEPTION or WARN
  }

  // return the voxel data
  return ((void*)(this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x())));
}

template <class DataType>
inline DataType * vlVolDataLayout<DataType, vlLayout::Linear>::getVoxelPtr(const vlPoint3ui & position) const
{
  if (position.x() > m_dimLimit.x() || position.y() > m_dimLimit.y() || position.z() > m_dimLimit.z()) {
    // given position is outside the dimensions of the data
    return 0;
    // THROW EXCEPTION or WARN
  }

  // return the voxel data
  return (this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x()));
}





#endif // _vlVolDataLinear_h
