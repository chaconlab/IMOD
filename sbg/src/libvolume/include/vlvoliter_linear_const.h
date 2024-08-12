/***************************************************************************
                          vlvoliter_linear_const.h  -  description
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


 #ifndef _vlVolIterLinearConst_h
#define _vlVolIterLinearConst_h

#include "tripletypes.h"
#include "vloffset.h"
#include "vlneighborhood.h"
#include "vlvolume.h"
#include "vlvoliterbase.h"
#include "vlvoldata_linear.h"
#include "vlinterpolator.h"
#include "vlvoxeloperator.h"


/**
 * This class defines an iterator for "linear" data volumes. It is derived
 * from vlVolIterator.
 * I  does not provide modifications methods
 *
 * @todo Make & use local copies of m_step, m_dimLimit etc
 * @todo Replace m_rData with m_pData and check performance
 * @todo Add robust error checks - easy if moved to m_pData
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <class DataType>
class vlVolIterConst<DataType, vlLayout::Linear>
	: public vlVolIterBaseConst<DataType>
{
public:

	/// Constructor when volume pointer is provided
	vlVolIterConst<DataType, vlLayout::Linear>(const vlVolume * vol)
		: vlVolIterBaseConst<DataType>(vol),
      m_endOfNeighborhood(true),
			m_outside(false),
			m_valid(false)
	{
		if(vol->dataLayout() == vlLayout::Linear) {
    	this->m_data = dynamic_cast< vlVolDataLayout<DataType, vlLayout::Linear> * >(this->m_volData);
			if(!m_data)
      {
        return;
      }
      m_pFirstVoxel = m_pCurrVoxel = m_data->dataPtr();
			m_voxelCount = m_data->voxelCount();
			m_step = m_data->stepping();
	    m_pLastVoxel = m_pFirstVoxel+m_voxelCount-1;
			m_dim = m_data->dim();
			m_dimLimit = vlDim(m_dim.x()-1, m_dim.y()-1, m_dim.z()-1);
			m_valid = true;
      m_nativeIpolator=NULL;
    } else {
    	fprintf(stderr,"ERROR : vlVolIterConst<DataType, vlLayout::Linear> : Trying to use wrong layout.\n");
			m_data = 0L;
		}
	};

	/// Constructor when volume data pointer is provided
	vlVolIterConst<DataType, vlLayout::Linear>(const vlVolData * data)
		: vlVolIterBaseConst<DataType>(data),
      m_endOfNeighborhood(true),
			m_outside(false),
			m_valid(false)
	{
		if(data->layout() == vlLayout::Linear) {
    	m_data = dynamic_cast< vlVolDataLayout<DataType, vlLayout::Linear> * >(this->m_volData);
			if(!m_data)
				return;
	    m_pFirstVoxel = m_pCurrVoxel = m_data->dataPtr();
			m_voxelCount = m_data->voxelCount();
			m_step = m_data->stepping();
	    m_pLastVoxel = m_pFirstVoxel+m_voxelCount-1;
			m_dim = m_data->dim();
			m_dimLimit = vlDim(m_dim.x()-1, m_dim.y()-1, m_dim.z()-1);
			m_valid = true;
      m_nativeIpolator=NULL;
		} else {
			fprintf(stderr,"ERROR : vlVolIterConst<DataType, vlLayout::Linear> : Trying to use wrong layout.\n");
			m_data = 0L;
		}
	};

	/// Constructor when volume data in linear layout pointer is provided
	vlVolIterConst<DataType, vlLayout::Linear>(const vlVolDataLayout<DataType, vlLayout::Linear> * data)
		: vlVolIterBaseConst<DataType>((vlVolDataLayoutBase<DataType>*)data),
      m_endOfNeighborhood(true),
			m_outside(false),
			m_valid(false)
	{
  	m_data = const_cast<vlVolDataLayout<DataType, vlLayout::Linear> *>(data);
    m_pFirstVoxel = m_pCurrVoxel = m_data->dataPtr();
		m_voxelCount = m_data->voxelCount();
		m_step = m_data->stepping();
    m_pLastVoxel = m_pFirstVoxel+m_voxelCount-1;
  	m_dim = m_data->dim();
		m_dimLimit = vlDim(m_dim.x()-1, m_dim.y()-1, m_dim.z()-1);
		m_valid = true;
    m_nativeIpolator=NULL;
	};

  /// default destructor
  virtual ~vlVolIterConst<DataType, vlLayout::Linear>()
  {
    if(m_nativeIpolator)
    {
        //delete m_nativeIpolator;
        //m_nativeIpolator
        m_nativeIpolator->~vlInterpolator<DataType, vlLayout::Linear>() ;
    }
  };

  /// Returns the data volume Buffer
  vlVolDataLayoutBase<DataType> *getVolume() const
    {

    return this->m_volData;
  };


  ///Return the valid flag
  bool isValid()
	{
  	return (m_valid);
	};

   /// move the internal iterator to the first voxel
  bool begin() {
    m_pCurrVoxel=m_pFirstVoxel;
    return true;
  };

  /// Returns true if end of volume reached
  bool end() {
    if(m_pLastVoxel == m_pCurrVoxel) return (true);
    return (false);
  };

  /// Return the position of the internal iterator
  vlPoint3ui pos();

  /// move the iterator to voxel at newPosition
  bool moveTo(const vlPoint3ui & newPosition);

  /// move the iterator by offset
  bool moveRelative(const vlOffset & offset);

  /// Get the value of voxel at current position
  DataType get();

  /**
   * Get the value of voxel at the given position.
   * NOTE : Not implemented yet.
   * \todo Implement getValueAt(vlPoint3ui pos)
   */
  DataType getValueAt(const vlPoint3ui & pos);

  /**
   * Get the value of voxel at the given position - interpolation will be used.
   * Set the type of interpolation to use using setInterpolation().
   * NOTE : Not implmented as of this update of the documentation.
   * \todo Implement interpolation.
   */
  DataType getValueAt(const vlPoint3f & pos);

  /// move the iterator to the next position (no order assumed)
  bool operator++();

  /// move the iterator to the previous position (no order assumed)
  bool operator--();

  /// move the iterator to the next voxel - NOTE : No order assumed - fastest traversal
  bool next();

  /// move the iterator to the next voxel - WARNING : No bounds checking
  void nextNBC();

  /// move the iterator to the next voxel - move first along X, then Y and then Z
  bool nextXYZ();

  /// move the iterator to the next voxel - move first along Y, then Z and then X
  bool nextYZX();

  /// move the iterator to the next voxel - move first along Z, then X and then Y
  bool nextZXY();

  /// move the iterator to the previous voxel - NOTE : No order assumed - fastest traversal
  bool prev();

  /// move the iterator to the previous voxel - WARNING : No bounds checking
  void prevNBC();

  /// move the iterator to the previous voxel - move first along X, then Y and then Z
  bool prevXYZ();

  /// move the iterator to the previous voxel - move first along Y, then Z and then X
  bool prevYZX();

  /// move the iterator to the previous voxel - move first along Z, then X and then Y
  bool prevZXY();

  /// get value at a voxel with position delta from the current one
  DataType getRelative(const vlOffset & offset);

  /**
  * Get value at a point delta from the curret position of the iterator.
  * This will be an iterpolated value. Set the type of interpolation
  * using setIterpolation()
  * NOTE : Not implmented as of this update of the documentation.
  * \todo Implement interpolation.
  */
  DataType getRelative(const vlPoint3f & delta);

  /// Same as getRelative but with NBC = No Boundary Check, thus faster.
  DataType getRelativeNBC(const vlOffset & offset);

  /// Get the voxel along X axis
  DataType getRelativeX(const int32 offset);

  /// Get the voxel along Y axis
  DataType getRelativeY(const int32 offset);

  /// Get the voxel along Z axis
  DataType getRelativeZ(const int32 offset);

  ///Set the template to select the enviorement voxels of a voxel
  bool setNeighborhood(const vlNeighborhood & neighborhood);

  ///Return the current neighbor voxel
  DataType getNeighbor();

  ///Return the Offset of the current neighbor to the current voxel
  vlOffset getNeighborOffset();

  ///Return the id of the current neighbor
  uint16 getNeighborId();

  ///Step top the next neighbot
  bool nextNeighbor();

  ///Set the current neighbor to the first voxel of the enviorement
  bool firstNeighbor();

  ///Return true if the nextNeoighbor method has been invoked over the last neighbor
  bool endOfNeighborhood() { return (m_endOfNeighborhood); };

  /// Returns the status of outside flag
  bool outside() { return (m_outside); };

  /// Explicitly reset the outside flag
  void resetOutside() { m_outside=false; };

 /**
 *  Set the interpolator
 */
  bool setInterpolation(vlInterpolatorBase<DataType> *interP);

  /// Return a pointer to the interpolator
  vlInterpolatorBase<DataType> * getInterpolator();

  /**
  * Set the operator of the Iterator
  */
  bool setVoxelOperation(vlVoxelOpBase<DataType> *op);

  /**
  *  Return pointer to the Operator
  */
  vlVoxelOpBase<DataType> const * getVoxelOperator();

  /**
  * Return the value given by the operator over the Volume
  */
  bool getOp(vlVoxelOpValue & value);

  /// Write the volume in a file
  // bool write(char *fileName,vlPoint3f *position) ;

  /**
   * Get value at a voxel with position deltaOffset from the current one.
   */
  DataType getRelative(const int32 deltaOffset);

  /**
   * Same as getRelative(offset) but with NBC = No Boundary Check.
   */
  DataType getRelativeNBC(const int32 deltaOffset);

protected:
  /**
   * Stores the end of neighborhood marker. Set to true when
   * nextNeighbor() is called while at the last voxel in then
   * neighborhood list.
   */
  bool m_endOfNeighborhood;

  /// the total voxels in the volume
  unsigned long m_voxelCount;

  /// this will always point to the first voxel
  DataType *m_pFirstVoxel;

  /// this will always point to the current voxel
  DataType *m_pCurrVoxel;

  /// Pointer to the last voxel of the volume
  DataType *m_pLastVoxel;

  /// the volume data on which the iterator is defined
  vlVolDataLayout<DataType, vlLayout::Linear> * m_data;

  /// Store the stepping distance along the three axis
  vlStep m_step;

  /// Store the dimensions of the volume
  vlDim m_dim;

  /// Limit for dimensions
  vlDim m_dimLimit;

  /// The type of interpolation to use
  int m_interpolation;

  /// Set to true if get*() accessed outside
  bool m_outside;

  /// True if iterator is valid - initialized properly
  bool m_valid;

  /// Stores the neighborhood being used
  vlNeighborhood m_neighborhood;

  /// Stores pointer to the current neighbor
  std::vector<vlOffset>::const_iterator m_currNeighbor;

  /// Stores pointer to the end of offset list
  std::vector<vlOffset>::const_iterator m_endNeighbor;

  /// Stores the ID of the neighbor.
  uint16 m_neighborId;

  ///Pointer to the Interpolator
   vlInterpolator<DataType, vlLayout::Linear> *m_nativeIpolator;

  /// Stores the pointer to native voxel operator
  vlVoxelOperator<DataType, vlLayout::Linear> *m_nativeVoxelOp;


};


/**
 * Return the position of the iterator
 *
 * @param DataType
 * @return
 */
template <class DataType>
vlPoint3ui vlVolIterConst<DataType, vlLayout::Linear>::pos()
{
  uint32 x_pos, y_pos, z_pos, offset;
  offset = m_pCurrVoxel - m_pFirstVoxel;
  z_pos = offset/m_step.z();
  y_pos = (offset%m_step.z())/m_step.y();
  x_pos = (offset%m_step.z())%m_step.y();
  return vlPoint3ui(x_pos, y_pos, z_pos);
}

/**
 * moves the iterator to newPosition. Also checks if the given position is inside
 * the volume. If outside, keeps the current position and returns false.
 *
 * @param newPosition     the new position to move the iterator to.
 * @return true if the move was successful, false otherwise.
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::moveTo(const vlPoint3ui & newPosition)
{
  // check if the new position is within limits
  if( newPosition.x() > m_dimLimit.x()
     || newPosition.y() > m_dimLimit.y()
     || newPosition.z() > m_dimLimit.z()) {
    return (false);
  }

  // compute offset and set current voxel pointer
  m_pCurrVoxel = m_pFirstVoxel + newPosition.x() + newPosition.y()*m_step.y()
    + newPosition.z()*m_step.z();

  return (true);
}


/**
 * moves the iterator by offset. Also checks if the given position is inside
 * the volume. If outside, keeps the current position and returns false.
 *
 * @param offset is the offset by which to move the iterator.
 * @return true if the move was successful, false otherwise.
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::moveRelative(const vlOffset & offset)
{
  // get the current position
  vlPoint3ui newPos(pos());

  // compute the new position
  if((int32)newPos.x()+offset.x() >= 0)
    newPos.x(newPos.x()+offset.x());
  else
    return false;
  if((int32)newPos.y()+offset.y() >= 0)
    newPos.y(newPos.y()+offset.y());
  else
    return false;
  if((int32)newPos.z()+offset.z() >= 0)
    newPos.z(newPos.z()+offset.z());
  else
    return false;

  // check if the new position is within limits
  if( newPos.x() > m_dimLimit.x()
     || newPos.y() > m_dimLimit.y()
     || newPos.z() > m_dimLimit.z()) {
    return false;
  }

  // compute offset and set current voxel pointer
  m_pCurrVoxel = m_pFirstVoxel + newPos.x() + newPos.y()*m_step.y()
    + newPos.z()*m_step.z();

  return true;
}


/**
 * Get the value of voxel at current position
 *
 * @param DataType
 * @return
 */
template <class DataType>
inline DataType vlVolIterConst<DataType, vlLayout::Linear>::get()
{
  return (*m_pCurrVoxel);
}


/**
 * move the iterator to the next position (no order assumed)
 *
 * @param DataType
 * @return
 */
template <class DataType>
inline bool vlVolIterConst<DataType, vlLayout::Linear>::operator++()
{
  if(m_pCurrVoxel < m_pLastVoxel) {
    ++m_pCurrVoxel;
    return (true);
  }

  return (false);
}


/**
 * move the iterator to the previous position (no order assumed)
 *
 * @param DataType
 * @return
 */
template <class DataType>
inline bool vlVolIterConst<DataType, vlLayout::Linear>::operator--()
{
  // check if going outside the volume
  if(m_pCurrVoxel > m_pFirstVoxel)
    return (false);

  --m_pCurrVoxel;

  return (true);
}



/**
 * move the iterator to the next voxel in the dataset.
 * WARNING : No order is assumed for the traversal.. i.e. its not guarenteed
 * that the next voxel will be neighbor of the current one and so on..
 * This is simply the fastest way to access all voxels in the dataset.
 *
 * @param DataType
 * @return
 * @see nextXYZ()
 * @see nextYZX()
 * @see nextZXY()
 */
template <class DataType>
inline bool vlVolIterConst<DataType, vlLayout::Linear>::next()
{
  if(m_pCurrVoxel < m_pLastVoxel) {
    ++m_pCurrVoxel;
    return (true);
  }

  return (false);
}

/**
 * move the iterator to the next voxel in the dataset w/o checking for any bounds.
 * WARNING : Using this w/o thinking is dangerous. You have been warned!
 * This is simply the fastest way to move to the next voxel. You need to do
 * bounds-checking yourself, else the pointer can easily go outside the data range.
 *
 * @see next()
 */
template <class DataType>
inline void vlVolIterConst<DataType, vlLayout::Linear>::nextNBC()
{
	++m_pCurrVoxel;
}


/**
 * move the iterator to the next voxel - move first along X, then Y and then Z
 *
 * @param DataType
 * @return
 */
template <class DataType>
inline bool vlVolIterConst<DataType, vlLayout::Linear>::nextXYZ()
{
  // check if going outside the volume
  if(m_pCurrVoxel < m_pLastVoxel) {
    m_pCurrVoxel+=m_step.x();
    return (true);
  }

  return (false);
}


/**
 * move the iterator to the next voxel - move first along Y, then Z and then X
 *
 * @param DataType
 * @return
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::nextYZX()
{
  uint32 x_pos, y_pos, z_pos;

  y_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())/m_step.y();
  if((int)(y_pos) == (int)(m_dim.y()-1)) {
    // cannot move further along y...
    y_pos = 0;
  } else {
    m_pCurrVoxel+=m_step.y();
    return (true);
  }

  x_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())%m_step.y();
  z_pos = (m_pCurrVoxel-m_pFirstVoxel)/m_step.z();

  if((int)(z_pos) == (int)(m_dim.z()-1)) {
    // cannot move further along z...
    z_pos = 0;
  } else {
    ++z_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }

  if((int)(x_pos) == (int)(m_dim.x()-1)) {
    // cannot move further atall! we are at the end!
    return (false);
  } else {
    ++x_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }
}


/**
 * move the iterator to the next voxel - move first along Z, then X and then Y
 *
 * @param DataType
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::nextZXY()
{
  uint32 x_pos, y_pos, z_pos;

  z_pos = (m_pCurrVoxel-m_pFirstVoxel)/m_step.z();

  if((int)(z_pos) == (int)(m_dim.z()-1)) {
    // cannot move further along z...
    z_pos = 0;
  } else {
    m_pCurrVoxel+=m_step.z();
    return (true);
  }

  x_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())%m_step.y();
  y_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())/m_step.y();

  if((int)(x_pos) == (int)(m_dim.x()-1)) {
    // cannot move further along x...
    x_pos = 0;
  } else {
    ++x_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }

  if((int)(y_pos) == (int)(m_dim.y()-1)) {
    // cannot move further at all..
    return (false);
  } else {
    ++y_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }
}


/**
 * move the iterator to the previous voxel in the dataset.
 * WARNING : No order is assumed for the traversal.. i.e. its not guarenteed
 * that the previous voxel will be neighbor of the current one and so on..
 * This is simply the fastest way to access all voxels in the dataset.
 *
 * @param DataType
 * @return
 * @see nextXYZ()
 * @see nextYZX()
 * @see nextZXY()
 */
template <class DataType>
inline bool vlVolIterConst<DataType, vlLayout::Linear>::prev()
{
  // check if going outside the volume
  if(m_pCurrVoxel > m_pFirstVoxel) {
    return false;
  }

  --m_pCurrVoxel;
  return true;
}

/// move the iterator to the previous voxel - WARNING : No bounds checking
template <class DataType>
inline void vlVolIterConst<DataType, vlLayout::Linear>::prevNBC()
{
	--m_pCurrVoxel;
}


/**
 * move the iterator to the previous voxel - move first along X, then Y and then Z
 *
 * @param DataType
 * @return
 */
template <class DataType>
inline bool vlVolIterConst<DataType, vlLayout::Linear>::prevXYZ()
{
  // check if going outside the volume
  if(m_pCurrVoxel > m_pFirstVoxel) {
    return false;
  }

  m_pCurrVoxel-=m_step.x();
  return true;

}


/**
 * move the iterator to the previous voxel - move first along Y, then Z and then X
 *
 * @param DataType
 * @return
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::prevYZX()
{
  uint32 x_pos, y_pos, z_pos;

  y_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())/m_step.y();
  if(y_pos == 0) {
    // cannot move below y==0...
    y_pos = m_dim.y()-1;
  } else {
    m_pCurrVoxel-=m_step.y();
    return (true);
  }

  x_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())%m_step.y();
  z_pos = (m_pCurrVoxel-m_pFirstVoxel)/m_step.z();

  if(z_pos == 0) {
    // cannot move below z==0...
    z_pos = m_dim.z()-1;
  } else {
    --z_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }

  if(x_pos == 0) {
    // cannot go any more behind.. we are at 0,0,0!!
    return (false);
  } else {
    --x_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }
}


/**
 * move the iterator to the previous voxel - move first along Z, then X and then Y
 *
 * @param DataType
 * @return
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::prevZXY()
{
  uint32 x_pos, y_pos, z_pos;

  z_pos = (m_pCurrVoxel-m_pFirstVoxel)/m_step.z();

  if(z_pos == 0) {
    // cannot move below z==0...
    z_pos = m_dim.z()-1;
  } else {
    m_pCurrVoxel-=m_step.z();
    return (true);
  }

  x_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())%m_step.y();
  y_pos = ((m_pCurrVoxel-m_pFirstVoxel)%m_step.z())/m_step.y();

  if(x_pos == 0) {
    // cannot move below x==0...
    x_pos = m_dim.x()-1;
  } else {
    --x_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }

  if(y_pos == 0) {
    // cannot move at all! we are at 0,0,0!!
    return (false);
  } else {
    --y_pos;
    return (moveTo(vlPoint3ui(x_pos, y_pos, z_pos)));
  }
}


/**
 * move the iterator a relative offset
 *
 * @param DataType
 * @return
 */
template <class DataType>
inline DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelative(const int32 deltaOffset)
{
  // check if going outside the volume
  if(0 <= (m_pCurrVoxel-m_pFirstVoxel)+deltaOffset < m_voxelCount) {
		m_outside = false;
    return (*(m_pCurrVoxel+deltaOffset));
  }

//  cout << "WARNING : Peeking outside the volume!" << endl;
	m_outside = true;
  return (0);
}

/**
* Move the iterator a relative offset, It does not check if the bounds are not checked.
* Youi must be carefoul when using this method
*
* @param DataType
* @return
*/
template <class DataType>
inline DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelativeNBC(const int32 deltaOffset)
{
	return (*(m_pCurrVoxel+deltaOffset));
}


/**
* Return Voxel in a relative position from the current Voxel
*
* @param offset: Relative position
*/
template <class DataType>
inline DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelative(const vlOffset & offset)
{
	vlPoint3ui newPos(pos());
  if(offset.x()+newPos.x() <= m_dimLimit.x() && offset.x()+newPos.x() >= 0
		 && offset.y()+newPos.y() <= m_dimLimit.y() && offset.y()+newPos.y() >= 0
		 && offset.z()+newPos.z() <= m_dimLimit.z() && offset.z()+newPos.z() >= 0	) {

	 	int32 deltaOffset = offset.x()*m_step.x() + offset.y()*m_step.y() + offset.z()*m_step.z();
		m_outside = false;
    return (*(m_pCurrVoxel+deltaOffset));
	}

//  cout << "WARNING : Peeking outside the volume!: " << m_dimLimit.x() << endl;

	m_outside = true;
  return (0);
}


/**
* Return Voxel in a relative position from the current Voxel
*  There is not bounds checking
* @param offset: Relative position
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelativeNBC(const vlOffset & offset)
{
 	int32 deltaOffset = offset.x()*m_step.x() + offset.y()*m_step.y() + offset.z()*m_step.z();

  return (*(m_pCurrVoxel+deltaOffset));
}


/**
* Return Voxel in a relative position from the current Voxel
* The offset is only applied along the X axe
*
* @param: offset: Offset along the X axe
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelativeX(const int32 offset)
{
  int32 deltaOffset = offset*m_step.x();
  // check if going outside the volume
  if((int)((m_pCurrVoxel-m_pFirstVoxel)+deltaOffset) < (int)(m_voxelCount)) {
		m_outside = false;
    return (*(m_pCurrVoxel+deltaOffset));
  }

//  cout << "WARNING : Peeking outside the volume!" << endl;
	m_outside = true;
  return (0);
}

/**
* Return Voxel in a relative position from the current Voxel
* The offset is only applied along the Y axe
*
* @param: offset: Offset along the Y axe
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelativeY(const int32 offset)
{
  int32 deltaOffset = offset*m_step.y();
  // check if going outside the volume
  if((int)((m_pCurrVoxel-m_pFirstVoxel)+deltaOffset) < (int)(m_voxelCount)) {
		m_outside = false;
    return (*(m_pCurrVoxel+deltaOffset));
  }

// std::cout << "WARNING : Peeking outside the volume!" << std::endl;
	m_outside = true;
  return (0);
}

/**
* Return Voxel in a relative position from the current Voxel
* The offset is only applied along the Z axe
*
* @param: offset: Offset along the Z axe
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelativeZ(const int32 offset)
{
  int32 deltaOffset = offset*m_step.z();
  // check if going outside the volume
  if((int)((m_pCurrVoxel-m_pFirstVoxel)+deltaOffset) < (int)(m_voxelCount)) {
		m_outside = false;
    return (*(m_pCurrVoxel+deltaOffset));
  }

//  cout << "WARNING : Peeking outside the volume!" << endl;
	m_outside = true;
  return (0);
}



/**
* Get the value of voxel at the given position
*
* @param: pos: Position of the returned voxel
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getValueAt(const vlPoint3ui & pos)
{
  // bounds checking
  if(pos.x() > m_dimLimit.x() || pos.y() > m_dimLimit.y() || pos.z() > m_dimLimit.z()) {
    m_outside = true;
    return (0);
  }
  // get offset to voxel
  int32 offset = pos.x()*m_step.x() + pos.y()*m_step.y() + pos.z()*m_step.z();
  m_outside = false;
  return (*(m_pFirstVoxel+offset));
}

/**
* Get the value of voxel at the given position
* An interpolation from the adjacent voxels is made to compute
* the value of the real position
*
* @param: pos: Position of the returned voxel
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getValueAt(const vlPoint3f & pos)
{
   // bounds checking
  if(pos.x() > m_dimLimit.x() || pos.y() > m_dimLimit.y() || pos.z() > m_dimLimit.z()) {
     m_outside = true;
     return (0);
   }

  if(m_nativeIpolator)
    return(m_nativeIpolator->getValueAt(*this, pos));
  else
	  fprintf(stderr,"ERROR : getValueAt(..) : No valid interpolator set. Use setInterpolation(...).\n");

  return (0);
}

/**
* Get the value of voxel at the given offset from the current voxel
* An interpolation from the adjacent voxels is made to compute
* the value of the real position
*
* @param: pos: offset from the current voxel
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getRelative(const vlPoint3f & delta)
{

  vlPoint3ui oldPos = pos();
  vlPoint3ui newpos(oldPos.x()+(int)delta.x(), oldPos.y()+(int)delta.y(), oldPos.z()+(int)delta.z());
  if(!moveTo(newpos)) {
   m_outside = true;
   fprintf(stderr,"ERROR : getRelative moving iterator outside the volume [ %d %d %d ]\n", newpos.x(),newpos.y(),newpos.z());
     return (0);
   }

   vlPoint3f newDelta(delta.x()-(int)delta.x(), delta.y()-(int)delta.y(), delta.z()-(int)delta.z());


   DataType value;
   if(m_nativeIpolator) {
     value = m_nativeIpolator->getValueAtOffset(*this, newDelta);
   } else {
	   fprintf(stderr,"ERROR : getRelative(..) : No valid interpolator set. Use setInterpolation(...).\n");
     return (0);
   }

   moveTo(oldPos);
   return (value);
}

/**
* Set the set of voxels around the current voxel that will from its neighborhood
*
* @param neighborhood: Object that describes the offset of the neighbors ofthe current voxel
*/
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::setNeighborhood(const vlNeighborhood & neighborhood)
{
  m_neighborhood = neighborhood;
  m_currNeighbor = neighborhood.getOffsets().begin();
  m_endNeighbor = neighborhood.getOffsets().end();
  m_neighborId = 0;
  if(neighborhood.getOffsets().size())
    m_endOfNeighborhood = false;
  else
    m_endOfNeighborhood = true;
  return (true);
}

/**
* Return the current neighbor of the current voxel
*/
template <class DataType>
DataType vlVolIterConst<DataType, vlLayout::Linear>::getNeighbor()
{
  if(m_neighborhood.getOffsets().empty()) {
    return (0);
  }
  return (getRelative(*m_currNeighbor));
}

/**
* Return the offset of the current neighbor from the current voxel
*/
template <class DataType>
vlOffset vlVolIterConst<DataType, vlLayout::Linear>::getNeighborOffset()
{
  return (*m_currNeighbor);
}

/**
* Return the id of the current neighbor
*/
template <class DataType>
uint16 vlVolIterConst<DataType, vlLayout::Linear>::getNeighborId()
{
  return (m_neighborId);
}

/**
* Step to the next neighbor
*/
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::nextNeighbor()
{
  if(m_neighborhood.getOffsets().empty()) {
    return (false);
  }
  if(m_neighborId == (uint16)(m_neighborhood.getOffsets().size()-1)) {
    // at the last voxel.. cannot go ahead!
    // set end of neighborhood to true
    m_endOfNeighborhood = true;
    return (false);
  }
  ++m_neighborId;
  ++m_currNeighbor;
  return (true);
}


/**
* Set the current neighbor to the first one in the neighborhood
*/
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::firstNeighbor()
{
  if(m_neighborhood.getOffsets().empty()) {
    m_endOfNeighborhood = true;
    return (false);
  }
  m_currNeighbor = m_neighborhood.getOffsets().begin();
  m_neighborId = 0;
  m_endOfNeighborhood = false;
  return (true);
}


/**
* Set the interpolation to apply for compute values in real positions
*
* @param *iter: Pointer to object Interpolation that implements the interpolation
*/
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::setInterpolation(vlInterpolatorBase<DataType> *interP)
{
  m_nativeIpolator=(vlInterpolator<DataType, vlLayout::Linear>*)interP;
  return true ;
}


/**
* Return the interpolation object
*/
template <class DataType>
vlInterpolatorBase<DataType> * vlVolIterConst<DataType, vlLayout::Linear>::getInterpolator()
  {
    if(m_nativeIpolator)
      return (m_nativeIpolator);
    else
      return (0L);
  };


/**
* Set the internal operator
*
* @param op: Object that implements an operation over the data buffer
*/
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::setVoxelOperation(vlVoxelOpBase<DataType> *op)
{
  m_nativeVoxelOp = (vlVoxelOperator<DataType, vlLayout::Linear>*)op;
  return true;
}

/**
 * Return pointer to the object Operator
 */
template <class DataType>
vlVoxelOpBase<DataType> const * vlVolIterConst<DataType, vlLayout::Linear>::getVoxelOperator()
{
  if(m_nativeVoxelOp)
    return (m_nativeVoxelOp);
  else
    return(0L);
}

/**
 * Return the computed value from the operator over the data
 *
 * @param value: Returned value
 */
template <class DataType>
bool vlVolIterConst<DataType, vlLayout::Linear>::getOp(vlVoxelOpValue & value)
{
      if(m_nativeVoxelOp)
        return (m_nativeVoxelOp->getValue(*this, value));
      else
        return(false);
}



/**
 * Write volume to a file
 *
 * @param fileName: name of the file
 * @param Position: Reference position of the volume in the real space
 * @return true if operation is sucesssful
 */

// removed by Pablo avoiding gcc errors....
//template <class DataType>
//bool vlVolIterConst<DataType, vlLayout::Linear>::write(char *fileName,vlPoint3f *position)
//{
//  int count=0;
//  FILE *fout;
//  vlVolDataLayoutBase<DataType> *vol;
//
//  fout=fopen(fileName,"w");
//  if(fout==NULL)
//  {
//	  fprintf(stderr,"Error: Can't open file\n");
//    return false;
//  }
//
//  vol=getVolume();
//
//  fprintf(fout," %f %f %f %f %d %d %d\n", vol->units().x(), position->x(), position->y(), position->z(), vol->dim().x(), vol->dim().y(), vol->dim().z());
//  fprintf(fout,"\n");
//
//  m_pCurrVoxel = m_pFirstVoxel;
//
//  DataType a;
//  while(m_pLastVoxel >= m_pCurrVoxel)
//  {
//
//    if((count+1)%10==0)
//      fprintf(fout," %10.6f \n",(*m_pCurrVoxel));
//
//    else
//      fprintf(fout," %10.6f ",(*m_pCurrVoxel));
//
//      m_pCurrVoxel+=m_step.x();
//
//  }
//
//
//
//  fclose(fout);
//
//
//return(true);
//}




#endif // _vlVolIterLinearConst_h
