/***************************************************************************
                          vlvolume.cpp  -  description
                             -------------------
    begin                : Mon Apr 5 2004
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

//
#include <stdio.h>

#include <string.h>
#include "vlmacros.h"
#include "vlvoldata_linear.h"
#include "vlvolume.h"



/*Template de asignacion de memoria para el Layout*/
template <class DataType>
vlVolData * vlVolume::getLayoutT(DataType & dummy, const vlDim & dim, const vlUnit & units)
{
   vlVolDataLayout<DataType,vlLayout::Linear>* layout = new vlVolDataLayout<DataType,vlLayout::Linear>(dim, units);
   return (layout);
}



/**
 * Default constructor for vlVolume.
 *
 * @param dim             dimensions of the volume to create (default = [1,1,1]).
 * @param dataType        data type of the volume data required.
 * @param layout          specify the layout in which the volume data is to be stored (defailt = Linear).
 */
vlVolume::vlVolume(const vlDim & dim, const vlDataType dataType, const vlUnit & units)
{

  m_pVolData = 0L;
  vlVolData *res=NULL;
  callFunctionOnDataType(dataType, res, getLayoutT, dim, units);
  m_pVolData=res;



}

/**
 * Copy Constructor
 * @param old: Reference Volume to be copied
 */
vlVolume::vlVolume(vlVolume *old)
{
  int i,num_voxel;
  vlPoint3f oldPosition, newPosition;

  m_pVolData = 0L;
  vlVolData *res=NULL;
  callFunctionOnDataType(old->dataType(), res, getLayoutT, old->dim(), old->units());
  m_pVolData=res;

  old->getPosition(&oldPosition);
  newPosition.x(oldPosition.x());
  newPosition.y(oldPosition.y());
  newPosition.z(oldPosition.z());

  setPosition(newPosition);

  num_voxel=m_pVolData->voxelCount();
  for(i=0;i<num_voxel;i++)
  {
    setVoxelSpecial(i,old->getVoxelSpecial(i));
  }

}


/**
 * Constructor to created a padded volume from a previous volume
 */
vlVolume::vlVolume(vlVolume *old,vlDim margin)
{
  vlPoint3f oldPosition, newPosition;

  vlDim total(1,1,1);
  total.x(old->dim().x()+(margin.x()*2));
  total.y(old->dim().y()+(margin.y()*2));
  total.z(old->dim().z()+(margin.z()*2));

  m_pVolData = 0L;
  vlVolData *res=NULL;
  callFunctionOnDataType(old->dataType(), res, getLayoutT, total, old->units());
  m_pVolData=res;

  old->getPosition(&oldPosition);
  newPosition.x(oldPosition.x()-(float)(margin.x()* old->units().x() ));
  newPosition.y(oldPosition.y()-(float)(margin.y()* old->units().y() ));
  newPosition.z(oldPosition.z()-(float)(margin.z()* old->units().z() ));

  setPosition(newPosition);
}

/// Constructor fto created a incremented volume from a previous volume
vlVolume::vlVolume(vlVolume *old,vlDim inc, int dummy)
{
  vlPoint3f aux;

  dummy++;
  vlDim total;
  total.x(old->dim().x()+inc.x());
  total.y(old->dim().y()+inc.y());
  total.z(old->dim().z()+inc.z());

  m_pVolData = 0L;
  vlVolData *res=NULL;
  callFunctionOnDataType(old->dataType(), res, getLayoutT, total, old->units());
  m_pVolData=res;
  old->getPosition(&aux);
  setPosition(aux);
}


/**
 * Default destructor.
 */
vlVolume::~vlVolume()
{
  // release data
	if(m_pVolData)
	  delete m_pVolData;
}


bool vlVolume::clear(const uint8 data)
{
	if(m_pVolData)
		return (m_pVolData->clear(data));
	else
		return (false);
}


vlDim vlVolume::dim() const
{
	if(m_pVolData)
		return (m_pVolData->dim());
	else
		return (vlDim(0,0,0));
}


vlStep vlVolume::stepping() const
{
	if(m_pVolData)
		return (m_pVolData->stepping());
	else
		return (vlStep(0,0,0));
}


vlUnit vlVolume::units() const
{
	if(m_pVolData)
		return (m_pVolData->units());
	else
		return (vlUnit(0,0,0));
}


uint16 vlVolume::bitsPerVoxel() const
{
	if(m_pVolData)
		return (m_pVolData->bitsPerVoxel());
	else
		return (0);
}


uint16 vlVolume::bytesPerVoxel() const
{
	if(m_pVolData)
		return (m_pVolData->bytesPerVoxel());
	else
		return (0);
}


uint64 vlVolume::voxelCount() const
{
	if(m_pVolData)
		return (m_pVolData->voxelCount());
	else
		return (0);
}


vlDataType vlVolume::dataType() const
{
	if(m_pVolData)
		return (m_pVolData->dataType());
	else
		return (UnknownDataType);
}

/// Returns the type info of the stored data
const char * vlVolume::typeInfo() const
{
	if(m_pVolData)
		return (m_pVolData->typeInfo());
	else
		return (0L);
}


vlLayoutType vlVolume::dataLayout() const
{
	if(m_pVolData)
		return (m_pVolData->layout());
	else
		return (vlLayout::Linear);
};


bool vlVolume::isValid() const
{
	if(m_pVolData) {
		if(!m_pVolData->isValid())
			return (false);
		else
			return (true);
	} else
		return (false);
}

/// Returns the state of the dirty flag
bool vlVolume::isDirty() const
{
	if(m_pVolData)
		return m_pVolData->isDirty();
	else
		return (false);
}

/// Set/Reset the dirty flag
void vlVolume::setDirty(bool dirty)
{
	if(m_pVolData)
		m_pVolData->setDirty(dirty);
}

void * vlVolume::getVoxelVoidPtr(const vlPoint3ui & position) const
{
	if(m_pVolData)
		return (m_pVolData->getVoxelVoidPtr(position));
	else
		return 0L;
}

void vlVolume::getPosition(vlPoint3f  *p)
{
  *p=position;
}
void vlVolume::setPosition(vlPoint3f p)
{
 position=p;
}

vlPoint3f vlVolume::getCenter()
{
  float x,y,z;
  x=(float) ( (float)(dim().x())/2.0   -0.5);
  y=(float) ( (float)(dim().y())/2.0   -0.5);
  z=(float) ( (float)(dim().z())/2.0   -0.5);
  return vlPoint3f(x*units().x()+position.x(),y*units().y()+position.y(),z*units().z()+position.z());
}
