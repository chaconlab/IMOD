/***************************************************************************
                          vlipolator_nneighbor.h  -  description
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


 #ifndef _vlInterpolatorNearestNeighbor_h
#define _vlInterpolatorNearestNeighbor_h

#include <stdio.h>

#include "vloffset.h"
#include "vlinterpolator.h"

/**
 * Instance of vlInterpolator with nearest neighborhood operator
 *
 * @author Sarang Lakare <sarang#users.sourceforge.net>
 * @see vlInterpolator
 */
template <typename DataType, vlLayoutType Layout>
class vlInterpolatorNearestNeighbor : public vlInterpolator<DataType, Layout>
{
public:
  vlInterpolationType type() { return vlInterpolation::NearestNeighbor; };

  std::string name() { return ("NearestNeighbor"); };

  /// Returns the layout for which this interpolator is implemented
  vlLayoutType layout() { return Layout; };

  DataType getValueAt(vlVolIterConst<DataType, Layout> & iter, const vlPoint3f & position, bool check=true)
  {
    vlPoint3ui oldPos = iter.pos();
    vlPoint3ui pos((int)(position.x()+0.5), (int)(position.y()+0.5), (int)(position.z()+0.5));
    iter.moveTo(pos);
    DataType value = iter.get();
    iter.moveTo(oldPos);
    return (value);
  };

  /// Returns the value at offset with respect to the iterator. If check is true, check for errors/bounds etc.
  DataType getValueAtOffset(vlVolIterConst<DataType, Layout> & iter, const vlPoint3f & offset, bool check=true)
  {
    vlOffset relativePos;
    if(offset.x() >= 0)
      relativePos.x((int)(offset.x()+0.5));
    else
      relativePos.x((int)(offset.x()-0.5));

    if(offset.y() >= 0)
      relativePos.y((int)(offset.y()+0.5));
    else
      relativePos.y((int)(offset.y()-0.5));

    if(offset.z() >= 0)
      relativePos.z((int)(offset.z()+0.5));
    else
      relativePos.z((int)(offset.z()-0.5));

    return (iter.getRelative(relativePos));
  };
};

#endif // _vlInterpolatorNearestNeighbor_h
