/***************************************************************************
                          vlneighborhood3d.h  -  description
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



#ifndef _vlKernel_h
#define _vlKernel_h

#include "vlconstants.h"
#include "vlneighborhood.h"
             
/**
 * This class creates a neighborhood which includes voxels in the immediate
 * neighborhood in 3D. The neighborhood can be either 6 (along X,Y and Z direction)
 * or 18 (along X,Y,Z + along minor digonals) or 26 (all the previous ones + along major
 * digonals). This type of neighborhood is often needed in image processing.
 *
 * To create this neighborhood, just call the constructor with the count that you need.
 * The count can only have three possible values : 6, 18 or 26.
 *
 * @see vlNeighborhood
 */
class vlKernel
  : public vlNeighborhood
{
public:
  /**
   * Default constructor. Takes count as argument. Count is the number
   * of neighbors to consider. It can only be 6, 18 or 26.
   */
  vlKernel(uint16 count=6);

  ~vlKernel();
};


#endif // _vlNeighborhood3D_h

