/***************************************************************************
                          vlneighborhood3d.cpp  -  description
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


#include "vlkernel.h"


vlKernel::vlKernel(uint16 count)
    : vlNeighborhood(3)
{
  int max,min,i,j,k;

  if(count != 27 && count != 125) {
	  fprintf(stderr,"vlKernel::ERROR : count for vlKernel has to be 27 or 125.\n");
  } else {
    if(count==27)
    {
      min=-1;
      max=1;
    }
    if(count==125)
    {
      min=-2;
      max=2;
    }

    for(i=min;i<=max;i++)
      for(j=min;j<=max;j++)
        for(k=min;k<=max;k++)
        {
         m_offsets.push_back(vlOffset(i,j,k)); 
        }
        
  }
}

vlKernel::~vlKernel()
{

};

