/***************************************************************************
                          vlneighborhood.h  -  description
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



 #ifndef _vlNeighborhood_h
#define _vlNeighborhood_h

#include <stdio.h>

#include <vector>

#include "vlconstants.h"
#include "vloffset.h"

/**
 * A class to define the neighborhood of a voxel. This is a collection of
 * offsets which defines the neighborhood of a voxel.
 */
class vlNeighborhood {

public:
  vlNeighborhood(uint16 dimension=3)
    : m_dimension(dimension)
  {

  };

  virtual ~vlNeighborhood() { };

  /// Returns the dimension of the neighborhood
  uint16 dimension() const { return (m_dimension); };

  /// Get the offsets
  const std::vector<vlOffset> & getOffsets() const { return m_offsets; };

  /// Add an offset to the neighborhood
  bool add(const vlOffset & offset)
  {
    if(offset.dim() != m_dimension) {
    	fprintf(stderr,"vlNeighborhood :: ERROR : Setting offset of wrong dimension (should be %d but is %d\n",m_dimension,offset.dim());
      return (false);
    }
    m_offsets.push_back(offset);
    return (true);
  }

  /// Remove all offsets from the neighborhood
  bool clear()
  {
    m_offsets.clear();
    return (true);
  }

  /// overriding << to enable writing neighborhood to a stream
  friend std::ostream & operator << (std::ostream & os, vlNeighborhood const * const n) {
    return (os << *n);
  }

  /// overriding << to enable writing neighborhood to a stream
  friend std::ostream & operator << (std::ostream & os, vlNeighborhood const & n) {
    std::vector<vlOffset>::const_iterator iter=n.m_offsets.begin(), end=n.m_offsets.end();
    while(iter != end) {
      os << *iter;
      ++iter;
    }
    return os;
  }

protected:
  /// Stores the offsets to the neighborhood voxels from the central voxel
  std::vector<vlOffset> m_offsets;

  /// Stores the dimension of this neighborhood
  uint16 m_dimension;

};




#endif // _vlNeighborhood_h
