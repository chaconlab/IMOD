/***************************************************************************
                          vlconstants.h  -  description
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



/*
 * This file will include constants used across STK
 */

#ifndef _vlConstants_h
#define _vlConstants_h


/*
 * Defining data types for cross-platform compatibility.
 *
 */

//#if defined(_x86_)

  typedef char int8;
  typedef unsigned char uint8;
  typedef short int16;
  typedef unsigned short uint16;
  typedef int int32;
  typedef unsigned int uint32;
  typedef long int int64;
  typedef unsigned long int uint64;


 #define VL_DEFAULT_FILE_EXTENSION "slc"
 #define VL_DEFAULT_FILE_FORMAT "slc"

// add support for other architectures here

//#endif

/// @todo automatic version numbers from those defined in configure.in
// API version of OpenVL
  #define OPENVL_API_MAJOR 0
  #define OPENVL_API_MINOR 3
  #define OPENVL_API_REVISION 0

// Version of OpenVL
  #define OPENVL_VERSION_MAJOR 0
  #define OPENVL_VERSION_MINOR 3
  #define OPENVL_VERSION_REVISION 0

#endif // _vlConstants_h
