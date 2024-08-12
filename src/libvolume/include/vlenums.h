/***************************************************************************
                          vlenums.h  -  description
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


#ifndef _vlEnums_h
#define _vlEnums_h


#include <vector>
#include "vlconstants.h"

/**
 * Data type of a possible volume
 */
typedef enum
{
SignedInt8,
UnsignedInt8,
SignedInt16,
UnsignedInt16,
SignedInt32,
UnsignedInt32,
Float,
Double,
TripleUInt8,
UnknownDataType,
Complexf 
/*
  VL_UNSIGNED_BYTE,
  VL_BYTE,
  VL_UNSIGNED_SHORT,
  VL_SHORT,
  VL_UNSIGNED_INT,
  VL_INT,
  VL_FLOAT,
  VL_DOUBLE,
  VL_COLOR_FLOAT,
  VL_COLOR_24BITRGB,
  VL_OTHER
*/
} vlDataType;

typedef uint16 vlLayoutType;
/**
 * The various data layouts supported by OpenVL.
 * The range 0 -- 9 is reserved for internal use.
 * 10 -- 99 will be used for layouts that will be built into OpenVL.
 * 100 -- 999 and above will be used for dynamic-plugin based layouts
 * which are known to OpenVL.
 * 1000-1009 is reserved for dynamic-plugin based layouts which are unknown to OpenVL.
 *
 * If you want to write your own layout plugin (dynamic - i.e., as a shared library),
 * then to start off you can choose one of
 * the 10 layout numbers named "Experimental 3rd party". Choose any number and use it
 * for your plugin. As long as no other plugin uses that same layout type number, you
 * are safe. Once you have a working layout and would like to get a layout type number
 * for your layout, send email to the OpenVL mailing list and you will get a reserved
 * layout number. This layout type will be reserved for your layout only and will not
 * clash with any other layout.
 */
class vlLayout {
	public:
	/// This will use virtual call to call appropriate layout iterator
	static const vlLayoutType VirtualCall = 1;

	/// This is for generic iterators which can work on any data layout - will be slow
	static const vlLayoutType Generic	= 2;

	/// Linear layout
	static const vlLayoutType Linear = 10;

	/// Experimental 3rd party plugin-based layouts - unknown to OpenVL
	static const vlLayoutType Unknown1 = 1000;
	static const vlLayoutType Unknown2 = 1001;
	static const vlLayoutType Unknown3 = 1002;
	static const vlLayoutType Unknown4 = 1003;
	static const vlLayoutType Unknown5 = 1004;
	static const vlLayoutType Unknown6 = 1005;
	static const vlLayoutType Unknown7 = 1006;
	static const vlLayoutType Unknown8 = 1007;
	static const vlLayoutType Unknown9 = 1008;
	static const vlLayoutType Unknown10 = 1009;

	/**
   * This is a special default value - indicates which type of iterator to use
   * when no layout is specified. i.e., when a call like vlVolIter<datatype> is made.
   */
	static const vlLayoutType DefaultLayout = VirtualCall;

  /// Returns the layouts supported by OpenVL
  static std::vector<vlLayoutType> supported() { return m_supportedLayouts; };

protected:
  /**
   * Stores the supported layouts. Add layouts here as and when they are supported.
   * A supported layout means it is complete and the iterator header file is accessible
   * and known to OpenVL. Before adding a new layout here, make sure an entry for that
   * layout is added to the callFuncOnLayout* macros in vlmacros.h.
   * Note1 : Generic and VirtualCall are not real layouts and hence not added.
   * Note2 : The actual definition is in vlenums.cpp file.
   */
  static const std::vector<vlLayoutType> m_supportedLayouts;
};


  typedef uint16 vlInterpolationType;

  /**
  * Generic interpolation types
  */
  class vlInterpolation {
   public:
   // "0" is reserved for unknown. Do not use it.
   static const vlInterpolationType NearestNeighbor = 1;
   static const vlInterpolationType TriLinear = 2;
  };

  typedef uint16 vlVoxelOpType;

  /**
  * Generic Operational types
  */
  class vlVoxelOp {
  public:
  // "0" is reserved for unknown. Do not use it.
  static const vlVoxelOpType Unknown = 0;
  static const vlVoxelOpType CentralDiff = 10;
  static const vlVoxelOpType Sobel = 11;
  static const vlVoxelOpType Fft = 12;
  static const vlVoxelOpType Anonymous = 10000;
  static const vlVoxelOpType AnonymousMax = 20000;
 };
#endif
