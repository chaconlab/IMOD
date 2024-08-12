/***************************************************************************
                          vlmacros.h  -  description
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
  #ifndef vlMacros_h_
#define vlMacros_h_

#include <complex>
#include "vlenums.h"
#include "vltriple.h"
using namespace std;

/**
 * Macro to make the task of calling a templatized function for a given
 * data type, damn easy. This macro should always be used instead of
 * explicitly doing the switch statement. Any new datatype will be reflected
 * in the macro and the user will not have to change a single line of code!
 * NOTE : If you do not want to pass any arguments to the called function,
 * use the macro callFunctionOnDataTypeNoArgs() instead of this one.
 *
 * Example :  A call,
 *
 * callFunctionOnDataType(dataType, retValue, resizeDataT, newDim, dataType, layout);
 *
 * will expand to:
 *
 * switch(dataType) {
 *   case UnsignedInt8: {
 *     unsigned char dummy;
 *     retValue = resizeDataT(dummy, newDim, dataType, layout);
 *   }
 *   break;
 *   case SignedInt8: {
 *     char dummy;
 *     retValue = resizeDataT(dummy, newDim, dataType, layout);
 *   }
 *   break;
 *   .
 *   .  [All possible Data Types]
 *   .
 *   case Double: {
 *     double dummy;
 *     retValue = resizeDataT(dummy, newDim, dataType, layout);
 *   }
 *   break;
 *   default:
 *     cout << "Unsupported datatype!" << endl;
 *   break;
 * }
 *
 *
 * @param dataType the data type variable to look at
 * @param ret the variable into which to store the return value
 * @param func the function name which is to be called
 * @param args the arguments to the function (other than the dummy variable
 *             that is passed as the first argument)
 */


#define callFunctionOnDataType(dataType, ret, func, dim, units) \
        switch(dataType) { \
    case UnsignedInt8: { unsigned char dummy; ret = func(dummy, dim, units); } break; \
    case SignedInt8: { char dummy; ret = func(dummy, dim, units); } break;\
    case UnsignedInt16: { unsigned short dummy; ret = func(dummy, dim, units); } break;\
    case SignedInt16: { short dummy; ret = func(dummy, dim, units); } break; \
    case UnsignedInt32: { unsigned int dummy; ret = func(dummy, dim, units); } break;\
    case SignedInt32: { int dummy; ret = func(dummy, dim, units); } break;\
    case Float: { float dummy; ret = func(dummy, dim, units); } break;\
    case Double: { double dummy; ret = func(dummy, dim, units); } break; \
    case TripleUInt8: { vlTriple<uint8> dummy; ret = func(dummy, dim, units); } break; \
        case Complexf: { complex<float> dummy; ret = func(dummy, dim, units); } break; \
    default: fprintf(stderr,"Unsupported datatype!\n"); break;\
  }





/**
 * This is a utility macro to return the vlDataType for a given variable. For an
 * unknown datatype UnknownDataType is returned.
 */
#define getVariableDataType(variable, datatype) \
  if(typeid(variable) == typeid(uint8)) datatype = UnsignedInt8; \
  else if (typeid(variable) == typeid(int8)) datatype = SignedInt8; \
  else if (typeid(variable) == typeid(uint16)) datatype = UnsignedInt16; \
  else if (typeid(variable) == typeid(int16)) datatype = SignedInt16; \
  else if (typeid(variable) == typeid(uint32)) datatype = UnsignedInt32; \
  else if (typeid(variable) == typeid(int32)) datatype = SignedInt32; \
  else if (typeid(variable) == typeid(float)) datatype = Float; \
  else if (typeid(variable) == typeid(double)) datatype = Double; \
  else if (typeid(variable) == typeid(vlTriple<uint8>)) datatype = TripleUInt8; \
  else datatype = UnknownDataType;

#endif // vlMacros_h_
