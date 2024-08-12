/***************************************************************************
                          vlipolator_trilinear.h  -  description
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


#ifndef _vlInterpolatorTriLinear_h
#define _vlInterpolatorTriLinear_h

#include <stdio.h>

#include "vloffset.h"
#include "vlinterpolator.h"

/*    Trilinear interpolation: x,y,z are the three weights for  */
/*    the x, y, and z directions, and a,b,c,d,e,f,g, and h are  */
/*    the 8 vertices (values) to interpolate between.           */
/*    (x, y, and z are the distances from point a along the x,  */
/*    y, and z axis respectively. )                             */
/*                                                              */
/*                                                              */
/*                         g----------------h                   */
/*            Y           /|               /|                   */
/*                       / |              / |                   */
/*                      /  |             /  |                   */
/*            |        c----------------d   |                   */
/*            |        |   |            |   |                   */
/*            |        |   |            |   |                   */
/*            |        |   |            |   |                   */
/*            |      Z |   |            |   |                   */
/*            |        |   e------------|---f                   */
/*            |     /  |  /             |  /                    */
/*            |    /   | /              | /                     */
/*            |   /    |/               |/                      */
/*            |  /     a----------------b                       */
/*            | /                                               */
/*            |/                                                */
/*            o------------------  X                            */
/*                                                              */
/*                                                              */
/// Trilinear interpolation operation
#define vlTriLinear(x,y,z,a,b,c,d,e,f,g,h)                 \
        ((((a)*(1.0 - (z))               +               \
           (e)*(z))*(1.0 - (y))          +               \
          ((c)*(1.0 - (z))               +               \
           (g)*(z))*(y))*(1.0-(x))       +               \
         (((b)*(1.0 - (z))               +               \
           (f)*(z))*(1.0 - (y))          +               \
          ((d)*(1.0 - (z))               +               \
          (h)*(z))*(y))*(x))


/**
 * Instance of vlInterpolator with trilinear operator
 *
 * @author Sarang Lakare <sarang#users.sourceforge.net>
 * @see vlInterpolator
 */
 template <typename DataType, vlLayoutType Layout>
 class vlInterpolatorTriLinear : public vlInterpolator<DataType, Layout>
 {
 public:
   vlInterpolationType type() { return vlInterpolation::TriLinear; };

   std::string name() { return ("TriLinear"); };

   vlLayoutType layout() { return Layout; };

   DataType getValueAt(vlVolIterConst<DataType, Layout> & iter, const vlPoint3f & position, bool check=true)
   {
     vlPoint3ui oldPos = iter.pos();
     vlPoint3ui pos((int)position.x(), (int)position.y(), (int)position.z());
     iter.moveTo(pos);

     DataType a,b,c,d,e,f,g,h;
     a = iter.get();
     b = iter.getRelativeX(1);
     c = iter.getRelativeY(1);
     d = iter.getRelative(vlOffset(vlTriple<int16>(1,1,0)));
     e = iter.getRelativeZ(1);
     f = iter.getRelative(vlOffset(vlTriple<int16>(1,0,1)));
     g = iter.getRelative(vlOffset(vlTriple<int16>(0,1,1)));
     h = iter.getRelative(vlOffset(vlTriple<int16>(1,1,1)));

 //    std::cout << "abcd.. " << a << "," << b << "," << c << "," << d << ","
 //                           << e << "," << f << "," << g << "," << h << std::endl;
     float x,y,z;
     x = position.x() - (int)position.x();
     y = position.y() - (int)position.y();
     z = position.z() - (int)position.z();

 //    std::cout << "xyz.. " << x << "," << y << "," << z << std::endl;

     DataType value = vlTriLinear(x,y,z,a,b,c,d,e,f,g,h);

     iter.moveTo(oldPos);
     return (value);
   };

   DataType getValueAtOffset(vlVolIterConst<DataType, Layout> & iter, const vlPoint3f & offset, bool check=true)
   {
     int16 ofX = (offset.x() < 0.0)?-1:1;
     int16 ofY = (offset.y() < 0.0)?-1:1;
     int16 ofZ = (offset.z() < 0.0)?-1:1;


     DataType a,b,c,d,e,f,g,h;
     a = iter.get();
     b = iter.getRelativeX(ofX);
     c = iter.getRelativeY(ofY);
     d = iter.getRelative(vlOffset(vlTriple<int16>(ofX,ofY,0)));
     e = iter.getRelativeZ(ofZ);
     f = iter.getRelative(vlOffset(vlTriple<int16>(ofX,0,ofZ)));
     g = iter.getRelative(vlOffset(vlTriple<int16>(0,ofY,ofZ)));
     h = iter.getRelative(vlOffset(vlTriple<int16>(ofX,ofY,ofZ)));

 //    std::cout << "abcd.. " << a << "," << b << "," << c << "," << d << ","
 //                           << e << "," << f << "," << g << "," << h << std::endl;

     DataType temp;
     float x,y,z;
     x = offset.x() - (int)offset.x();
     y = offset.y() - (int)offset.y();
     z = offset.z() - (int)offset.z();

 //    std::cout << "xyz.. " << x << "," << y << "," << z << std::endl;

 //    std::cout << "About to start swapping.. " << std::endl;

     // Perform swaps so that "a" represents the base voxel and offsets are all +ve
     if(ofX < 0) {
 //      std::cout << "Swapping x" << std::endl;
       temp = b;
       b = a;
       a = temp;

       temp = d;
       d = c;
       c = temp;

       temp = f;
       f = e;
       e = temp;

       temp = h;
       h = g;
       g = temp;

       x = x+1.0;
     }

     if(ofY < 0) {
 //      std::cout << "Swapping y" << std::endl;
       temp = c;
       c = a;
       a = temp;

       temp = d;
       d = b;
       b = temp;

       temp = g;
       g = e;
       e = temp;

       temp = h;
       h = f;
       f = temp;

       y = y+1.0;
     }

     if(ofZ < 0) {
 //      std::cout << "Swapping z" << std::endl;
       temp = e;
      e = a;
       a = temp;

       temp = f;
       f = b;
       b = temp;

       temp = g;
       g = c;
       c = temp;

       temp = h;
       h = d;
       d = temp;

       z = z+1.0;
     }

 //    std::cout << "abcd.. " << a << "," << b << "," << c << "," << d << ","
 //                           << e << "," << f << "," << g << "," << h << std::endl;

 //    std::cout << "xyz.. " << x << "," << y << "," << z << std::endl;

     DataType value = vlTriLinear(x,y,z,a,b,c,d,e,f,g,h);

     return (value);
   };
 };

#endif // _vlInterpolatorTriLinear_h
