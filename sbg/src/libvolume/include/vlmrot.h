/***************************************************************************
                          mrot.h  -  description
                             -------------------
    begin                : Tue Sep 7 2004
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

#ifndef MROT_H
#define MROT_H


/**
  *@author Jose Ignacio Garzon
  */


/**
 *Auxiliary structure. It stores shears to rotate and traslate a volume
 */
typedef struct {
   int    ax[4] ;
   double scl[4][3] , sft[4] ;
   int flip0,flip1;
} Tshear ;

#define TRACE(A) ( (A).getP(0,0) + (A).getP(1,1) + (A).getP(2,2) )

/**
 * Rotational matrix class
 */
class mrot {
public:
/**
 * Constructor. Empty matrix
 */
mrot();
/**
 * Constructor. Direct construction
 *
 * @param a11,a12,a13,a21,a22,a23,a31,a32,a33: Values of the rotational matrix
 */
mrot(float a11, float a12, float a13, float a21,float a22, float a23, float a31,float a32, float a33);
/**
 * Constructor. Euler angle construction
 *
 * @param psi,theta,phi: Euler angles
 * @param opcion: Euler angle convention. 0 (ZXZ) 1 (XYZ) 2 (ZYZ).
 */
mrot(float psi, float theta, float phi,int opcion=0);
/**
 * Destructor
 */
~mrot();

/**
 * Introduction of values in the rotational matrix
 *
 * @param a11,a12,a13,a21,a22,a23,a31,a32,a33: Values of the rotational matrix
 */
bool fill(float a11, float a12, float a13, float a21,float a22, float a23, float a31,float a32, float a33);

/**
 * Returns the value of a given position in the rotational matrix
 *
 * @param x,y: matrix coordinates [0,2]
 * @output value of the position
 */
float getP(int x, int y);
/**
 * Inserts the value in a given position in the rotational matrix
 *
 * @param x,y: matrix coordinates [0,2]
 * @param values: value to introduce
 * @output True if position is valid
 */
bool putP(int x, int y, float value);

/**
 * Shear decomposition of the matrix. Permutation followed: shearX, shearZ, shearY, shearX
 *
 * @param shift: x,y,z shifts
 * @param out: output decomposition
 * @output True if the decompositions has been successfully created
 */
bool decomposition(float shift[],Tshear *out);

/**
 * rotational matrix multiplication
 *
 * @param m2: Matrix to multiplicate with
 * @output The new rotational matrix
 */
mrot mul(mrot m2);

/**
 *Returns the transposed?? matrix
 */
mrot trans();
/**
 * Returns the matrix determinant
 */
float det();

/**
 * Returns the matrix trace
 */
float trace();

/**
 * Returns the angular distance between two matrix
 *
 * @param m2: the second matrix
 */
float distance(mrot m2);

private:
	/**
	 * Rotational matrix
	 */
	float m[3][3];
};


#endif
