/***************************************************************************
                          mrot.cpp  -  description
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

#include "vlmrot.h"
#include <stdio.h>

#include <math.h>
#include "vldefine.h"





/**
 * Constructor del objeto matriz
 */
mrot::mrot(){
}


/**
 * Constructor del objeto matriz introduciendo los valores de la matriz
 */
mrot::mrot(float a11, float a12, float a13, float a21,float a22, float a23, float a31,float a32, float a33)
{
  m[0][0]=a11; m[0][1]=a12; m[0][2]=a13;
  m[1][0]=a21; m[1][1]=a22; m[1][2]=a23;
  m[2][0]=a31; m[2][1]=a32; m[2][2]=a33;
}

/**
 * Constructor del objeto matriz indicando los angulos de Euler que definen la matriz
 * de rotaciones. Criterio X,Y,Z
 */
/*mrot::mrot(float psi, float theta, float phi){


  //Calculo de la matriz de rotaciones
  float sin_psi = sin(psi);
  float cos_psi = cos(psi);
  float sin_theta= sin(theta);
  float cos_theta= cos(theta);
  float sin_phi=sin(phi);
  float cos_phi=cos(phi);
  m[0][0]= cos_theta*cos_phi;
  m[0][1]= cos_theta*sin_phi;
  m[0][2]= -sin_theta;
  m[1][0]= sin_psi*sin_theta*cos_phi - cos_psi*sin_phi;
  m[1][1]= sin_psi*sin_theta*sin_phi + cos_phi*cos_psi;
  m[1][2]= cos_theta*sin_psi;
  m[2][0]= cos_psi*sin_theta*cos_phi + sin_psi*sin_phi;
  m[2][1]= cos_psi*sin_theta*sin_phi - sin_psi*cos_phi;
  m[2][2]= cos_theta*cos_psi;
}*/
/**
 * Constructor del objeto matriz indicando los angulos de Euler que definen la matriz
 * de rotaciones. opcion=0 => convencion ZXZ. opcion=1 =>convencion= XYZ.
 * opcion=2 convencion ZYZ
 */
mrot::mrot(float psi, float theta, float phi,int opcion){


  //Calculo de la matriz de rotaciones
  float sin_psi = sin(psi);
  float cos_psi = cos(psi);
  float sin_theta= sin(theta);
  float cos_theta= cos(theta);
  float sin_phi=sin(phi);
  float cos_phi=cos(phi);

  if(opcion==0)
  {
    m[0][0]= cos_psi*cos_phi - cos_theta*sin_phi*sin_psi;
    m[0][1]= cos_psi*sin_phi + cos_theta*cos_phi*sin_psi;
    m[0][2]= sin_psi*sin_theta;
    m[1][0]= -sin_psi*cos_phi - cos_theta*sin_phi*cos_psi;
    m[1][1]= -sin_psi*sin_phi + cos_theta*cos_phi*cos_psi;
    m[1][2]= cos_psi*sin_theta;
    m[2][0]= sin_theta*sin_phi;
    m[2][1]= -sin_theta*cos_phi;
    m[2][2]= cos_theta;
    return;
  }
  if(opcion==1)
  {
    m[0][0]= cos_theta*cos_phi;
    m[0][1]= cos_theta*sin_phi;
    m[0][2]= -sin_theta;
    m[1][0]= sin_psi*sin_theta*cos_phi - cos_psi*sin_phi;
    m[1][1]= sin_psi*sin_theta*sin_phi + cos_phi*cos_psi;
    m[1][2]= cos_theta*sin_psi;
    m[2][0]= cos_psi*sin_theta*cos_phi + sin_psi*sin_phi;
    m[2][1]= cos_psi*sin_theta*sin_phi - sin_psi*cos_phi;
    m[2][2]= cos_theta*cos_psi;
    return;
  }
  if(opcion==2)
  {
    m[0][0]= -sin_psi*sin_phi+cos_theta*cos_phi*cos_psi;
    m[0][1]= sin_psi*cos_phi+cos_theta*sin_phi*cos_psi;
    m[0][2]= -cos_psi*sin_theta;
    m[1][0]= -cos_psi*sin_phi-cos_theta*cos_phi*sin_psi;
    m[1][1]= cos_psi*cos_phi-cos_theta*sin_phi*sin_psi;
    m[1][2]= sin_psi*sin_theta;
    m[2][0]= sin_theta*cos_phi;
    m[2][1]= sin_theta*sin_phi;
    m[2][2]= cos_theta;
    return;
  }
  printf("mrot> Error: Bad option creating matrix\n");
}



/**
 * Destructor de la clase
 */
mrot::~mrot()
{}

/**
 * Obtencion del valor de la matriz de rotacion en la posicion indicada
 */
float mrot::getP(int x, int y)
{
  if((x>=0 && x<=2) && (y>=0 && y<=2))
    return m[x][y];
  else
    return WRONG;
}

/**
 * Introduccion del valor en la posicion indicada de la matriz de rotacion
 */
bool mrot::putP(int x, int y, float value)
{
  if((x>=0 && x<=2) && (y>=0 && y<=2))
  {
    m[x][y]=value;
    return true;
  }
  else
    return false;
}


/**
 * Descomposicion de la matriz de rotaciones en un conjunto de Shears siguiendo
 * siempre la permutacion shearX, shearZ, shearY, ShearX
 */
bool mrot::decomposition(float shift[],Tshear *out)
{
   /* input variables */

   double q11,q12,q13,q21,q22,q23,q31,q32,q33 , xdel,ydel,zdel ;

   /* computed parameters */

   double f,bx2,cx2,az,bz,ay,cy,bx1,cx1 , dx,dy,dz ;


   /* internals (created by Maple) */

   double t1, t3, t4, t5, t6, t7, t8, t9, t10, t11,
          t12, t13, t15, t16, t17, t18, t19, t20, t22, t23,
          t24, t25, t26, t27, t28, t29, t30, t32, t34, t35,
          t36, t37, t38, t44, t45, t47, t50, t51, t53, t54,
          t55, t57, t61, t62, t64, t66, t67, t68, t69, t70,
          t73, t75, t77, t78, t79, t80, t81, t84, t85, t86,
          t87, t89, t90, t92, t94, t96, t102, t107, t109, t113,
          t118, t119, t121, t123, t125, t127, t129, t131, t132, t134,
          t141, t145, t148, t150, t151, t157, t160, t163, t164, t167,
          t185, t190, t193, t194, t195, t203, t206, t207, t210, t220,
          t221, t224, t230, t233, t238, t240, t241, t252, t264, t267,
          t269, t275, t292;


   /* load inputs into local variables */

   q11=m[0][0];
   q12=m[0][1];
   q13=m[0][2];
   q21=m[1][0];
   q22=m[1][1];
   q23=m[1][2];
   q31=m[2][0];
   q32=m[2][1];
   q33=m[2][2];
   xdel = (double)shift[0] ; ydel = (double)shift[1] ; zdel = (double)shift[2] ;

   /* the code generated by Maple, slightly massaged */

      ay = q21;
      dy = ydel;
      t1 = q21*q12;
      t3 = q13*q22;
      t4 = t3*q31;
      t5 = q21*q13;
      t6 = t5*q32;
      t7 = q23*q11;
      t8 = t7*q32;
      t9 = q12*q23;
      t10 = t9*q31;
      t11 = q22*q11;
      t12 = t11*q33;
      t13 = t1*q33+t4-t6+t8-t10-t12;
      t15 = q32*q32;
      t16 = t15*q32;
      t17 = q21*q21;
      t18 = t17*q21;
      t19 = t16*t18;
      t20 = q22*q22;
      t22 = q31*q31;
      t23 = t22*q32;
      t24 = q21*t20*t23;
      t25 = t20*q22;
      t26 = t22*q31;
      t27 = t25*t26;
      t28 = t15*t17;
      t29 = q22*q31;
      t30 = t28*t29;
      t32 = t13*t13;

      t34 = (-t19-3.0*t24+t27+3.0*t30)*t32 ;
           if( t34 > 0.0 ) t34 =   pow(  t34 , 0.333333333333333 ) ;
      else if( t34 < 0.0 ) t34 = - pow( -t34 , 0.333333333333333 ) ;
      else                 t34 = 0.0 ;

      if( t13 == 0.0 ) return false ;

      t35 = 1/t13*t34;
      t36 = t35+q31;
      t37 = t36*q21;
      t38 = q12*q33;
      t44 = t36*q23;
      t45 = q11*q32;
      t47 = t36*q12;
      t50 = t36*q22;
      t51 = q11*q33;
      t53 = q32*t17;
      t54 = t53*q12;
      t55 = q32*q21;
      t57 = q32*q31;
      t61 = q32*q23*q11*q31;
      t62 = q31*q21;
      t64 = q22*q12;
      t66 = t22*q23;
      t67 = t66*q12;
      t68 = t22*q13;
      t69 = t68*q22;
      t70 = t29*t51;
      t73 = -t37*t38-t36*q13*t29+t37*q13*q32-t44*t45+t47*q23*q31+t50*t51+t54-
             t55*t11-t57*t5+t61+t62*t38-t62*t64-t67+t69-t70+q31*t20*q11;
      t75 = t20*t22;

      t77 = (t28-2.0*t29*t55+t75) ;
      if( t77 == 0.0 ) return false ;
      t77 = 1/t77 ;

      cx2 = t73*t77;
      t78 = t44*q31;
      t79 = t62*q22;
      t80 = t62*q33;
      t81 = t37*q33;
      t84 = t34*t34;

      if( t84 == 0.0 ) return false ;
      t85 = 1/t84;

      cy = (-t78+t79-t80+t81-t53+t66)*t32*t85;
      t86 = q21*t22;
      t87 = t64*t36;
      t89 = t17*q12;
      t90 = t89*t36;
      t92 = t51*t22;
      t94 = t68*q32;
      t96 = t36*t36;
      t102 = t51*q31;
      t107 = t11*t36;
      t109 = t38*t22;
      t113 = t86*t87-t57*t90+2.0*t50*t92+2.0*t37*t94+t96*t22*t3+t96*q31*t8-3.0*
             t24-t96*q22*t102+3.0*t30+t27-2.0*t36*t26*t3-t19-t62*q32*t107-2.0*t37*t109-t26*
             q21*t64;
      t118 = q32*q22;
      t119 = t118*q11;
      t121 = q11*t36;
      t123 = t26*q13;
      t125 = t26*q23;
      t127 = q33*t26;
      t129 = t96*q12;
      t131 = t96*q21;
      t132 = t38*q31;
      t134 = t22*t22;
      t141 = q31*q13*q32;
      t145 = -q11*t15*t17*q31+t23*t89+t86*t119+t121*t28-t123*t55+t125*t45+t1*
             t127-t129*t66+t131*t132+t134*q13*q22-t9*t134-2.0*t36*t22*t8-t131*t141-t11*t127+
             2.0*t47*t125;

      if( t34 == 0.0 ) return false ;
      t148 = 1/t34;

      if( q21 == 0.0 ) return false ;
      t150 = 1/q21;

      t151 = t148*t77*t150;
      bx2 = (t113+t145)*t13*t151;
      az = -t35;
      f = (-t29+t55)*t13*t148;
      t157 = ydel*q12;
      t160 = zdel*t17;
      t163 = ydel*t22;
      t164 = t163*q21;
      t167 = ydel*t26;
      t185 = xdel*q21;
      t190 = ydel*q11;
      t193 = -ydel*q22*t51*t26-t157*q23*t134-t160*t129*q33+t164*t119+t163*t54-
             t167*q21*q22*q12+ydel*t134*t3+t157*q21*q33*t26-t167*t6-3.0*ydel*q21*t75*q32+
             t167*t8+3.0*ydel*t15*t17*q22*q31+t185*t20*t26-ydel*t16*t18-t190*t28*q31;
      t194 = zdel*q21;
      t195 = t125*q12;
      t203 = xdel*t18;
      t206 = xdel*t17;
      t207 = q22*t22;
      t210 = t160*q32;
      t220 = zdel*t18;
      t221 = q32*q12;
      t224 = t123*q22;
      t230 = t194*t96;
      t233 = t194*t195+t194*t22*t12+t160*t94-t194*q32*t7*t22+t203*q31*t15-2.0*
             t206*t207*q32-t210*t107-t194*t75*q11+t194*q31*t20*q11*t36+t210*t11*q31-t220*
             t221*q31-t194*t224+t160*t207*q12+ydel*t25*t26-t230*t8-t160*t109;
      t238 = t194*t36;
      t240 = ydel*t96;
      t241 = t240*q21;
      t252 = ydel*t36;
      t264 = -t203*t36*t15+t230*t12+2.0*t238*t61-t241*t141-t240*q22*t102+t220*
             t221*t36-t160*q31*t87+2.0*t206*t36*t29*q32-2.0*t252*t22*t8+t164*t87+t190*t36*
             t17*t15-t240*t67-2.0*t252*t224+2.0*t252*q22*t92+t241*t132;
      t267 = t160*t36;
      t269 = ydel*q31;
      t275 = t252*q21;
      t292 = -2.0*t238*t67+t240*t69+2.0*t267*t132-t269*q32*t90-t185*t36*t20*t22
             +2.0*t275*t94+t230*t10+2.0*t238*t69+2.0*t252*t195-t230*t4-2.0*t275*t109-2.0*
             t267*t141-2.0*t238*t70+t240*q31*t8-t269*q21*t118*t121+t160*t96*q13*q32;
      dx = -(t193+t233+t264+t292)*t13*t151;
      bz = t36*t150;
      cx1 = -(t78+t79-t80-t96*q23+t81-t53)*t150*t32*t85;
      dz = (-t252+t194)*t150;
      bx1 = -(-t50+t55)*t150*t13*t148;

    out->ax[3] = 0; out->scl[3][0] = f; out->scl[3][1] = bx2; out->scl[3][2] = cx2; out->sft[3] = dx;
    out->ax[2] = 2; out->scl[2][2] = f; out->scl[2][0] = az ; out->scl[2][1] = bz ; out->sft[2] = dz;
    out->ax[1] = 1; out->scl[1][1] = f; out->scl[1][0] = ay ; out->scl[1][2] = cy ; out->sft[1] = dy;
    out->ax[0] = 0; out->scl[0][0] = 1; out->scl[0][1] = bx1; out->scl[0][2] = cx1; out->sft[0] = 0 ;

    return true;
}

/**
 * Introduccion de los valores de la matriz de rotacion
 */
bool mrot::fill(float a11, float a12, float a13, float a21,float a22, float a23, float a31,float a32, float a33)
{
  m[0][0]=a11; m[0][1]=a12; m[0][2]=a13;
  m[1][0]=a21; m[1][1]=a22; m[1][2]=a23;
  m[2][0]=a31; m[2][1]=a32; m[2][2]=a33;
  return true;
}

/**
 * Multiplicacion de la matriz de rotaciones por otra matriz
 */
mrot mrot::mul(mrot m2)
{
  int i,j,k;
  double total;
  mrot out;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      total=0.0;
      for(k=0;k<3;k++)
        total+= (double)m[i][k]*(double)m2.getP(k,j);

      out.putP(i,j,(float)total);
    }
  return out;
}

/**
 * Determinante de la matriz
 */
float mrot::det()
{
 return(m[0][0]*m[1][1]*m[2][2] -m[0][0]*m[1][2]*m[2][1] -m[1][0]*m[0][1]*m[2][2] \
 + m[1][0]*m[0][2]*m[2][1]+ m[2][0]*m[0][1]*m[1][2] -m[2][0]*m[0][2]*m[1][1]);
}

/**
 *Devolucion de la matriz de rotaciones transpuesta
 */
mrot mrot::trans()
{
  mrot out;

   out.m[0][0] =  m[0][0];
   out.m[1][0] =  m[0][1];
   out.m[2][0] =  m[0][2];
   out.m[0][1] =  m[1][0];
   out.m[1][1] =  m[1][1];
   out.m[2][1] =  m[1][2];
   out.m[0][2] =  m[2][0];
   out.m[1][2] =  m[2][1];
   out.m[2][2] =  m[2][2];
   return out;
}

/**
 *Devolucion de la traza de la matriz
 */
float mrot::trace()
{
  return(m[0][0]+m[1][1]+m[2][2]);
}

/**
* Distancia angular entre 2 matrices
*/
float mrot::distance(mrot m2)
{
   int i,j;
   mrot m2T,mult;
   float trace;

  m2T=m2.trans();
  //Atencion: la multiplicacion parece introducir errores
  mult=mul(m2T);
  trace=0.5*(mult.trace()-1.0);


  if(trace>1.0)
    trace=1.0;
  if(trace<-1.0)
    trace=-1.0;

  trace=acos(trace);

  return(trace*180.0/M_PI);
}
