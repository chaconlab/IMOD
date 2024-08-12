/***************************************************************************
                          frotate.cpp  -  description
                             -------------------
    begin                : Wed Dec 1 2004
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



#include <fftw3.h>
#include "floatops.h"
#include "vlmrot.h"
#include "vldefine.h"


 static fftwf_complex MRI_cxb;
 /*#define CMULT(u,v) ( MRI_cxb[0] = u[0] * v[0] - u[1] * v[1] , \
                    MRI_cxb[1] = u[0] * v[1] + u[1] * v[0] , MRI_cxb )
*/
#define CMULT2(u,v,r) ( MRI_cxb[0] = u[0] * v[0] - u[1] * v[1] , \
                    MRI_cxb[1] = u[0] * v[1] + u[1] * v[0] , \
                    r[0]=MRI_cxb[0],r[1]=MRI_cxb[1] )



 /*****************fumciones auxiliares de rotiacion**********/
/*Funcion auxiliar*/
  float normShear( Tshear sh )
  {
   float top=0.0 , val ;
   int ii , jj ;

   for( ii=0 ; ii < 3 ; ii++ ){
     jj  = sh.ax[ii] ;
     val = fabs( sh.scl[ii][(jj+1)%3] ) ; if( val > top ) top = val ;
     val = fabs( sh.scl[ii][(jj+2)%3] ) ; if( val > top ) top = val ;
   }

   return top ;
  }

  mrot * permutarMatriz(mrot p , float *shift,int ox1,int ox2,int ox3)
  {
   int ii , jj , pi[3] ;
   mrot *m= new mrot();;
   float shift2[3];

   pi[0] = ox1 ; pi[1] = ox2 ; pi[2] = ox3 ;


   for( ii=0 ; ii < 3 ; ii++ )
   {
      for( jj=0 ; jj < 3 ; jj++ )
         m->putP(ii,jj,p.getP(pi[ii],pi[jj]));

      shift2[ii] = shift[ pi[ii] ] ;
   }
   for( ii=0 ; ii < 3 ; ii++ )
      shift[ii]=shift2[ii];


   return m ;
  }

  Tshear permutarShear(Tshear shin,int ox1, int ox2, int ox3)
  {
    Tshear shout ;
    int ii , ain,aout , pi[3] ;


    pi[0] = ox1 ; pi[1] = ox2 ; pi[2] = ox3 ;

    for( ii=0 ; ii < 4 ; ii++ ){

      ain  = shin.ax[ii] ;   /* axis of input */
      aout = pi[ain] ;       /* axis of output */

      shout.ax[ii]         = aout ;             /* store new axis */
      shout.scl[ii][pi[0]] = shin.scl[ii][0] ;  /* permuted scalings */
      shout.scl[ii][pi[1]] = shin.scl[ii][1] ;
      shout.scl[ii][pi[2]] = shin.scl[ii][2] ;
      shout.sft[ii]        = shin.sft[ii] ;     /* copy shift */
    }

    return shout ;
  }

  /*Funcion auxiliar*/
  Tshear apply_shear( mrot m , float *shift,int ox1,int ox2,int ox3 )
  {
   mrot *mm ;
   Tshear out;

   /* permute the input matrix and vector to the desired order */

//   std::cout<<"apply_shear>Permutacion de matriz"<<std::endl;

   mm=permutarMatriz(m,shift,ox1,ox2,ox3);

   /* compute the Sx Sz Sy Sx factorization */

   //std::cout<<"apply_shear>Descomposicion de matriz"<<std::endl;
   if(!mm->decomposition(shift,&out))
   {
      out.ax[0]=-10;
      out.ax[1]=-10;
      out.ax[2]=-10;
      out.ax[3]=-10;
      return out;
   }

   /* permute the shear factorization back */

   //std::cout<<"apply_shear>Permutacion de shear "<<out.ax[1]<<std::endl;
   out=permutarShear(out,ox1,ox2,ox3);

   return out ;
  }

  Tshear bestShear(mrot rot,float *shift)
  {
    int ii,jbest=-1;
    Tshear sol[6];
    float val,best=1.e+38;


    sol[0]=apply_shear( rot ,shift,0,1,2 );
    sol[1]=apply_shear( rot ,shift,0,2,1 );
    sol[2]=apply_shear( rot ,shift,1,0,2 );
    sol[3]=apply_shear( rot ,shift,1,2,0 );
    sol[4]=apply_shear( rot ,shift,2,0,1 );
    sol[5]=apply_shear( rot ,shift,2,1,0 );


    for(ii=0;ii<6;ii++)
    {
      if(sol[ii].ax[0]!=-10)
      {
        val=normShear(sol[ii]);

        if(val<best)
        {
          best=val;
          jbest=ii;
        }
      }
    }
    if(jbest==-1)
    {
      jbest=0;
      sol[0].ax[0]=-10;
    }
    return sol[jbest];
  }



  #define RMAX 3
  #define N35 ((RMAX+1)*(RMAX+1))

  int nextup( int idim )
  {
   static int * tf = NULL , * dn = NULL ;
   int ibase , ii ;

   /*-- build table of allowable powers of 3 and 5 [tf],
        the powers of 2 just less than them        [dn],
        and their ratios tf/dn                     [rat].
        Then sort tf and dn to be increasing in rat.     --*/

   if( tf == NULL ){
      int p , q , tt,ff , i=0 , j ; float * rat ; float r ;

      tf  = (int *)   malloc(sizeof(int)  *N35) ;
      dn  = (int *)   malloc(sizeof(int)  *N35) ;
      rat = (float *) malloc(sizeof(float)*N35) ;

      /* create tables */

      for( p=0,tt=1 ; p <= RMAX ; p++,tt*=3 ){         /* tt = 3^p */
         for( q=0,ff=1 ; q <= RMAX ; q++,ff*=5,i++ ){  /* ff = 5^q */
            tf[i] = tt * ff ;

            j=2; while( j < tf[i] ){ j*=2; } /* j = power of 2 just >= tf */
            dn[i] = j/2 ;                    /* dn = power of 2 just < tf */

            rat[i] = (float)(tf[i]) / (float)(dn[i]) ;
         }
      }

      /* sort on rat (crude, but it works) */

      do{
         for( i=1,p=0 ; i < N35 ; i++ ){
            if( rat[i-1] > rat[i] ){
               r = rat[i-1] ; rat[i-1] = rat[i] ; rat[i] = r ;
               q = tf [i-1] ; tf [i-1] = tf [i] ; tf [i] = q ;
               q = dn [i-1] ; dn [i-1] = dn [i] ; dn [i] = q ;
               p++ ;
            }
         }
      } while( p > 0 ) ;

      free(rat) ;
   }

   /*-- loop on powers of 2 (ibase);
        we can do FFTs of size = tf*ibase/dn (note 1 < tf/dn < 2);
        sinc tf/dn is sorted, we're scanning in increasing sizes  --*/

   ibase = 1 ;
   while(1){
      if( idim <= ibase ) return ibase ;

      for( ii=0 ; ii < N35 ; ii++ )
         if( dn[ii] <= ibase && idim <= tf[ii]*ibase/dn[ii] )
            return tf[ii]*ibase/dn[ii] ;

      ibase *= 2 ;
   }
  }

/*----------------------------------------------------------------------
   return the next legal value that has only 1 power of 3 and/or 5,
   and that is also even
------------------------------------------------------------------------*/

  int nextup_one35( int idim )
  {
   int jj = idim ;
   do{
      jj = nextup(jj) ;
      if( jj%9 == 0 || jj%25 == 0 || jj%2 == 1 ) jj++ ;
      else                                       return jj ;
   } while(1) ;
   return 0 ; /* cannot be reached */
  }

/*----------------------------------------------------------------------
   return the next legal value that is even [13 May 2003 - RWCox]
------------------------------------------------------------------------*/

  int nextup_even( int idim )
  {
   int jj = idim ;
   do{
      jj = nextup(jj) ;
      if( jj%2 == 1 ) jj++ ;
      else            return jj ;
   } while(1) ;
   return 0 ; /* cannot be reached */
  }

/***************FIM FUNCIONES AUXILIARES ROTACION******************/


 namespace FOPS
 {
 //void Xshear(vlVolume *vol, float a, float b);

 // **********************FUNCIONES DE ROTACION********************************************
  /**
   *Aplicacion de un shear dependiente de a y b en el eje x al volumen vol
   */
  bool xshear(vlVolume *vol, float a, float b, float s)
  {
    float *fj1,*fj0;

    int nx=vol->dim().x();
    int ny=vol->dim().y();
    int nz=vol->dim().z();

    int ny1=ny-1,nz1=nz-1;

    float ny2=0.5*ny1, nz2=0.5*nz1;
    int ii,jj,kk,nup,nst;
    float a0,a1,st;


    //Calculo de nup
    st= fabs(a)*ny2+fabs(b)*nz2+fabs(s);
    if(st<1.e-3) {/*std::cout<<"No shear"<<std::endl;*/return false;}
    nst= nx + 0.5*nx;
    nup=nextup_one35(nst);

    fj0= (float*) malloc (sizeof(float) * nx);
    fj1= (float*) malloc (sizeof(float) * nx);


    for(kk=0;kk<nz;kk++)
    {
      for(jj=0;jj<ny1;jj+=2)
      {
        //extraccion de dos columnas de datos en x del volumen
        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel(vlPoint3ui(ii,jj,kk),fj0[ii]);
          vol->getVoxel(vlPoint3ui(ii,jj+1,kk),fj1[ii]);
        }
        a0= a*(jj-ny2)+b*(kk-nz2)+s;
        a1= a0 +a;

        //Calculo para las dos columnas de datos de los valores de interpolacion
        shifter(nx,nup,a0,fj0,a1,fj1);

        //introducion de dos columnas de datos en x del volumen
        for(ii=0; ii<nx;ii++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),fj0[ii]);
          vol->setVoxel(vlPoint3ui(ii,jj+1,kk),fj1[ii]);
        }
      }

      //Idem a lo anterior para el caso de que quede una sola columna de datos
      if(jj==ny1)
      {
        for(ii=0;ii<nx;ii++)
          vol->getVoxel(vlPoint3ui(ii,jj,kk),fj0[ii]);
        a0=a*(jj-ny2) + b*(kk-nz2) +s;

        shifter(nx,nup,a0,fj0,a1,NULL);


        for(ii=0; ii<nx; ii++)
          vol->setVoxel(vlPoint3ui(ii,jj,kk),fj0[ii]);

      }
    }

    free(fj0);
    free(fj1);
    return true;
  }

  /**
   *Aplicacion de un shear dependiente de a y b en el eje y al volumen vol
   */
  bool yshear(vlVolume *vol, float a, float b, float s)
  {
    float *fj1,*fj0;

    int nx=vol->dim().x();
    int ny=vol->dim().y();
    int nz=vol->dim().z();

    int nx1=nx-1,nz1=nz-1;

    float nx2=0.5 * nx1;
    float nz2=0.5*nz1;
    int ii,jj,kk,nup,nst;
    float a0,a1,st;


    //Calculo de nup
    st= fabs(a)*nx2+fabs(b)*nz2+fabs(s);
    if(st<1.e-3) {/*std::cout<<"No shear"<<std::endl;*/return false;}
    nst= ny + 0.5*ny;
    nup=nextup_one35(nst);

    fj0= (float*) malloc (sizeof(float) * ny);
    fj1= (float*) malloc (sizeof(float) * ny);


    for(kk=0;kk<nz;kk++)
    {
      for(ii=0;ii<nx1;ii+=2)
      {
        //extracion de dos columnas de datos en y del volumen
        for(jj=0;jj<ny;jj++)
        {
          vol->getVoxel(vlPoint3ui(ii,jj,kk),fj0[jj]);
          vol->getVoxel(vlPoint3ui(ii+1,jj,kk),fj1[jj]);
        }
        a0= a*(ii-nx2)+b*(kk-nz2)+s;
        a1= a0 +a;

        //Calculo para las dos columnas de datos de los valores de interpolacion
        shifter(ny,nup,a0,fj0,a1,fj1);

        //introducion de dos columnas de datos en y del volumen
        for(jj=0; jj<ny;jj++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),fj0[jj]);
          vol->setVoxel(vlPoint3ui(ii+1,jj,kk),fj1[jj]);
        }
      }

      //Idem a lo anterior para el caso de que quede una sola columna de datos
      if(ii==nx1)
      {
        for(jj=0;jj<ny;jj++)
          vol->getVoxel(vlPoint3ui(ii,jj,kk),fj0[jj]);
        a0=a*(ii-nx2) + b*(kk-nz2) +s;

        shifter(ny,nup,a0,fj0,a1,NULL);


        for(jj=0; jj<ny; jj++)
          vol->setVoxel(vlPoint3ui(ii,jj,kk),fj0[jj]);

      }
    }

    free(fj0);
    free(fj1);
    return true;
  }


  /**
   *Aplicacion de un shear dependiente de a y b en el eje z al volumen vol
   */
 bool zshear(vlVolume *vol, float a, float b, float s)
  {
    float *fj1,*fj0;

    int nx=vol->dim().x();
    int ny=vol->dim().y();
    int nz=vol->dim().z();

    int nx1=nx-1,ny1=ny-1;

    float nx2=0.5 * nx1, ny2=0.5*ny1;
    int ii,jj,kk,nup,nst;
    float a0,a1,st;


    //Calculo de nup
    st= fabs(a)*nx2+fabs(b)*ny2+fabs(s);
    if(st<1.e-3) {/*std::cout<<"No shear"<<std::endl;*/return false;}
    nst= nz + 0.5*nz;
    nup=nextup_one35(nst);


    fj0= (float*) malloc (sizeof(float) * nz);
    fj1= (float*) malloc (sizeof(float) * nz);


    for(jj=0;jj<ny;jj++)
    {
      for(ii=0;ii<nx1;ii+=2)
      {
        //extracion de dos columnas de datos en z del volumen
        for(kk=0;kk<nz;kk++)
        {
          vol->getVoxel(vlPoint3ui(ii,jj,kk),fj0[kk]);
          vol->getVoxel(vlPoint3ui(ii+1,jj,kk),fj1[kk]);
        }
        a0= a*(ii-nx2)+b*(jj-ny2)+s;
        a1= a0 +a;

        //Calculo para las dos columnas de datos de los valores de interpolacion
        shifter(nz,nup,a0,fj0,a1,fj1);

        //introducion de dos columnas de datos en z del volumen
        for(kk=0; kk<nz;kk++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),fj0[kk]);
          vol->setVoxel(vlPoint3ui(ii+1,jj,kk),fj1[kk]);
        }
      }

      //Idem a lo anterior para el caso de que quede una sola columna de datos
      if(ii==nx1)
      {
        for(kk=0;kk<nz;kk++)
          vol->getVoxel(vlPoint3ui(ii,jj,kk),fj0[kk]);
        a0=a*(ii-nx2) + b*(jj-ny2) +s;

        shifter(nz,nup,a0,fj0,a1,NULL);


        for(kk=0; kk<nz; kk++)
          vol->setVoxel(vlPoint3ui(ii,jj,kk),fj0[kk]);

      }
    }

    free(fj0);
    free(fj1);
    return true;
  }

  /**
  * Shear auxiliary variable. Fourier transform plan
  * fpaln,fplan2: planes para realizacion de la transformada y la transformada inversa de Fourier
  * row,cf,cg: buffer para operar en el dominio de Fourier
  */
  fftwf_plan fplan;
  /**
  * Shear auxiliary variable. Inverse Fourier transform plan
  */
  fftwf_plan fplan2;
  /**
   * Shear Auxiliary variable
   */
  int oldnup=0;
  /**
   * Shear auxiliary variable for operations on the Foutrier domain
   */
  fftwf_complex *row=NULL;
  /**
   * Shear auxiliary variable for operations on the Foutrier domain
   */
  fftwf_complex *cf=NULL;
  /**
   * Shear auxiliary variable for operations on the Foutrier domain
   */
  fftwf_complex *cg=NULL;


  /**
   *Funcion que interpola los valores de dos columnas de datos en un rotacion
   */
  bool shifter(int n,int nup,float af,float *f,float ag,float *g)
  {
    //int nupold=0,nuptop=0;
    int ii, nby2=nup/2,n21=nby2+1;
    fftwf_complex fac,gac;
    fftwf_complex csf,csg;

    double sf,sg,dk;

    //Si el shift es muy grande, no hacer nada
    if ( (af<-n || af>n) && (ag<-n || ag>n) )
      return false;


    //Si el valor de nup ha cambiado, es necesario recalcular los planes de transformadas
    // y reservar nuevo espacio de memoria
    if(nup!=oldnup)
    {
      free(row);
      free(cf);
      free(cg);
      row= (fftwf_complex*) malloc(sizeof(fftwf_complex)*nup);
      cf= (fftwf_complex*) malloc(sizeof(fftwf_complex)*n21);
      cg = (fftwf_complex*) malloc(sizeof(fftwf_complex)*n21);

      if(oldnup!=0)
      {
      fftwf_destroy_plan(fplan2);
      fftwf_destroy_plan(fplan);
      }
      fplan=fftwf_plan_dft_1d(nup,row,row,FFTW_FORWARD,FFTW_ESTIMATE);
      fplan2=fftwf_plan_dft_1d(nup,row,row,FFTW_BACKWARD,FFTW_ESTIMATE);
      oldnup=nup;
    }



    //Introduccion de valores de las columenas de datos en buffer de complejos para
    //realizar la transformada
    if(g!=NULL)
      for(ii=0;ii<n;ii++)
      {
        row[ii][0]=f[ii];
        row[ii][1]=g[ii];
      }
    else
      for(ii=0;ii<n;ii++)
      {
        row[ii][0]=f[ii];
        row[ii][1]=0.0;
      }

    for(ii=n;ii<nup;ii++)
    {
        row[ii][0]=0.0;
        row[ii][1]=0.0;
    }

    //llamada a la transformada de fourier
    fftwf_execute(fplan);


    //Extraccion de valores en el dominio de fourier
    cf[0][0]=2.0*row[0][0];
    cg[0][0]=2.0*row[0][1];
    cf[0][1]=cg[0][1]=0.0;

    for(ii=1;ii< nby2;ii++)
    {
      cf[ii][0]=row[ii][0] + row[nup-ii][0];
      cf[ii][1]=row[ii][1] - row[nup-ii][1];
      cg[ii][0]=row[ii][1] + row[nup-ii][1];
      cg[ii][1]=-row[ii][0] + row[nup-ii][0];
    }

    cf[nby2][0]= 2.0 * row[nby2][0];
    cf[nby2][1]= cg[nby2][1]=0.0;
    cg[nby2][0]= 2.0 * row[nby2][1];


    //Alteracion de valores en el dominio de Fourier para realizar la transformada
    dk= (2.0*M_PI)/nup;
    sf=-af*dk;
    sg=-ag*dk;

    csf[0]=cos(sf);
    csf[1]=sin(sf);
    csg[0]=cos(sg);
    csg[1]=sin(sg);

    fac[0]=gac[0]=1.0;
    fac[1]=gac[1]=0.0;

    for(ii=1;ii<=nby2;ii++)
    {
      //fac=CMULT(csf,fac);
      CMULT2(csf,fac,fac);
      //cf[ii]=CMULT(fac,cf[ii]);
      CMULT2(fac,cf[ii],cf[ii]);

      //gac=CMULT(csg,gac);
      CMULT2(csg,gac,gac);
      //cg[ii]=CMULT(gac,cg[ii]);
      CMULT2(gac,cg[ii],cg[ii]);
    }

    cf[nby2][1]=cg[nby2][1]=0.0;

    row[0][0]=cf[0][0];
    row[0][1]=cg[0][0];
    for(ii=1;ii<nby2;ii++)
    {
      row[ii][0]=cf[ii][0] - cg[ii][1];
      row[ii][1]=cf[ii][1] + cg[ii][0];
      row[nup-ii][0]=cf[ii][0]+cg[ii][1];
      row[nup-ii][1]=-cf[ii][1]+cg[ii][0];
    }

    row[nby2][0]=cf[nby2][0];
    row[nby2][1]=cg[nby2][0];

    //transformada inversa para retornar al dominio real
    fftwf_execute(fplan2);


    //Extraccion de las columnas de datos en los buffers iniciales
    sf=0.5/nup;
    if(g!=NULL)
      for(ii=0;ii<n;ii++)
      {
        f[ii]=sf*row[ii][0];
        g[ii]=sf*row[ii][1];
      }
    else
      for(ii=0;ii<n;ii++)
        f[ii]=sf*row[ii][0];
    return true;
  }


  /**
   * Rotacion y dezplazamiento de un volumen en funcion de los angulos y distancia en voxels
   * indicados como argumentos
   */
  void rotate(vlVolume *vol,float psi, float theta, float phi,float sh1,float sh2,float sh3,int opcion)
  {

    int i;
    int flip0=-1,flip1=-1;
    float shift[3];
    shift[0]=sh1;
    shift[1]=sh2;
    shift[2]=sh3;
    float a,b,s;
    Tshear shear;

    //Calculo de la matriz de rotaciones
    mrot rot(psi,theta,phi,opcion);
    mrot rot2;

    //Si la traza de la matriz de rotaciones es menor que 1
    //quizas es conveniente realizar uno o varios flips sobre el volumen antes de
    //aplicar interpolacion
    if(TRACE(rot)<1.0)
    {
      double top=rot.getP(0,0);
      int itop=0, i1,i2;

      //Calculo de matriz de rotaciones con flip
      if( top<rot.getP(1,1))
      { top=rot.getP(1,1); itop=1; }
      if( top<rot.getP(2,2))
      { top=rot.getP(2,2); itop=2; }
      switch(itop)
      {
        case 0: i1=1; i2=2;rot2.putP(0,0,1);rot2.putP(1,1,-1);rot2.putP(2,2,-1);
        break;
        case 1: i1=0; i2=2;rot2.putP(0,0,-1);rot2.putP(1,1,1);rot2.putP(2,2,-1);
        break;
        case 2: i1=0; i2=1;rot2.putP(0,0,-1);rot2.putP(1,1,-1);rot2.putP(2,2,1);
        break;

      }
      rot2.putP(0,1,0);rot2.putP(0,2,0);rot2.putP(1,0,0);rot2.putP(1,2,0);
      rot2.putP(2,0,0);rot2.putP(2,1,0);

      //si merece la pena, modificar la matriz de rotaciones para realizar flips
      if(rot.getP(i1,i1)+rot.getP(i2,i2)<-0.02)
      {
        rot=rot.mul(rot2);
        flip0=i1;flip1=i2;
      }

    }

    //Buscar la mejor decomposicion de la matriz de rotaciones en shears en los
    //ejes x,y,z para realizar la rotacion
    shear=bestShear(rot,shift);

    //Si no existe descomposicion posible, modificar ligeramente la matriz de rotaciones
    if (shear.ax[0]==-10)
    {

      mrot rotaux1(1.09e-6,0.0,0.0);
      mrot rotaux2=rotaux1.mul( mrot(0.0,1.22e-6,0.0));
      mrot rotaux3=mrot(0.0,0.0,1.37e-6);
      rotaux3=rotaux3.mul(rotaux2);
      rot=rot.mul(rotaux3);
      mrot trans=rotaux3.trans();
      rot=trans.mul(rot);
      shear=bestShear(rot,shift);
    }


    //Si se ha encontrado una descomposicion posible
    if(shear.ax[0]!=-10)
    {
      //Se realiza flip sobre volumen si es necesario
      if(flip0>=0)
      {

        switch(flip0+flip1)
        {
          case 1: flipxy(vol);break;
          case 2: flipxz(vol);break;
          case 3: flipyz(vol);break;
        }
      }

     /*for(i=0;i<4;i++)
      {
    	  if(shear.ax[i]==0)
    	  {
    		  a=shear.scl[i][1];
    	      b=shear.scl[i][2];
    	      s=shear.sft[i];
    	  }
      }
      fprintf(stderr,"%f %f %f\n",a,b,s);
      vlVolume *v=new vlVolume(vol);
      xshear(v,a,b,s);
      FOPS::writeFile(v,"kk.ccp4");
      Xshear(vol,a,b);
      FOPS::writeFile(vol,"kk2.ccp4");
      exit(1);*/

      //Se realizan los diferentes shears sobre el volumen en el orden adecuado
      for(i=0;i<4;i++)
      {
        switch( shear.ax[i]){
        case 0:

          a=shear.scl[i][1];
          b=shear.scl[i][2];
          s=shear.sft[i];
          xshear(vol,a,b,s);

        break;
        case 1:

          a=shear.scl[i][0];
          b=shear.scl[i][2];
          s=shear.sft[i];
          yshear(vol,a,b,s);
        break;
        case 2:

          a=shear.scl[i][0];
          b=shear.scl[i][1];
          s=shear.sft[i];
          zshear(vol,a,b,s);
        break;
        }


      }

      return ;
    }

    //Si no se ha conseguido descomponer la matriz de rotaciones,
    //la rotacion no es posible
    else
    {
    	fprintf(stderr,"No possible decomposition\n");
      // delete vol;
      return ;
    }


  }


  //Realizacion del giro de 180 grados sobre el eje z
  bool flipxy(vlVolume *vol)
  {
    int ii,jj,kk;

    int nx=vol->dim().x(), ny=vol->dim().y(), nz=vol->dim().z();
    int nx1=nx-1;
    int ny1=ny-1, ny2=ny/2;
    float *r1;
    float aux;

    //buffer auxiliar
    r1=(float*)malloc(sizeof(float)*nx);

    for(kk=0;kk<nz;kk++)
    {
      for(jj=0;jj<ny2;jj++)
      {
        //realizacion del intercambio por filas en x
        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel( vlPoint3ui(nx1-ii,ny1-jj,kk),aux);
          r1[ii]=aux;
        }

        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel(vlPoint3ui(nx1-ii,jj,kk),aux);
          vol->setVoxel(vlPoint3ui(ii,ny1-jj,kk),aux);
        }

        for(ii=0;ii<nx;ii++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),r1[ii]);
        }
      }

      if(ny%2==1)
      {

        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel(vlPoint3ui(nx1-ii,jj,kk),aux);
          r1[ii]=aux;
        }

        for(ii=0;ii<nx;ii++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),r1[ii]);
        }
      }
    }

    free(r1);
    return true;
  }

  //Realizacion del giro de 180 grados sobre el eje x
 bool flipyz(vlVolume *vol)
  {
    int ii,jj,kk;

    int nx=vol->dim().x(), ny=vol->dim().y(), nz=vol->dim().z();
    int ny1=ny-1;
    int nz1=nz-1, nz2=nz/2;
    float *r1;
    float aux;

    //buffer auxiliar
    r1=(float*)malloc(sizeof(float)*ny);

    for(ii=0;ii<nx;ii++)
    {
      for(kk=0;kk<nz2;kk++)
      {
        //realizacion del intercambio por filas en y
        for(jj=0;jj<ny;jj++)
        {
          vol->getVoxel( vlPoint3ui(ii,ny1-jj,nz1-kk),aux);
          r1[jj]=aux;
        }
        for(jj=0;jj<ny;jj++)
        {
          vol->getVoxel(vlPoint3ui(ii,ny1-jj,kk),aux);
          vol->setVoxel(vlPoint3ui(ii,jj,nz1-kk),aux);
        }
        for(jj=0;jj<ny;jj++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),r1[jj]);
        }
      }

      if(nz%2==1)
      {

        for(jj=0;jj<ny;jj++)
        {
          vol->getVoxel(vlPoint3ui(ii,ny1-jj,kk),aux);
          r1[jj]=aux;
        }
        for(jj=0;jj<ny;jj++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),r1[jj]);
        }
      }
    }

    free(r1);
    return true;
  }

 //Realizacion del giro de 180 grados sobre el eje y
 bool flipxz(vlVolume *vol)
  {
    int ii,jj,kk;

    int nx=vol->dim().x(), ny=vol->dim().y(), nz=vol->dim().z();
    int nx1=nx-1;
    int nz1=nz-1, nz2=nz/2;
    float *r1;
    float aux;

    //buffer auxiliar
    r1=(float*)malloc(sizeof(float)*nx);

    for(jj=0;jj<ny;jj++)
    {
      for(kk=0;kk<nz2;kk++)
      {
        //realizacion del intercambio por filas en x
        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel( vlPoint3ui(nx1-ii,jj,nz1-kk),aux);
          r1[ii]=aux;
        }
        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel(vlPoint3ui(nx1-ii,jj,kk),aux);
          vol->setVoxel(vlPoint3ui(ii,jj,nz1-kk),aux);
        }
        for(ii=0;ii<nx;ii++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),r1[ii]);
        }
      }

      if(nz%2==1)
      {

        for(ii=0;ii<nx;ii++)
        {
          vol->getVoxel(vlPoint3ui(nx1-ii,jj,kk),aux);
          r1[ii]=aux;
        }
        for(ii=0;ii<nx;ii++)
        {
          vol->setVoxel(vlPoint3ui(ii,jj,kk),r1[ii]);
        }
      }
    }

    free(r1);
    return true;
  }


void new_invert(float a[3][3], float b[3][3]) ;

vlVolume* rotate_simple(vlVolume *vol, float psi, float theta, float phi, int conv)
{
	float sx, cxx;
	float sy, cy;
	float sz, cz;
	int i,j,k;
	float x,y,z;
	//int x,y,z;
	float value=1.0;
	float c[3];

	float a[3][3],b[3][3];
	vlPoint3f auxZ,auxY,auxX;
	vlPoint3f B1,B2,B3;
	vlPoint3f currPos;

	c[0]=floor(vol->dim().x()/2.0);
	c[1]=floor(vol->dim().y()/2.0);
	c[2]=floor(vol->dim().z()/2.0);


//	sx = SinCos(psi, &cxx);
	sx = sin(psi); cxx=cos(psi);
//	sy = SinCos(theta, &cy);
	sy = sin(theta); cy=cos(theta);
//	sz = SinCos(phi, &cz);
	sz = sin(phi); cz=cos(phi);


	switch(conv)
	{
	case 0:
		a[0][0]= cxx*cz - cy*sz*sx;
		a[0][1]= cxx*sz + cy*cz*sx;
		a[0][2]= sx*sy;
		a[1][0]= -sx*cz - cy*sz*cxx;
		a[1][1]= -sx*sz + cy*cz*cxx;
		a[1][2]= cxx*sy;
		a[2][0]= sy*sz;
		a[2][1]= -sy*cz;
		a[2][2]= cy;
		break;
	case 1:
		a[0][0] = cy * cz;
		a[0][1] = cy * sz;
		a[0][2] = -sy;

		a[1][0] = cz * sx * sy - cxx * sz;
		a[1][1] = cxx * cz + sx * sy * sz;
		a[1][2] = cy * sx;

		a[2][0] = cxx * cz * sy + sx * sz;
		a[2][1] = cxx * sy * sz - cz * sx;
		a[2][2] = cxx * cy;
		break;
	case 2:
		a[0][0]= -sx*sz+cy*cz*cxx;
		a[0][1]= sx*cz+cy*sz*cxx;
		a[0][2]= -cxx*sy;
		a[1][0]= -cxx*sz-cy*cz*sx;
		a[1][1]= cxx*cz-cy*sz*sx;
		a[1][2]= sx*sy;
		a[2][0]= sy*cz;
		a[2][1]= sy*sz;
		a[2][2]= cy;
		break;
	default:
		fprintf(stderr,"Error: Bad rotational convention (0,1,2) %d\n",conv);
		exit(1);
		break;
	}




	new_invert(a,b);


	vlVolume *vol2=new vlVolume(vol);
	vlPoint3f pos;
	vol->getPosition(&pos);
	vol2->setPosition(pos);
	vol2->clear(0.0);




	B1.x(b[0][0]);
	B1.y(b[1][0]);
	B1.z(b[2][0]);
	B2.x(b[0][1]);
	B2.y(b[1][1]);
	B2.z(b[2][1]);
	B3.x(b[0][2]);
	B3.y(b[1][2]);
	B3.z(b[2][2]);

	x=-c[0]-1.0;
	y=-c[1]-1.0;
	z=-c[2]-1.0;

	auxZ.x(b[0][0] * x +b[0][1] * y + b[0][2] * z + c[0]);
	auxZ.y(b[1][0] * x +b[1][1] * y + b[1][2] * z + c[1]);
	auxZ.z(b[2][0] * x +b[2][1] * y + b[2][2] * z + c[2]);

	//z=-c[2];
	for(k=0;k<vol2->dim().z();k++,z+=1.0)
	{
		//y   = -c[1];
		auxZ+=B3;
		auxY=auxZ;
		for(j=0;j<vol2->dim().y();j++,y+=1.0)
		{

			auxY+=B2;
			//x   = -c[0];
			currPos=auxY;

			for(i=0;i<vol2->dim().x();i++,x+=1.0)
			{
				currPos+=B1;

			//	fprintf(stderr,"currPos1= %f %f %f\n",currPos.x(),currPos.y(),currPos.z());
				//currPos.x(b[0][0] * x +b[0][1] * y + b[0][2] * z + c[0]);
				//currPos.y(b[1][0] * x +b[1][1] * y + b[1][2] * z + c[1]);
				//currPos.z(b[2][0] * x +b[2][1] * y + b[2][2] * z + c[2]);
			//	fprintf(stderr,"currPos2= %d %d %d\n",i,j,k);
				if(currPos.x()<0 || currPos.x()>=vol2->dim().x() ||
						currPos.y()<0 || currPos.y()>=vol2->dim().y() ||
						currPos.z()<0 || currPos.z()>=vol2->dim().z())
					value=0.0;
				else
					vol->getVoxel(currPos,value);


				if(value!=0.0)
				{
					//fprintf(stderr,"%f\n",value);
					vol2->setVoxel(vlPoint3ui(i,j,k),value);

				}
				//fprintf(stderr,"flag3\n");
			}
		}
	}






	return vol2;
}


void rotate_interp(vlVolume *vol, vlVolume *vol2, float psi, float theta, float phi, int conv)
{
	float sx, cxx;
	float sy, cy;
	float sz, cz;
	int i,j,k;
	float x,y,z;

	float value=1.0;
	float c[3];

	float a[3][3],b[3][3];
	vlPoint3f auxZ,auxY,auxX;
	vlPoint3f B1,B2,B3;
	vlPoint3f currPos;

	c[0]=floor(vol->dim().x()/2.0);
	c[1]=floor(vol->dim().y()/2.0);
	c[2]=floor(vol->dim().z()/2.0);


	//	sx = SinCos(psi, &cxx);
		sx = sin(psi); cxx=cos(psi);
	//	sy = SinCos(theta, &cy);
		sy = sin(theta); cy=cos(theta);
	//	sz = SinCos(phi, &cz);
		sz = sin(phi); cz=cos(phi);

	//fprintf(stderr,"%f  %3.3f %3.3f\n",psi, theta, phi);

	switch(conv)
	{
	case 0:
		a[0][0]= cxx*cz - cy*sz*sx;
		a[0][1]= cxx*sz + cy*cz*sx;
		a[0][2]= sx*sy;
		a[1][0]= -sx*cz - cy*sz*cxx;
		a[1][1]= -sx*sz + cy*cz*cxx;
		a[1][2]= cxx*sy;
		a[2][0]= sy*sz;
		a[2][1]= -sy*cz;
		a[2][2]= cy;
		break;
	case 1:
		a[0][0] = cy * cz;
		a[0][1] = cy * sz;
		a[0][2] = -sy;

		a[1][0] = cz * sx * sy - cxx * sz;
		a[1][1] = cxx * cz + sx * sy * sz;
		a[1][2] = cy * sx;

		a[2][0] = cxx * cz * sy + sx * sz;
		a[2][1] = cxx * sy * sz - cz * sx;
		a[2][2] = cxx * cy;
		break;
	case 2:
		a[0][0]= -sx*sz+cy*cz*cxx;
		a[0][1]= sx*cz+cy*sz*cxx;
		a[0][2]= -cxx*sy;
		a[1][0]= -cxx*sz-cy*cz*sx;
		a[1][1]= cxx*cz-cy*sz*sx;
		a[1][2]= sx*sy;
		a[2][0]= sy*cz;
		a[2][1]= sy*sz;
		a[2][2]= cy;
		break;
	default:
		fprintf(stderr,"Error: Bad rotational convention (0,1,2) %d\n",conv);
		exit(1);
		break;
	}


	new_invert(a,b);



	vol2->clear(0.0);

	B1.x(b[0][0]);
	B1.y(b[1][0]);
	B1.z(b[2][0]);
	B2.x(b[0][1]);
	B2.y(b[1][1]);
	B2.z(b[2][1]);
	B3.x(b[0][2]);
	B3.y(b[1][2]);
	B3.z(b[2][2]);

	x=-c[0]-1.0;
	y=-c[1]-1.0;
	z=-c[2]-1.0;

	auxZ.x(b[0][0] * x +b[0][1] * y + b[0][2] * z + c[0]);
	auxZ.y(b[1][0] * x +b[1][1] * y + b[1][2] * z + c[1]);
	auxZ.z(b[2][0] * x +b[2][1] * y + b[2][2] * z + c[2]);



	//z=-c[2];
	for(unsigned k=0;k<vol2->dim().z();k++,z+=1.0)
	{
		//y   = -c[1];
		auxZ+=B3;
		auxY=auxZ;
		for(unsigned j=0;j<vol2->dim().y();j++,y+=1.0)
		{

			auxY+=B2;
			//x   = -c[0];
			currPos=auxY;

			for(unsigned i=0;i<vol2->dim().x();i++,x+=1.0)
			{
				currPos+=B1;


			if(currPos.x()<0 || currPos.x()>=vol2->dim().x() ||
						currPos.y()<0 || currPos.y()>=vol2->dim().y() ||
						currPos.z()<0 || currPos.z()>=vol2->dim().z())
					value=0.0;
				else
					vol->getVoxel(currPos,value);



					vol2->setVoxel(vlPoint3ui(i,j,k),value);


			}
		}
	}



}




void mult(float a[3][3],float b[3][3],float c[3][3])
{
	int i,j,k;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			c[i][j]=0;
			for(k=0;k<3;k++)
			{
				c[i][j]+=a[k][i]*b[j][k];
			}
		}
}


void new_invert(float a[3][3], float  b[3][3]) {

        float m00=a[0][0];
        float m01=a[0][1];
        float m02=a[0][2];

        float m10=a[1][0];
        float m11=a[1][1];
        float m12=a[1][2];

        float m20=a[2][0];
        float m21=a[2][1];
        float m22=a[2][2];


        float cof00 = m11*m22-m12*m21;
        float cof11 = m22*m00-m20*m02;
        float cof22 = m00*m11-m01*m10;
        float cof01 = m10*m22-m20*m12;
        float cof02 = m10*m21-m20*m11;
        float cof12 = m00*m21-m01*m20;
        float cof10 = m01*m22-m02*m21;
        float cof20 = m01*m12-m02*m11;
        float cof21 = m00*m12-m10*m02;

        float det = m00* cof00 + m02* cof02 -m01*cof01;

        b[0][0] =   cof00/det;
        b[0][1] = - cof10/det;
        b[0][2] =   cof20/det;
        b[1][0] = - cof01/det;
        b[1][1] =   cof11/det;
        b[1][2] = - cof21/det;
        b[2][0] =   cof02/det;
        b[2][1] = - cof12/det;
        b[2][2] =   cof22/det;

}
/*
void Xshear(vlVolume *vol, float a, float b)
{
	float *line;
	int x,y,z;
	float xx;
	float half=(vol->dim().x()/2.0)-2.5;
	float auxY=-a,auxZ=-b;
	float value;
	line=(float*)malloc(sizeof(float)*vol->dim().x());

	for(y=0;y<vol->dim().y();y++)
	{
		auxY+=a;
		for(z=0;z<vol->dim().z();z++)
		{
			if(z==0)
				auxZ=0;
			else
				auxZ+=b;

			for(x=0;x<vol->dim().x();x++)
			{
				line[x]=0.0;
				xx=(float)x-auxY-auxZ-2.0;//+half;

				if(xx>=0 && xx<vol->dim().x())
				{
					vol->getVoxel(vlPoint3f(xx,(float)y,(float)z),value);

					if(value!=0)
						line[x]=value;
				}
			}

			for(x=0;x<vol->dim().x();x++)
			{

					vol->setVoxel(vlPoint3ui(x,y,z),line[x]);
			}

		}
	}
	free(line);
}*/
// **********************FIN FUNCIONES DE ROTACION********************************************


  }
