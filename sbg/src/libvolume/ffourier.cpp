/***************************************************************************
                          ffourier.cpp  -  description
                             -------------------
    begin                : Tue Jun 15 2004
    copyright            : (C) 2004 by
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



#include "fftw3.h"
#include "floatops.h"
#include "vldefine.h"




namespace FOPS
 {

  // ************FILTROS EN EL DOMINIO DE FOURIER***********************************

 void GaussFilter(vlVolume *fourier,float res, int option)
  {

    vlPoint3f radio, sigma;
    float valuex,valuey,valuez;
    complex<float> value;
    int dimfx,dimfy,dimfz;

    // Calculo de la posicion de  las frecuencias limite en cada una de las direcciones
    // en el espacio de Fourier
    // Para cada direccion = mitad de la longitud del volumen partido por la resolucion
    radio.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res*2) );
    radio.y( (float)(fourier->units().y() * fourier->dim().y()) / (res*2) );
    radio.z( (float)(fourier->units().z() * fourier->dim().z()) / (res*2) );

    //Calculo del sigma para cada direccion, depende de la posicion de la frecuencia
    sigma.x(radio.x()/(sqrt(LN2)));
    sigma.y(radio.y()/(sqrt(LN2)));
    sigma.z(radio.z()/(sqrt(LN2)));

//    std::cout <<" Dimensiones del volumen: " << fourier->dim() << std::endl;
//    std::cout <<" Sigmas del volumen: " << sigma << std::endl;
//    std::cout <<" Numero de voxels del volumen: " << fourier->voxelCount() << std::endl;

    // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
    //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
    vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);


//    std::cout <<" centro: " << center << std::endl;

    // Creacion de un iterador para recorrer el espacioo de fourier
    vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

    dimfx=fourier->dim().x();
    dimfy=fourier->dim().y();
    dimfz=fourier->dim().z();

//    int cont=0;
    iterador.begin();
//    std::cout << "Al principio: "<< iterador.get()<< std::endl;
    do
    {
      //Calculo de la posicion de cada voxel teniendo en cuenta que los cuadrantes estan
      //desordenados
      vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

      if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
      if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}

      //Calculo de la funcion de gauss para el voxel y modificacion del valor del
      //voxel en funcion de esta funcion
      valuex= (float)((pos.x()*pos.x())) / ((sigma.x()*sigma.x()+(sigma.y()*sigma.y())+(sigma.z()*sigma.z())) );
      valuey= (float)((pos.y()*pos.y())) / ((sigma.x()*sigma.x()+(sigma.y()*sigma.y())+(sigma.z()*sigma.z())) );
      valuez= (float)((pos.z()*pos.z())) / ((sigma.x()*sigma.x()+(sigma.y()*sigma.y())+(sigma.z()*sigma.z())) );
      value= iterador.get()*complex<float>(expf(-0.5*(valuex+valuey+valuez)),0);
      /*Si el filtro es de paso bajo, el valor calculado es el del voxel*/
      if(option==0)
      {
        iterador.set(value);
      }
      /*Si el filtro es de paso alto, el valor calculado es restado al valor del voxel*/
      else
      {  iterador.set(iterador.get()-value);
      }


    }while(iterador.next());


    /*Fin del filtro*/
  }



/**
* Pablo's
*the constant a for a Gaussian y=exp(-a*x^2) is
*        a=ln(2)/res^2   for res defined at y=0.5
*        a=1/res^2       for res defined at y=1/e
*
*  here we use the res defined in Fourier space
*  eg. Fourier transform of the above function: Y=exp(-pi^2*k^2/a)
*   a=ln(2)*pi^2/res^2   for res defined at Y=0.5
*   a=pi^2/res^2         for res defined at Y=1/e
*                          reciprocal of the 1/2 width of a Gaussian in Fourier space
*
*   Notes:  the gaussian exp 3/2 factor vanishes due to
*                          sigma_3D=sigma_1D*sqrt(3)
*                          sigma1D=(res*1/2)/width
*
*           res_0.5=1/sqrt(ln2)*res_1/e = 1.20112240878*res_1/e
*/
   void GaussFilter(vlVolume *vlm,
                float res)
  {

    vlVolIter<float, vlLayout::Linear> iter(vlm);


    float *vol;
    int dimx,dimy,dimz;
    int dimfx,dimfy,dimfz;
    vlVolDataLayoutBase<float> *dataLayout;

    fftwf_complex *fout;
    fftwf_plan p;
    vlVolume *fourier;

    // get pointer to volume
    dataLayout=iter.getVolume();

    // get pointer to densitiy data
    vol=dataLayout->getVoxelPtr(vlPoint3ui(0,0,0));

    //Dimensiones del array de datos
    dimx=dataLayout->dim().x();
    dimy=dataLayout->dim().y();
    dimz=dataLayout->dim().z();

    // new volume for storing the FT
    fourier= new vlVolume(vlDim((dimx/2)+1,dimy,dimz),Complexf,dataLayout->units());

    // get FT volume pointer
    fout=(fftwf_complex*)((vlVolDataLayoutBase<complex<float> >*)fourier->volData())->getVoxelPtr(vlPoint3ui(0,0,0));

    // Perform FT with FFTW
    p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
                            vol,fout,FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);

    /* fttw fourier constant normalization factor
       <> to FOPS::mul(vlm,1.0/(dimx*dimy*dimz)); after FFT */
    float fttw_cte=log((float)(dimx*dimy*dimz));

    vlPoint3f sigma;
    float g_exp;
    complex<float> value;

    // Calculo de la posicion de  las frecuencias limite en cada una de las direcciones
    // en el espacio de Fourier
    // Para cada direccion = mitad de la longitud del volumen partido por la resolucion
    sigma.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res) );
    sigma.y( (float)(fourier->units().y() * fourier->dim().y()) / (res) );
    sigma.z( (float)(fourier->units().z() * fourier->dim().z()) / (res) );

    //Calculo del sigma para cada direccion, depende de la posicion de la frecuencia
    sigma.x(sigma.x()*sigma.x());
    sigma.y(sigma.y()*sigma.y());
    sigma.z(sigma.z()*sigma.z());

    // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
    //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
    vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);

    // std::cout <<" Dimensiones del volumen: " << fourier->dim() << std::endl;
    // std::cout <<" Sigmas del volumen: " << sigma << std::endl;
    // std::cout <<" Numero de voxels del volumen: " << fourier->voxelCount() << std::endl;
    // std::cout <<" centro: " << center << std::endl;

    // Creacion de un iterador para recorrer todo el espacioo de fourier
    vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

    dimfx=fourier->dim().x();
    dimfy=fourier->dim().y();
    dimfz=fourier->dim().z();

    iterador.begin();
//    std::cout << "Al principio: "<< iterador.get()<< std::endl;
    do
    {
      // frequency distance
      vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

      if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
      if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}

      g_exp = expf( -(  (float)((pos.x()*pos.x())) / (sigma.x()) +
              (float)((pos.y()*pos.y())) / (sigma.y()) +
              (float)((pos.z()*pos.z())) / (sigma.z())) - fttw_cte );

     // multiply gaussian function
      value= iterador.get()*complex<float>(g_exp,0.0);
     iterador.set(value);

    }while(iterador.next());

    // back FFT transform
    p=fftwf_plan_dft_c2r_3d(dimz,dimy,dimx,
                            fout,vol,FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    delete fourier;
  }

   // Mon made (13/2/2009)
   // Map filtering in Fourier Space.
   // Needs: two "fftw plans" (direct and inverse FFT), a Fourier map, the resolution,
   // and the fourier normalization factor "fttw_cte".
   // Outputs: the filtered map will be placed in the map pointed by the "fftw plan"
   void GaussFilter(void *ppv, void *ppv2, vlVolume *fourier, float res, float fttw_cte)
   // void GaussFilter(fftwf_plan *pp, fftwf_plan *pp2, vlVolume *fourier, float res, float fttw_cte)
   {
     fftwf_plan *pp = (fftwf_plan *) ppv;
     fftwf_plan *pp2 = (fftwf_plan *) ppv2;

     int dimfx,dimfy,dimfz;

     fftwf_execute(*pp); // Perform FT with FFTW
     vlPoint3f sigma;
     float g_exp;
     complex<float> value;

     // Calculo de la posicion de  las frecuencias limite en cada una de las direcciones
     // en el espacio de Fourier
     // Para cada direccion = mitad de la longitud del volumen partido por la resolucion
     sigma.x( (float)(fourier->units().x() * fourier->dim().x()*2) / res );
     sigma.y( (float)(fourier->units().y() * fourier->dim().y()) / res );
     sigma.z( (float)(fourier->units().z() * fourier->dim().z()) / res );

     //Calculo del sigma para cada direccion, depende de la posicion de la frecuencia
     sigma.x(sigma.x()*sigma.x());
     sigma.y(sigma.y()*sigma.y());
     sigma.z(sigma.z()*sigma.z());

     // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
     //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
     vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);

     // std::cout <<" Dimensiones del volumen: " << fourier->dim() << std::endl;
     // std::cout <<" Sigmas del volumen: " << sigma << std::endl;
     // std::cout <<" Numero de voxels del volumen: " << fourier->voxelCount() << std::endl;
     // std::cout <<" centro: " << center << std::endl;

     // Creacion de un iterador para recorrer todo el espacio de fourier
     vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

     dimfx=fourier->dim().x();
     dimfy=fourier->dim().y();
     dimfz=fourier->dim().z();

     iterador.begin();
   //    std::cout << "Al principio: "<< iterador.get()<< std::endl;
     do
     {
       // frequency distance
       vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

       if((pos.y()+center.y())>=dimfy)
       	pos.y(pos.y()-dimfy);
       if((pos.z()+center.z())>=dimfz)
       	pos.z(pos.z()-dimfz);

       g_exp = expf( -(  (float)((pos.x()*pos.x())) / (sigma.x()) +
               (float)((pos.y()*pos.y())) / (sigma.y()) +
               (float)((pos.z()*pos.z())) / (sigma.z())) - fttw_cte );

       // multiply gaussian function
       value= iterador.get()*complex<float>(g_exp,0.0);
       iterador.set(value); // modifies fourier volume

     }while(iterador.next());

     fftwf_execute(*pp2); // back FFT transform

     iterador.~vlVolIter<complex<float>, vlLayout::Linear>();
   //  delete fourier;
   }

  /*Filtrado en paso banda de un volumen tridimensional mediante un filtro gaussiano.
    Entradas:
      iter: Iterador del volumen (mejor puntero a volumen?)
      res1: resolucion a la que se filtra el volumen en paso bajo. Este parametro esta limitado
      por los siguientes valores:
      Resolucion minima= Tama�o de Voxel
      Resolucion minima= Maxima longitud del volumen * Tama�o de voxel /2

      Donde maxima longitud es la longitud de la mayor diagonal que cruza el volumen

      res2: Idem en paso alto
      */
   void GaussBandFilter(vlVolume *vlm,
                float res1, float res2)
  {
   vlVolIter<float, vlLayout::Linear> iter(vlm);
    float *vol;
    int dimx,dimy,dimz;
    int dimfx,dimfy,dimfz;
    vlVolDataLayoutBase<float> *dataLayout;

    fftwf_complex *fout;
    fftwf_plan p;
    vlVolume *fourier;

    // Obtencion del puntero al volumen
    dataLayout=iter.getVolume();

    //Comprobacion de que el tama�o de la resolucion es adecuado
    float diagonal;
    diagonal = sqrt( pow((float)dataLayout->dim().x(),2) + pow((float)dataLayout->dim().y(),2) + pow((float)dataLayout->dim().z(),2));

    diagonal = diagonal*(float)(dataLayout->units().x())/2.0;

    if(res1 > res2)
    {
      fprintf(stderr, "GaussianBandFilter ERROR: First resolution must be lower than second.\n");
      //return ;
    }
    if( res1<dataLayout->units().x() || res1>diagonal || res2<dataLayout->units().x() || res2>diagonal)
    {
    	fprintf(stderr,"GaussianBandFilter ERROR: Resolution size out of limits.\n");
    	fprintf(stderr, "Diagonal: %f resolucion: %f\n",diagonal,res2);
      //return ;
    }
    // Fin de la comprobacion


    // Obtencion del puntero al array de datos del volumen
    vol=dataLayout->getVoxelPtr(vlPoint3ui(0,0,0));

    //Dimensiones del array de datos
    dimx=dataLayout->dim().x();
    dimy=dataLayout->dim().y();
    dimz=dataLayout->dim().z();

    //Creacion del volumen donde se alamacenara la transformada de Fourier
    fourier= new vlVolume(vlDim((dimx/2)+1,dimy,dimz),Complexf,dataLayout->units());
    //Obtencion del puntero al array de datos de la transformada de Fourier
    fout=(fftwf_complex*)((vlVolDataLayoutBase<complex<float> >*)fourier->volData())->getVoxelPtr(vlPoint3ui(0,0,0));

    //Creacion del plan para realizar la transformada. Se le pasa como argumentos
    // las dimensiones del array de entrada asi como las direcciones de memoria de los
    // array de entrada y salida
    p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
                            vol,fout,FFTW_ESTIMATE);

    // Calculo de la transformada de Fourier
    fftwf_execute(p);
    fftwf_destroy_plan(p);




    /*EN LA VARIABLE FOURIER ESTA EL VOLUMEN EN EL ESPACIO DE FOURIER*/
    vlPoint3f radio, sigma, radio2, sigma2;
    float valuex,valuey,valuez;
    complex<float> value;

    // Calculo de laposicion de  las frecuencias maximas y minimas en cada una de las direcciones
    // en el espacio de Fourier
    // Para cada direccion = mitad de la longitud del volumen partido por la resolucion
    radio.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res1*2) );
    radio.y( (float)(fourier->units().y() * fourier->dim().y()) / (res1*2) );
    radio.z( (float)(fourier->units().z() * fourier->dim().z()) / (res1*2) );

    radio2.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res2*2) );
    radio2.y( (float)(fourier->units().y() * fourier->dim().y()) / (res2*2) );
    radio2.z( (float)(fourier->units().z() * fourier->dim().z()) / (res2*2) );


    //Calculo del sigma para cada direccion y para cada resolucion, depende de la posicion de la frecuencia
    sigma.x(radio.x()/(sqrt(LN2)));
    sigma.y(radio.y()/(sqrt(LN2)));
    sigma.z(radio.z()/(sqrt(LN2)));

    sigma2.x(radio2.x()/(sqrt(LN2)));
    sigma2.y(radio2.y()/(sqrt(LN2)));
    sigma2.z(radio2.z()/(sqrt(LN2)));

//    std::cout <<" Dimensiones del volumen: " << fourier->dim() << std::endl;
//    std::cout <<" Sigmas del volumen: " << sigma << std::endl;
//    std::cout <<" Numero de voxels del volumen: " << fourier->voxelCount() << std::endl;

    // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
    //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
    vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);


//    std::cout <<" centro: " << center << std::endl;

    // Creacion de un iterador para recorrer todo el espacioo de fourier
    vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

    dimfx=fourier->dim().x();
    dimfy=fourier->dim().y();
    dimfz=fourier->dim().z();

//    int cont=0;
    iterador.begin();
    //    std::cout << "Al principio: "<< iterador.get()<< std::endl;
    do
    {
      //Calculo de la posicion de cada voxel teniendo en cuenta que los cuadrantes estan
      //desordenados
      vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

      if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
      if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}

      //Paso del filtro de paso bajo
      valuex= (float)((pos.x()*pos.x())) / ((sigma.x()*sigma.x()+(sigma.y()*sigma.y())+(sigma.z()*sigma.z())) );
      valuey= (float)((pos.y()*pos.y())) / ((sigma.x()*sigma.x()+(sigma.y()*sigma.y())+(sigma.z()*sigma.z())) );
      valuez= (float)((pos.z()*pos.z())) / ((sigma.x()*sigma.x()+(sigma.y()*sigma.y())+(sigma.z()*sigma.z())) );
      value= iterador.get()*complex<float>(expf(-0.5*(valuex+valuey+valuez)),0);
      iterador.set(value);

      // Paso del filtro de paso alto
      valuex= (float)((pos.x()*pos.x())) / ((sigma2.x()*sigma2.x()+(sigma2.y()*sigma2.y())+(sigma2.z()*sigma2.z())) );
      valuey= (float)((pos.y()*pos.y())) / ((sigma2.x()*sigma2.x()+(sigma2.y()*sigma2.y())+(sigma2.z()*sigma2.z())) );
      valuez= (float)((pos.z()*pos.z())) / ((sigma2.x()*sigma2.x()+(sigma2.y()*sigma2.y())+(sigma2.z()*sigma2.z())) );
      value= iterador.get()*complex<float>(expf(-0.5*(valuex+valuey+valuez)),0);
      iterador.set(iterador.get()-value);


    }while(iterador.next());


//std::cout << "Numero de voxels eliminados: "<< cont<< std::endl;


    //Fin del filtro

    //Recuperacion del espacio real deshaciendo la transformada de Fourier
    p=fftwf_plan_dft_c2r_3d(dimz,dimy,dimx,
                            fout,vol,FFTW_ESTIMATE);

    fftwf_execute(p);
    fftwf_destroy_plan(p);
    mul(vlm,1.0/(dimx*dimy*dimz));

  }





 /*Filtrado de un volumen tridimensional mediante un filtro de Butterworth.
    Entradas:
      iter: Iterador del volumen (mejor puntero a volumen?)
      n: orden del filtro. A mayor orden mas pronunciado es el cambio en el limite de
      frecuencia
      res: resolucion a la que se filtra el volumen. Este parametro esta limitado
      por los siguientes valores:
      Resolucion minima= Tama�o de Voxel
      Resolucion minima= Maxima longitud del volumen * Tama�o de voxel /2

      Donde maxima longitud es la longitud de la mayor diagonal que cruza el volumen

      option: si es 0 el filtrado sera de paso bajo, en caso contrario de paso alto
      */
   void ButterFilter(vlVolume *vlm,
                float res, int n, int option)
   {

   vlVolIter<float, vlLayout::Linear> iter(vlm);
    float *vol;
    int dimx,dimy,dimz;
    int dimfx,dimfy,dimfz;
    vlVolDataLayoutBase<float> *dataLayout;

    fftwf_complex *fout;
    fftwf_plan p;
    vlVolume *fourier;

    // Obtencion del puntero al volumen
    dataLayout=iter.getVolume();

    //Comprobacion de que el tama�o de la resolucion es adecuado
    float diagonal;
    diagonal = sqrt( pow((float)dataLayout->dim().x(),2) + pow((float)dataLayout->dim().y(),2) + pow((float)dataLayout->dim().z(),2));

    diagonal = diagonal*(float)(dataLayout->units().x())/2.0;

    if( res<dataLayout->units().x() || res>diagonal)
    {
    	fprintf(stderr,"ERROR: Resolution size out of limits.\n");
      return ;
    }
    // Fin de la comprobacion


    // Obtencion del puntero al array de datos del volumen
    vol=dataLayout->getVoxelPtr(vlPoint3ui(0,0,0));

    //Dimensiones del array de datos
    dimx=dataLayout->dim().x();
    dimy=dataLayout->dim().y();
    dimz=dataLayout->dim().z();

    //Creacion del volumen donde se alamacenara la transformada de Fourier
    fourier= new vlVolume(vlDim((dimx/2)+1,dimy,dimz),Complexf,dataLayout->units());
    //Obtencion del puntero al array de datos de la transformada de Fourier
    fout=(fftwf_complex*)((vlVolDataLayoutBase<complex<float> >*)fourier->volData())->getVoxelPtr(vlPoint3ui(0,0,0));

    //Creacion del plan para realizar la transformada. Se le pasa como argumentos
    // las dimensiones del array de entrada asi como las direcciones de memoria de los
    // array de entrada y salida
    p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
                            vol,fout,FFTW_ESTIMATE);

    // Calculo de la transformada de Fourier
    fftwf_execute(p);
    fftwf_destroy_plan(p);





    /*EN LA VARIABLE FOURIER ESTA EL VOLUMEN EN EL ESPACIO DE FOURIER*/

    vlPoint3f radio;
    float valuex,valuey,valuez;
    complex<float> value;

    // Calculo de laposicion de  las frecuencias limite en cada una de las direcciones
    // en el espacio de Fourier
    // Para cada direccion = mitad de la longitud del volumen partido por la resolucion
    radio.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res*2) *0.75);
    radio.y( (float)(fourier->units().y() * fourier->dim().y()) / (res*2) *0.75);
    radio.z( (float)(fourier->units().z() * fourier->dim().z()) / (res*2) *0.75);


//    std::cout <<" Dimensiones del volumen: " << fourier->dim() << std::endl;
//    std::cout <<" Sigmas del volumen: " << sigma << std::endl;
//    std::cout <<" Numero de voxels del volumen: " << fourier->voxelCount() << std::endl;

    // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
    //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
    vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);



    // Creacion de un iterador para recorrer todo el espacioo de fourier
    vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

    dimfx=fourier->dim().x();
    dimfy=fourier->dim().y();
    dimfz=fourier->dim().z();

//    int cont=0;
    iterador.begin();
//    std::cout << "Al principio: "<< iterador.get()<< std::endl;

    do
    {
      //Calculo de la posicion de cada voxel teniendo en cuenta que los cuadrantes estan
      //desordenados
      vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

      if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
      if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}

      //Calculo de la funcion de Butterworth para el voxel y modificacion del valor del
      //voxel en funcion de esta funcion. La funcion de Butterworth depende de las
      //posiciones de la frecuencia maxima en cada direccion
      valuex= (float)((pos.x()*pos.x())) / ((radio.x()*radio.x()+(radio.y()*radio.y())+(radio.z()*radio.z())) );
      valuey= (float)((pos.y()*pos.y())) / ((radio.x()*radio.x()+(radio.y()*radio.y())+(radio.z()*radio.z())) );
      valuez= (float)((pos.z()*pos.z())) / ((radio.x()*radio.x()+(radio.y()*radio.y())+(radio.z()*radio.z())) );
      valuex= pow((float)valuex+valuey+valuez,n);
      valuex=1.0+(sqrt(2.0)-1.0)*valuex;
      valuex=1.0/valuex;
      value= iterador.get()*complex<float>(valuex,0);
      // Si el filtro es de paso bajo, el valor calculado es asignado al voxel
      if(option==0)
        iterador.set(value);
      //Si el filtro es de paso alto, el valor es restado al valor del voxel
      else
        iterador.set(iterador.get()-value);

    }while(iterador.next());


//std::cout << "Numero de voxels eliminados: "<< cont<< std::endl;


    /*Fin del filtro*/

    //Recuperacion del espacio real deshaciendo la transformada de Fourier
    p=fftwf_plan_dft_c2r_3d(dimz,dimy,dimx,
                            fout,vol,FFTW_ESTIMATE);

    fftwf_execute(p);
    fftwf_destroy_plan(p);
    mul(vlm,1.0/(dimx*dimy*dimz));
  }



 /*Filtrado en paso banda de un volumen tridimensional mediante un filtro de Butterworth.
    Entradas:
      iter: Iterador del volumen (mejor puntero a volumen?)
      n: orden del filtro. A mayor orden mas pronunciado es el cambio en el limite de
      frecuencia
      res1: resolucion a la que se filtra el volumen en paso bajo. Este parametro esta limitado
      por los siguientes valores:
      Resolucion minima= Tama�o de Voxel
      Resolucion minima= Maxima longitud del volumen * Tama�o de voxel /2

      Donde maxima longitud es la longitud de la mayor diagonal que cruza el volumen

      res2: idem para resolucion en paso alto
      */
   void ButterBandFilter(vlVolume *vlm,
                float res1,float res2, int n)
   {

    vlVolIter<float, vlLayout::Linear> iter(vlm);
    float *vol;
    int dimx,dimy,dimz;
    int dimfx,dimfy,dimfz;
    vlVolDataLayoutBase<float> *dataLayout;

    fftwf_complex *fout;
    fftwf_plan p;
    vlVolume *fourier;

    // Obtencion del puntero al volumen
    dataLayout=iter.getVolume();

    //Comprobacion de que el tama�o de la resolucion es adecuado
    float diagonal;
    diagonal = sqrt( pow((float)dataLayout->dim().x(),2) + pow((float)dataLayout->dim().y(),2) + pow((float)dataLayout->dim().z(),2));

    diagonal = diagonal*(float)(dataLayout->units().x())/2.0;

    if(res2< res1)
    {
    	fprintf(stderr,"ButterBandFilter ERROR: First resolution must be lower than second.\n");
      return ;
    }
    if( res1<dataLayout->units().x() || res1>diagonal || res2<dataLayout->units().x() || res2>diagonal)
    {
    	fprintf(stderr,"ButterGaussianBandFilter ERROR: Resolution size out of limits.\n");
    	return ;
    }

    // Fin de la comprobacion


    // Obtencion del puntero al array de datos del volumen
    vol=dataLayout->getVoxelPtr(vlPoint3ui(0,0,0));

    //Dimensiones del array de datos
    dimx=dataLayout->dim().x();
    dimy=dataLayout->dim().y();
    dimz=dataLayout->dim().z();

    //Creacion del volumen donde se alamacenara la transformada de Fourier
    fourier= new vlVolume(vlDim((dimx/2)+1,dimy,dimz),Complexf,dataLayout->units());
    //Obtencion del puntero al array de datos de la transformada de Fourier
    fout=(fftwf_complex*)((vlVolDataLayoutBase<complex<float> >*)fourier->volData())->getVoxelPtr(vlPoint3ui(0,0,0));

    //Creacion del plan para realizar la transformada. Se le pasa como argumentos
    // las dimensiones del array de entrada asi como las direcciones de memoria de los
    // array de entrada y salida
    p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
                            vol,fout,FFTW_ESTIMATE);

    // Calculo de la transformada de Fourier
    fftwf_execute(p);
    fftwf_destroy_plan(p);



    //EN LA VARIABLE FOURIER ESTA EL VOLUMEN EN EL ESPACIO DE FOURIER

    vlPoint3f radio, radio2;
    float valuex,valuey,valuez;
    complex<float> value;

    // Calculo de laposicion de  las frecuencias maxima y minima en cada una de las direcciones
    // en el espacio de Fourier
    // Para cada direccion = mitad de la longitud del volumen partido por la resolucion
    radio.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res1*2) );
    radio.y( (float)(fourier->units().y() * fourier->dim().y()) / (res1*2) );
    radio.z( (float)(fourier->units().z() * fourier->dim().z()) / (res1*2) );

    radio2.x( (float)(fourier->units().x() * fourier->dim().x()*2) / (res2*2) );
    radio2.y( (float)(fourier->units().y() * fourier->dim().y()) / (res2*2) );
    radio2.z( (float)(fourier->units().z() * fourier->dim().z()) / (res2*2) );



//    std::cout <<" Dimensiones del volumen: " << fourier->dim() << std::endl;
//    std::cout <<" Sigmas del volumen: " << sigma << std::endl;
//    std::cout <<" Numero de voxels del volumen: " << fourier->voxelCount() << std::endl;

    // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
    //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
    vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);


//    std::cout <<" centro: " << center << std::endl;

    // Creacion de un iterador para recorrer todo el espacioo de fourier
    vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

    dimfx=fourier->dim().x();
    dimfy=fourier->dim().y();
    dimfz=fourier->dim().z();

//    int cont=0;
    iterador.begin();
//    std::cout << "Al principio: "<< iterador.get()<< std::endl;

    do
    {
      //Calculo de la posicion de cada voxel teniendo en cuenta que los cuadrantes estan
      //desordenados
      vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

      if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
      if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}

      //Paso del filtro en paso bajo
      valuex= (float)((pos.x()*pos.x())) / ((radio.x()*radio.x()+(radio.y()*radio.y())+(radio.z()*radio.z())) );
      valuey= (float)((pos.y()*pos.y())) / ((radio.x()*radio.x()+(radio.y()*radio.y())+(radio.z()*radio.z())) );
      valuez= (float)((pos.z()*pos.z())) / ((radio.x()*radio.x()+(radio.y()*radio.y())+(radio.z()*radio.z())) );
      valuex= pow((float)valuex+valuey+valuez,n);
      valuex=1.0+(sqrt(2.0)-1.0)*valuex;
      valuex=1.0/valuex;
      value= iterador.get()*complex<float>(valuex,0);
      iterador.set(value);

      //Paso del filtro en paso alto
      valuex= (float)((pos.x()*pos.x())) / ((radio2.x()*radio2.x()+(radio2.y()*radio2.y())+(radio2.z()*radio2.z())) );
      valuey= (float)((pos.y()*pos.y())) / ((radio2.x()*radio2.x()+(radio2.y()*radio2.y())+(radio2.z()*radio2.z())) );
      valuez= (float)((pos.z()*pos.z())) / ((radio2.x()*radio2.x()+(radio2.y()*radio2.y())+(radio2.z()*radio2.z())) );
      valuex= pow((float)valuex+valuey+valuez,n);
      valuex=1.0+(sqrt(2.0)-1.0)*valuex;
      valuex=1.0/valuex;
      value= iterador.get()*complex<float>(valuex,0);
      iterador.set(iterador.get()-value);


    }while(iterador.next());

//std::cout << "Numero de voxels eliminados: "<< cont<< std::endl;


    //Fin del filtro

    //Recuperacion del espacio real deshaciendo la transformada de Fourier
    p=fftwf_plan_dft_c2r_3d(dimz,dimy,dimx,
                            fout,vol,FFTW_ESTIMATE);

    fftwf_execute(p);
    fftwf_destroy_plan(p);
   mul(vlm,1.0/(dimx*dimy*dimz));
  }

/*Paso del filtro Laplaciano en el espacio de Fourier a un volumen*/
void laplacianFilter(vlVolume *vlm)
 {
   vlVolIter<float, vlLayout::Linear> iter(vlm);
   float *vol;
   int dimx,dimy,dimz;
   int dimfx,dimfy,dimfz;
   vlVolDataLayoutBase<float> *dataLayout;

   fftwf_complex *fout;
   fftwf_plan p;
   vlVolume *fourier;

   // Obtencion del puntero al volumen
   dataLayout=iter.getVolume();

   // Obtencion del puntero al array de datos del volumen
   vol=dataLayout->getVoxelPtr(vlPoint3ui(0,0,0));

   //Dimensiones del array de datos
   dimx=dataLayout->dim().x();
   dimy=dataLayout->dim().y();
   dimz=dataLayout->dim().z();

   //Creacion del volumen donde se alamacenara la transformada de Fourier
   fourier= new vlVolume(vlDim((dimx/2)+1,dimy,dimz),Complexf,dataLayout->units());
   //Obtencion del puntero al array de datos de la transformada de Fourier
   fout=(fftwf_complex*)((vlVolDataLayoutBase<complex<float> >*)fourier->volData())->getVoxelPtr(vlPoint3ui(0,0,0));

   //Creacion del plan para realizar la transformada. Se le pasa como argumentos
   // las dimensiones del array de entrada asi como las direcciones de memoria de los
   // array de entrada y salida
   p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
                           vol,fout,FFTW_ESTIMATE);

   // Calculo de la transformada de Fourier
   fftwf_execute(p);
   fftwf_destroy_plan(p);



   //En fourier esta el volumen en el espacio de fourier

   float valuef;
   complex<float> value;


   // Calculo del punto central teniendo en cuenta que solo se tiene la mitad del
   //volumen en Fourier y que se tiene en cuadrantes desordenados (FFT)
   vlPoint3ui center(0,fourier->dim().y()/2,fourier->dim().z()/2);


//    std::cout <<" centro: " << center << std::endl;

   // Creacion de un iterador para recorrer todo el espacio de fourier
   vlVolIter<complex<float>, vlLayout::Linear> iterador(fourier);

   dimfx=fourier->dim().x();
   dimfy=fourier->dim().y();
   dimfz=fourier->dim().z();

   // Modificacion para cada voxel
//    int cont=0;
   iterador.begin();
//    std::cout << "Al principio: "<< iterador.get()<< std::endl;
   while(!iterador.end())
   {
     //Calculo de la posicion de cada voxel teniendo en cuenta que los cuadrantes estan
     //desordenados
     vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

     if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
     if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}

     // Modificacion del valor del voxel en funcion del filtro Laplaciano
     valuef= -(pos.x()*pos.x()+pos.y()*pos.y()+pos.z()*pos.z());
     value= iterador.get()*complex<float>(valuef,0);
     iterador.set(value);

//      if(valuex<=0.70710678375225)
//        cont++;

     iterador.next();

   }

   vlPoint3f pos(iterador.pos().x(),iterador.pos().y(),iterador.pos().z());

   if((pos.y()+center.y())>=dimfy){pos.y(pos.y()-dimfy);}
   if((pos.z()+center.z())>=dimfz){pos.z(pos.z()-dimfz);}


     //Calculo de la funcion de Butterworth para el voxel y modificacion del valor del
     //voxel en funcion de esta funcion. La funcion de Butterworth depende de las
     //posiciones de la frecuencia maxima en cada direccion
     valuef= -(pos.x()*pos.x()+pos.y()*pos.y()+pos.z()*pos.z());
     value= iterador.get()*complex<float>(valuef,0);
     iterador.set(value);

//      if(valuex<=0.70710678375225)
//        cont++;


//std::cout << "Numero de voxels eliminados: "<< cont<< std::endl;


   //Fin del filtro

   //Recuperacion del espacio real deshaciendo la transformada de Fourier
   p=fftwf_plan_dft_c2r_3d(dimz,dimy,dimx,
                           fout,vol,FFTW_ESTIMATE);

   fftwf_execute(p);
   fftwf_destroy_plan(p);
  mul(vlm,1.0/(dimx*dimy*dimz));
 }


}
