/***************************************************************************
                          fmat.cpp  -  description
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


 #include "floatops.h"
 namespace FOPS
 {
// ***********OPERACIONES MATEMATICAS BASICAS********************

  /*Operacion de suma de una constante en todos los voxels de un volumen
    iter: Iterador sobre el volumen
    value: valor constante
  */
  void add(vlVolume *vol, float value)
  {
     vlVolIter<float, vlLayout::Linear> iter(vol);
    iter.begin();
    while(!iter.end())
    {
      iter.set(iter.get()+value);
      iter.next();
    }
      iter.set(iter.get()+value);
   }

/*Operacion de suma de dos volumenes
  */
  void add(vlVolume *vol, vlVolume *vol2)
  {

    vlVolIter<float, vlLayout::Linear> iter(vol);
    vlVolIter<float, vlLayout::Linear> iter2(vol2);
    iter.begin();
    iter2.begin();

    while(!iter.end())
    {
      iter.set(iter.get()+iter2.get());
      iter.next();
      iter2.next();
    }
      iter.set(iter.get()+iter2.get());
   }


  /*Operacion de multiplicacion de una constante en todos los voxels de un volumen
    iter: Iterador sobre el volumen
    value: valor constante
*/
  void mul(vlVolume *vol, float value)
  {
     vlVolIter<float, vlLayout::Linear> iter(vol);

    iter.begin();
    while(!iter.end())
    {
      iter.set(iter.get()*value);
      iter.next();
    }
      iter.set(iter.get()*value);
   }

 /*Calculo de la sumatoria de la masa total del volumen
  iter: iterador sobre el volumen*/
  float calc_total(vlVolume *vol)
  {
    double sum=0.0;
    double aux;
    float sum2;
    vlVolIter<float, vlLayout::Linear> iter(vol);

    iter.begin();


    do
    {
        if(iter.get()>=0.0)
        {
          aux=(double)iter.get();
          sum+=aux;
        }
    }while(iter.next());


    sum2=(float)sum;
    return sum2;

  }

  /*Calculo de la densidad del volumen
    iter: iterador sobre el volumen*/
    float calc_total_dens(vlVolume *vol)
    {
      double sum=0.0;
      double aux;
      float sum2;
      vlVolIter<float, vlLayout::Linear> iter(vol);

      iter.begin();

      do
      {
          if(iter.get()>=0.0)
          {
            aux=(double)iter.get();
            sum+=aux;
          }
      }while(iter.next());

      sum2=(float)sum/(vol->units().x()*vol->units().y()*vol->units().z() );

      return sum2;

    }




	float calc_total_posneg(vlVolume *vol)
  {
    double sum=0.0;
    double aux;
    float sum2;
    vlVolIter<float, vlLayout::Linear> iter(vol);

    iter.begin();


    do
    {
          aux=(double)iter.get();
          sum+=aux;
    }while(iter.next());


    sum2=(float)sum;
    return sum2;

  }

  /*Calculo del valor minimo de masa*/
  float min_mass(vlVolume *vol)
  {
    vlVolIter<float, vlLayout::Linear> iter(vol);
    iter.begin();
    float min=iter.get();
    do
    {
      iter.next();

      if(iter.get()<min)
        min=iter.get();

    }while(!iter.end());


    return min;
  }

  /*Calculo del valor maximo de masa*/
  float max_mass(vlVolume *vol)
  {
    vlVolIter<float, vlLayout::Linear> iter(vol);
    iter.begin();
    float min=iter.get();
    do
    {
      iter.next();

      if(iter.get()>min)
        min=iter.get();

    }while(!iter.end());

   return min;
  }

  /*Calculo del valor medio de los voxels de un volumen*/
  float calc_average(vlVolume *vol)
  {
    vlVolIter<float, vlLayout::Linear> iter(vol);
    float total=0,i=0;
    iter.begin();
    do{
      i++;
      total+=iter.get();
    }while (iter.next());
    return total/i;
  }

  /*Calculo de la desviacion media de los valores de los voxels de un volumen*/
  float sigma(vlVolume *vol)
  {
    vlVolIter<float, vlLayout::Linear>  iter(vol);
    float average,sigma=0,i=0;
    average=calc_average(vol);
    iter.begin();
    do{
      sigma+=pow((float)iter.get()-average,2);
      i++;
    }while(iter.next());
    return (sqrt(sigma/i));
  }

  /*Calculo de la normal de los valores de los voxels de un volumen*/
  float normal(vlVolume *vol)
  {
    vlVolIter<float, vlLayout::Linear>  iter(vol);
    float norm=0,i=0;

    iter.begin();
    do{
      norm+=pow((float)iter.get(),2);
      i++;
    }while(iter.next());
    return (sqrt(norm/i));
  }

  int max_lenght(vlVolume *vol)
  {
    vlPoint3ui maxCorner(0,0,0);
    vlPoint3ui center(vol->dim().x()/2,vol->dim().y()/2,vol->dim().z()/2);
    int aux2,max=0;
    float aux;

    vlVolIter<float, vlLayout::Linear> iter(vol);

    iter.begin();
    do{
      if(iter.get()>0)
      {


        aux=sqrt((float)(pow((float)center.x()-iter.pos().x(),2)+pow((float)center.y()-iter.pos().y(),2)+pow((float)center.z()-iter.pos().z(),2)));
        aux2=ceil(aux);
        if(aux2>max)
        {
          max=aux2;
        }
      }
  }while(iter.next());
  return max;

  }

  void center_masses(vlVolume *vol, float *cx,float *cy, float *cz)
  {

    double aux,mass=0.0,auxX=0.0,auxY=0.0,auxZ=0.0;
    vlVolIter<float, vlLayout::Linear> iter(vol);


    iter.begin();
    do{
      aux=(double)iter.get();
     if(aux>0)
     {
      auxX+=aux * (double)iter.pos().x();
      auxY+=aux * (double)iter.pos().y();
      auxZ+=aux * (double)iter.pos().z();
      mass+=aux;
     }
    }while(iter.next());

    *cx=(float)(auxX/mass);
    *cy=(float)(auxY/mass);
    *cz=(float)(auxZ/mass);

  }

 void center_geometric(vlVolume *vol, float *cx,float *cy, float *cz)
  {

    double aux,mass=0.0,auxX=0.0,auxY=0.0,auxZ=0.0;
    vlVolIter<float, vlLayout::Linear> iter(vol);


    iter.begin();
    do{
      aux=(double)iter.get();
     if(aux>0)
     {
      auxX+=(double)iter.pos().x();
      auxY+=(double)iter.pos().y();
      auxZ+=(double)iter.pos().z();
      mass++;
     }
    }while(iter.next());

    *cx=(float)(auxX/mass);
    *cy=(float)(auxY/mass);
    *cz=(float)(auxZ/mass);

 }

 void center_vol(vlVolume *vol, float *center)
  {
		center[0]=(float)(vol->dim().x())/2.0 - 0.5;
		center[1]=(float)(vol->dim().y())/2.0 - 0.5;
		center[2]=(float)(vol->dim().z())/2.0 - 0.5;
  }


  float radiusMax(vlVolume *vol)
  {
    float cx,cy,cz;
    float dist,distMax=0;
    center_masses(vol,&cx,&cy,&cz);
    vlVolIter<float, vlLayout::Linear> iter(vol);

    iter.begin();
    do{
      if(iter.get()>0)
      {
      dist=(iter.pos().x()-cx)*(iter.pos().x()-cx)
         +(iter.pos().y()-cy)*(iter.pos().y()-cy)
         +(iter.pos().z()-cz)*(iter.pos().z()-cz);
      if(dist>distMax)
        distMax=dist;
      }
    }while(iter.next());
    return sqrt(distMax);
  }

  // Mon made (6/2/2009)
  // Computes Average inside a frame (size = frame dimensions = padding)
  float avg_frame(vlVolume *vol, vlDim size)
  {
	  float average=0;
	  int i,j,z;
	  int nvox;
	  int stepx,stepy,stepz;

	  vlStep step=vol->stepping();
	  stepx=step.x();
	  stepy=step.y();
	  stepz=step.z();

	  vlDim dim=vol->dim();

	  float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	  for(i=size.x();i<dim.x()-size.x();i++)
		  for(j=size.y();j<dim.y()-size.y();j++)
			  for(z=size.z();z<dim.z()-size.z();z++)
				  average += *( container + (i*stepx)+(j*stepy)+(z*stepz) ); // 1st vol voxel
	  nvox = (dim.x() - 2*size.x() ) * (dim.y() - 2*size.y() ) * (dim.z() - 2*size.z() );
	  average /= nvox;
	  return average;
  }

  // Mon made (13/9/2010)
  // Computes Average (DOUBLE precision) inside a frame (size = frame dimensions = padding)
  double avg_frameD(vlVolume *vol, vlDim size)
  {
	  double average=0;
	  int i,j,z;
	  int nvox;
	  int stepx,stepy,stepz;

	  vlStep step=vol->stepping();
	  stepx=step.x();
	  stepy=step.y();
	  stepz=step.z();

	  vlDim dim=vol->dim();

	  float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	  for(i=size.x();i<dim.x()-size.x();i++)
		  for(j=size.y();j<dim.y()-size.y();j++)
			  for(z=size.z();z<dim.z()-size.z();z++)
				  average += *( container + (i*stepx)+(j*stepy)+(z*stepz) ); // 1st vol voxel
	  nvox = (dim.x() - 2*size.x() ) * (dim.y() - 2*size.y() ) * (dim.z() - 2*size.z() );
	  average /= nvox;
	  return average;
  }

  // Mon made (6/2/2009)
    // Computes Sigma inside a frame (size = frame dimensions = padding)
    float sig_frame(vlVolume *vol, vlDim size)
    {
  	  float average,sigma;
  	  int i,j,z;
  	  int nvox;
  	  int stepx,stepy,stepz;

  	  average = FOPS::avg_frame(vol,size); // compute average
  	  vlStep step=vol->stepping();
  	  stepx=step.x();
  	  stepy=step.y();
  	  stepz=step.z();
  	  vlDim dim=vol->dim();
  	  float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

  	  for(i=size.x();i<dim.x()-size.x();i++)
  		  for(j=size.y();j<dim.y()-size.y();j++)
  			  for(z=size.z();z<dim.z()-size.z();z++)
  				  sigma += pow( *( container + (i*stepx)+(j*stepy)+(z*stepz) ) - average, 2);

  	  nvox = (dim.x() - 2*size.x() ) * (dim.y() - 2*size.y() ) * (dim.z() - 2*size.z() );
  	  sigma = sqrt(sigma/nvox);
  	  return sigma;
    }

    // Mon made (13/9/2009)
    // Computes Sigma (DOUBLE precision) inside a frame (size = frame dimensions = padding)
    double sig_frameD(vlVolume *vol, vlDim size)
    {
    	double average,sigma=0;
    	int i,j,z;
    	int nvox;
    	int stepx,stepy,stepz;

    	average = FOPS::avg_frameD(vol,size); // compute average
    	vlStep step=vol->stepping();
    	stepx=step.x();
    	stepy=step.y();
    	stepz=step.z();
    	vlDim dim=vol->dim();
    	float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

    	for(i=size.x();i<dim.x()-size.x();i++)
    		for(j=size.y();j<dim.y()-size.y();j++)
    			for(z=size.z();z<dim.z()-size.z();z++)
    				sigma += pow( *( container + (i*stepx)+(j*stepy)+(z*stepz) ) - average, 2);

    	nvox = (dim.x() - 2*size.x() ) * (dim.y() - 2*size.y() ) * (dim.z() - 2*size.z() );
    	sigma = sqrt(sigma/nvox);
    	return sigma;
    }

// ***********FIN OPERACIONES MATEMATICAS BASICAS********************

}
