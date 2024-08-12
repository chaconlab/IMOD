/***************************************************************************
                          fkernel.cpp  -  description
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


#include <string.h>
#include "floatops.h"
#include "vlkernel.h"
#ifdef USE_TBB // Enables Intel's Threading Building Blocks parallel routines
#include "tbb/tbb.h"
#endif


 namespace FOPS
 {
  // ***************FILTROS EN EL DOMINIO REAL***************************************
   float correlation(vlVolume *vol, vlVolume *vol2)
   {
    //int cont=0;
    if(vol->voxelCount() != vol2->voxelCount())
    {
    	fprintf(stderr,"Correlation Error: Diferentes dimensiones\n");
    }

    vlVolIter<float, vlLayout::Linear> iter(vol);
    vlVolIter<float, vlLayout::Linear> iter2(vol2);
    int cont;
    float aux,corrTotal=0;
    float average1,average2,sigma1,sigma2;

    //map 1 average
		iter.begin();
		cont=0;
		average1=0;
		do{
			average1+=iter.get();
      cont++;

		}while(iter.next() );
		average1/=(float)cont;

    //sigma1
    iter.begin();
		cont=0;
		sigma1=0;
		do{
				sigma1+=pow((float)iter.get()-average1,2);
      	cont++;
		}while(iter.next() );
		sigma1=sqrt(sigma1/cont);

    //map 2 average
		iter2.begin();
		average2=0;
		cont=0;
		do{
			average2+=iter2.get();
      cont++;

		}while(iter2.next() );
		average2/=(float)cont;

    //sigma2
    iter2.begin();
		sigma2=0;
		cont=0;
		do{
				sigma2+=pow((float)iter2.get()-average2,2);
      	cont++;
		}while(iter2.next() );
		sigma2=sqrt(sigma2/cont);

    //corr
    iter.begin();
    cont=0;
    do{
      vol2->getVoxel(iter.pos(),aux);

      //if(iter.get()*aux!=0.0)
      //{
      	corrTotal+= iter.get()*aux;
      		cont++;
      //}
    }while(iter.next());
    corrTotal=(corrTotal-(cont*average1*average2))/((cont)*sigma1*sigma2);

		//corrTotal/=(float)(cont);
    return corrTotal;
   }

 float correlation_simple(vlVolume *vol, vlVolume *vol2)
    {
     //int cont=0;
     if(vol->voxelCount() != vol2->voxelCount())
     {
    	 fprintf(stderr,"Correlation Error: Different dimensions\n");
     }

     vlVolIter<float, vlLayout::Linear> iter(vol);
     vlVolIter<float, vlLayout::Linear> iter2(vol2);
     int cont;
     float corrTotal=0;


     //map 1 average
     iter.begin();
     iter2.begin();
      	cont=0;

 		do{
 			corrTotal+=iter.get()*iter2.get();
 			cont++;
 			iter2.next();
 		}while(iter.next() );


 	 corrTotal/=(float)(cont);
     return corrTotal;
    }


// Mon made (3/6/2008)
// Cross.Corr. inside a cubic-frame (nozero=false) or inside the final volume (nozero=true, default)
// size = kernel width
// vol = Target (Final) volume
float correlation_frame(vlVolume *vol2, vlVolume *vol,int size, bool nozero)
{
    if(vol->voxelCount() != vol2->voxelCount())
    {
    	fprintf(stderr,"Correlation Error: dimension mismatch\n");
    }

    float aux,corrTotal=0.0;
    float average1=0,average2=0,sigma1=0,sigma2=0;
    int i,j,z;
    int cont=0;
    float *p_voxel,*p_voxel2;
    float voxel,voxel2;
    int stepx,stepy,stepz;
    int offset;

    vlStep step=vol->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();

    vlDim dim=vol->dim();
    int border= (size-1)/2; // respect to the kernell size

    float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container2= (float*)vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

    if(nozero)
    { // Inside final map (non-zero voxels)
		// Map 1 & 2 average & dot product
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          if( voxel != 0.0 )
	          {
	          	average1 += voxel;
	          	average2 += voxel2;
	          	cont++;
	          	corrTotal += voxel * voxel2;
	          }
	        }
	    average1 /= cont;
	    average2 /= cont;

	    // Map 1 & 2 sigma
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          if( voxel != 0.0 )
	          {
	          	sigma1 += pow(voxel-average1,2);
	          	sigma2 += pow(voxel2-average2,2);
	          }
	        }
	     sigma1 = sqrt(sigma1/cont);
	     sigma2 = sqrt(sigma2/cont);
    }
    else
    { // Inside Cubic Frame
		// Map 1 & 2 average & dot product
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	       	  average1 += voxel;
	       	  average2 += voxel2;
	          corrTotal += voxel * voxel2;
	          cont++;
	        }
	    average1 /= cont;
	    average2 /= cont;

	    // Map 1 & 2 sigma
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          sigma1 += pow(voxel-average1,2);
	          sigma2 += pow(voxel2-average2,2);
	        }
	     sigma1 = sqrt(sigma1/cont);
	     sigma2 = sqrt(sigma2/cont);
    }

	// Normalized cross-correlation computation
    corrTotal = (corrTotal-(cont*average1*average2))/((cont)*sigma1*sigma2);
    return corrTotal;
}

// Mon made (3/6/2008)
 // Cross.Corr. inside a cubic-frame (nozero=false) or inside the final volume (nozero=true, default)
 // size = kernel width
 // vol = Target (Final) volume
 // If we have avg1 and sig1 (from: vol), we should provide them!
 float correlation_frame(vlVolume *vol2, vlVolume *vol, vlDim size, bool nozero, bool volcalc, float volavg, float volsig)
 {
     if(vol->voxelCount() != vol2->voxelCount())
     {
    	 fprintf(stderr,"Correlation Error: Dimension mismatch\n");
     }

     float corrTotal=0.0;
     float average1=0,average2=0,sigma1=0,sigma2=0;
     int i,j,z;
     int cont=0,null1=0,null2=0;
     float *p_voxel,*p_voxel2;
     float voxel,voxel2;
     int stepx,stepy,stepz;
     int offset;

     vlStep step=vol->stepping();
     stepx=step.x();
     stepy=step.y();
     stepz=step.z();

     vlDim dim=vol->dim();
 //    int border= (size-1)/2; // respect to the kernell size

     float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
     float *container2= (float*)vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

     if(nozero)
     { // Inside final map (>zero voxels)
 		// Map 1 & 2 average & dot product
     	if(volcalc) // then compute BOTH averages & sigmas
     	{
 		    for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          if( voxel > 0.0 ) // "vol" >0 mask
 		          {
 		            voxel2= *(container2+offset); // 2nd vol voxel
 	          	    average2 += voxel2;
 		          	average1 += voxel;
 		          	corrTotal += voxel * voxel2;
 		          	cont++;
 		          }
 		        }
 		    average1 /= cont;
 		    average2 /= cont;

 		    // Map 1 & 2 sigma
 		    for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          if( voxel > 0.0 )
 		          {
 			        voxel2= *(container2+offset); // 2nd vol voxel
 		          	sigma1 += pow(voxel-average1,2);
 		          	sigma2 += pow(voxel2-average2,2);
 		          }
 		        }
 		     sigma1 = sqrt(sigma1/cont);
 		     sigma2 = sqrt(sigma2/cont);
     	}
     	else // compute SINGLE average & sigma (vol2's avg.)
     	{
     		average1 = volavg;
     		sigma1 = volsig;

     		for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          if( voxel > 0.0 ) // "vol" >0 mask
 		          {
 		            voxel2= *(container2+offset); // 2nd vol voxel
 //		          	average1 += voxel;
 		          	average2 += voxel2;
 		          	corrTotal += voxel * voxel2;
 	          	    cont++;
 		          }
 		        }
 //		    average1 /= cont;
 		    average2 /= cont;

 		    // Map 1 & 2 sigma
 		    for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          if( voxel > 0.0 )
 		          {
 			        voxel2= *(container2+offset); // 2nd vol voxel
 //		          	sigma1 += pow(voxel-average1,2);
 		          	sigma2 += pow(voxel2-average2,2);
 		          }
 		        }
 //		     sigma1 = sqrt(sigma1/cont);
 		     sigma2 = sqrt(sigma2/cont);
     	}
     }
     else
     { // Inside CUBIC FRAME

  	// Number of voxels (without taking into account the frame)
      cont = (dim.x() - 2*size.x() ) * (dim.y() - 2*size.y() ) * (dim.z() - 2*size.z() );

      // Map 1 & 2 average & dot product
     	if(volcalc) // then compute BOTH averages & sigmas
     	{
 		    for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          voxel2= *(container2+offset); // 2nd vol voxel
 		          if(voxel != 0)
 		        	  average1 += voxel;
 		          else
 		        	  null1++; // counting null voxels (speeds-up sigma)
 		          if(voxel2 != 0)
 		          {
 		        	  average2 += voxel2;
 		        	  corrTotal += voxel * voxel2;
 		          }
 		          else
		        	  null2++; // counting null voxels (speeds-up sigma)
//   		          cont++;
 		        }
 		    average1 /= cont;
 		    average2 /= cont;
//   		    average2 = 0.220978; //

 		    // Map 1 & 2 sigma
 		    for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          if(voxel != 0)
 		        	  sigma1 += pow(voxel-average1,2);
	        	  voxel2= *(container2+offset); // 2nd vol voxel
 		          if(voxel2 != 0)
 		        	  sigma2 += pow(voxel2-average2,2);
 		        }
 		     sigma1 += null1 * pow(average1,2);
 		     sigma1 =  sqrt(sigma1/cont);
 		     sigma2 += null2 * pow(average2,2);
 		     sigma2 =  sqrt(sigma2/cont);
//   		     sigma2 = 0.694558; //
//   		     sigma2 = 0.8;
//   		     sigma2 = 2;
     	}
     	else // compute SINGLE average & sigma (vol2's avg.)
     	{
     		average1 = volavg;
     		sigma1 = volsig;

     		for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 		          voxel= *(container+offset); // 1st vol voxel
 		          voxel2= *(container2+offset); // 2nd vol voxel
 //		          average1 += voxel;
 		          if(voxel2 != 0)
 		          {
 		        	  average2 += voxel2;
 		        	  corrTotal += voxel * voxel2;
 		          }
 		          else
 		        	  null2++;
//             	  cont++;
 		        }
 //		    average1 /= cont;
 		    average2 /= cont;

 		    // Map 1 & 2 sigma
 		    for(i=size.x();i<dim.x()-size.x();i++)
 		      for(j=size.y();j<dim.y()-size.y();j++)
 		        for(z=size.z();z<dim.z()-size.z();z++)
 		        {
 		          offset = (i*stepx)+(j*stepy)+(z*stepz);
 			      voxel2= *(container2+offset); // 2nd vol voxel
 //		          sigma1 += pow(voxel-average1,2);
 			      if(voxel2 != 0)
 			    	  sigma2 += pow(voxel2-average2,2);

 		        }
 //		     sigma1 = sqrt(sigma1/cont);
		     sigma2 += null2 * pow(average2,2);
 		     sigma2 = sqrt(sigma2/cont);
     	}
     }

 	// Normalized cross-correlation computation
 //    printf("Msg(correlation_frame): dot= %f   avg1= %f   avg2= %f   sig1= %f   sig2= %f",corrTotal,average1,average2,sigma1,sigma2);
     corrTotal = (corrTotal-(cont*average1*average2))/(cont*sigma1*sigma2);
 //    printf("   corr= %f\n",corrTotal);
     return corrTotal;
 }

  // Mon made (13/9/2010) (DOUBLE precision in accumulations)
  // Cross.Corr. inside a cubic-frame (nozero=false) or inside the final volume (nozero=true, default)
  // size = kernel width
  // vol = Target (Final) volume
  // If we have avg1 and sig1 (from: vol), we should provide them!
  float correlation_frame(vlVolume *vol2, vlVolume *vol, vlDim size, bool nozero, bool volcalc, double volavg, double volsig)
  {
      if(vol->voxelCount() != vol2->voxelCount())
      {
    	  fprintf(stderr,"Correlation Error: Dimension mismatch\n");
      }

      double corrTotal=0.0;
      double average1=0.0,average2=0.0,sigma1=0.0,sigma2=0.0;
      int i,j,z;
      int cont=0,null1=0,null2=0;
      float *p_voxel,*p_voxel2;
      float voxel,voxel2;
      int stepx,stepy,stepz;
      int offset;

      vlStep step=vol->stepping();
      stepx=step.x();
      stepy=step.y();
      stepz=step.z();

      vlDim dim=vol->dim();

      float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
      float *container2= (float*)vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

      if(nozero)
      { // Inside final map (>zero voxels)
  		// Map 1 & 2 average & dot product
      	if(volcalc) // then compute BOTH averages & sigmas
      	{
  		    for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          if( voxel > 0.0 ) // "vol" >0 mask
  		          {
  		            voxel2= *(container2+offset); // 2nd vol voxel
  	          	    average2 += voxel2;
  		          	average1 += voxel;
  		          	corrTotal += voxel * voxel2;
  		          	cont++;
  		          }
  		        }
  		    average1 /= cont;
  		    average2 /= cont;

  		    // Map 1 & 2 sigma
  		    for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          if( voxel > 0.0 )
  		          {
  			        voxel2= *(container2+offset); // 2nd vol voxel
  		          	sigma1 += pow(voxel-average1,2);
  		          	sigma2 += pow(voxel2-average2,2);
  		          }
  		        }
  		     sigma1 = sqrt(sigma1/cont);
  		     sigma2 = sqrt(sigma2/cont);
      	}
      	else // compute SINGLE average & sigma (vol2's avg.)
      	{
      		average1 = volavg;
      		sigma1 = volsig;

      		for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          if( voxel > 0.0 ) // "vol" >0 mask
  		          {
  		            voxel2= *(container2+offset); // 2nd vol voxel
  		          	average2 += voxel2;
  		          	corrTotal += voxel * voxel2;
  	          	    cont++;
  		          }
  		        }
  		    average2 /= cont;

  		    // Map 1 & 2 sigma
  		    for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          if( voxel > 0.0 )
  		          {
  			        voxel2= *(container2+offset); // 2nd vol voxel
  		          	sigma2 += pow(voxel2-average2,2);
  		          }
  		        }
  		     sigma2 = sqrt(sigma2/cont);
      	}
      }
      else
      { // Inside CUBIC FRAME

   	// Number of voxels (without taking into account the frame)
       cont = (dim.x() - 2*size.x() ) * (dim.y() - 2*size.y() ) * (dim.z() - 2*size.z() );

       // Map 1 & 2 average & dot product
      	if(volcalc) // then compute BOTH averages & sigmas
      	{
  		    for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          voxel2= *(container2+offset); // 2nd vol voxel
  		          if(voxel != 0)
  		        	  average1 += voxel;
  		          else
  		        	  null1++; // counting null voxels (speeds-up sigma)
  		          if(voxel2 != 0)
  		          {
  		        	  average2 += voxel2;
  		        	  corrTotal += voxel * voxel2;
  		          }
  		          else
 		        	  null2++; // counting null voxels (speeds-up sigma)
  		        }
  		    average1 /= cont;
  		    average2 /= cont;

  		    // Map 1 & 2 sigma
  		    for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          if(voxel != 0)
  		        	  sigma1 += pow(voxel-average1,2);
 	        	  voxel2= *(container2+offset); // 2nd vol voxel
  		          if(voxel2 != 0)
  		        	  sigma2 += pow(voxel2-average2,2);
  		        }
  		     sigma1 += null1 * pow(average1,2);
  		     sigma1 =  sqrt(sigma1/cont);
  		     sigma2 += null2 * pow(average2,2);
  		     sigma2 =  sqrt(sigma2/cont);
      	}
      	else // compute SINGLE average & sigma (vol2's avg.)
      	{
      		average1 = volavg;
      		sigma1 = volsig;

      		for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  		          voxel= *(container+offset); // 1st vol voxel
  		          voxel2= *(container2+offset); // 2nd vol voxel
  		          if(voxel2 != 0)
  		          {
  		        	  average2 += voxel2;
  		        	  corrTotal += voxel * voxel2;
  		          }
  		          else
  		        	  null2++;
  		        }
  		    average2 /= cont;

  		    // Map 1 & 2 sigma
  		    for(i=size.x();i<dim.x()-size.x();i++)
  		      for(j=size.y();j<dim.y()-size.y();j++)
  		        for(z=size.z();z<dim.z()-size.z();z++)
  		        {
  		          offset = (i*stepx)+(j*stepy)+(z*stepz);
  			      voxel2= *(container2+offset); // 2nd vol voxel
  			      if(voxel2 != 0)
  			    	  sigma2 += pow(voxel2-average2,2);

  		        }
 		     sigma2 += null2 * pow(average2,2);
  		     sigma2 = sqrt(sigma2/cont);
      	}
      }

  	// Normalized cross-correlation computation
  //    printf("Msg(correlation_frame): dot= %f   avg1= %f   avg2= %f   sig1= %f   sig2= %f",corrTotal,average1,average2,sigma1,sigma2);
      corrTotal = (corrTotal-(cont*average1*average2))/(cont*sigma1*sigma2);
  //    printf("   corr= %f\n",corrTotal);
      return (float)corrTotal;
  }

// Mon made (3/7/2008)
// Local Cross.Corr. inside a cubic-frame (nozero=false) or inside the Final volume (vol)
// (nozero=true, default)--> Warning, only NON-ZERO voxels will be taken into account!
// size = kernel width
// vol = Final volume
float localcorr(vlVolume **p_corrmap, vlVolume *vol2, vlVolume *vol,int size, bool nozero)
{
    if(vol->voxelCount() != vol2->voxelCount())
    {
    	fprintf(stderr,"Local Correlation Error: Different map sizes\n");
      exit(1);
    }

    float aux,corrTotal=0.0;
    float average1=0,average2=0,sigma1=0,sigma2=0;
    int i,j,z;
    int cont=0;
    float *p_voxel,*p_voxel2;
    float voxel,voxel2,dot;
    int stepx,stepy,stepz;
    int stepx2,stepy2,stepz2,offset;

    vlVolume *corrmap;
    corrmap = *p_corrmap;
	if(corrmap == NULL) // allocates volume memory
		corrmap = new vlVolume( vol );
	corrmap->clear(0); // clears volume memory

    vlStep step=vol->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();

    vlDim dim=vol->dim();
    int border= (size-1)/2; // respect to the kernell size
    float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container2= (float*)vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container_corr= (float*)corrmap->getVoxelVoidPtr(vlPoint3ui(0,0,0));

    if(nozero)
    { // Inside Final map (non-zero voxels)
		// Map 1 & 2 average & dot product
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          if( voxel != 0.0 )
	          {
	            dot = voxel * voxel2;
	            *(container_corr+offset) = dot; // dot-product voxel
	          	corrTotal += dot;
	          	average1 += voxel;
	          	average2 += voxel2;
	          	cont++;
	          }
	        }
	    average1 /= cont;
	    average2 /= cont;

	    // Map 1 & 2 sigma
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          if( voxel != 0.0 )
	          {
	          	sigma1 += pow(voxel-average1,2);
	          	sigma2 += pow(voxel2-average2,2);
	          }
	        }
	     sigma1 = sqrt(sigma1/cont);
	     sigma2 = sqrt(sigma2/cont);
    }
    else
    { // Inside Cubic Frame
		// Map 1 & 2 average & dot product
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          dot = voxel * voxel2;
	          *(container_corr+offset) = dot; // dot-product voxel
	          corrTotal += dot;
	       	  average1 += voxel;
	       	  average2 += voxel2;
	          cont++;
	        }
	    average1 /= cont;
	    average2 /= cont;

	    // Map 1 & 2 sigma
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          sigma1 += pow(voxel-average1,2);
	          sigma2 += pow(voxel2-average2,2);
	        }
	     sigma1 = sqrt(sigma1/cont);
	     sigma2 = sqrt(sigma2/cont);
    }

	// Normalized cross-correlation computation
    corrTotal = (corrTotal-(cont*average1*average2))/((cont)*sigma1*sigma2);
    *p_corrmap = corrmap;
    return corrTotal;
}

// Mon made (15/7/2008)
// Volume Spherical Average
// size = kernel width (to avoid unnecessary computations)
// radius --> [1,(size-1)/2]
void spherical_avg(vlVolume *vol_in, vlVolume *vol_out, int size, int radius)
{
    if(vol_in->voxelCount() != vol_out->voxelCount())
    {
    	fprintf(stderr,"Spherical Average Error: Different map sizes\n");
      exit(1);
    }

    int i,j,z,a,b,c;
    int cont=0;
    float voxel;
    int stepx,stepy,stepz;
    int stepx2,stepy2,stepz2,offset;

    vlStep step=vol_in->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();

    vlDim dim=vol_in->dim();
    int border= (size-1)/2; // respect to the kernell size

    if(radius>border)
    {
      printf("Spherical Average Error: Averaging Radius bigger than Border! (%d > %d)",radius,border);
      exit(1);
    }

    float *container= (float*)vol_in->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	vol_out->clear(0); // clears volume memory
    float *container2= (float*)vol_out->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Map 1 & 2 average & dot product
    for(i=border;i<dim.x()-border;i++)
      for(j=border;j<dim.y()-border;j++)
        for(z=border;z<dim.z()-border;z++)
        {

          voxel=0;
          cont=0;
          for(a=i-radius; a<=i+radius; a++)
          	for(b=j-radius; b<=j+radius; b++)
          	  for(c=z-radius; c<=z+radius; c++)
          	  	if( sqrt( pow((float)(a-i),2) + pow((float)(b-j),2) + pow((float)(c-z),2) ) < radius )
				{  // whether its inside a "radius" sphere
          	  	  voxel+= *(container+a*stepx+b*stepy+c*stepz);
          	  	  cont++;
				}
          voxel /= cont; // computes average
		  *(container2 + (i*stepx)+(j*stepy)+(z*stepz)) = voxel;
        }
}

// Mon made (15/7/2008)
// Weighted Correlation Map Normalization
// p_norm --> [0,1]
void norm_wcorr(vlVolume *vol, float p_norm)
{
	// Weight-Map Normalization
	float max_weight = FOPS::max_mass(vol);
	FOPS::mul( vol, -p_norm / max_weight );
	FOPS::add( vol, 1 );
}

// Mon made (3/7/2008)
// Computes the Weighted Cross-Correlation between two maps given a correlation map
// Weight will be applied to the "vol" map (the 2nd one)
float correlation_weight(vlVolume *vol2, vlVolume *vol, vlVolume *volw, int size, bool nozero)
{
    if( !( vol->voxelCount() == vol2->voxelCount() && vol->voxelCount() == volw->voxelCount() ) )
    	fprintf(stderr,"Correlation Error: Dimension mismatch\n");

    float aux,corrTotal=0.0;
    float average1=0,average2=0,sigma1=0,sigma2=0;
    int i,j,z;
    int cont=0;
    float *p_voxel,*p_voxel2;
    float voxel,voxel2,weight;
    int stepx,stepy,stepz,offset;

    vlStep step=vol->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();
    vlDim dim=vol->dim();
    int border= (size-1)/2; // respect to the kernell size
    float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container2= (float*)vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container_weight = (float*) volw->getVoxelVoidPtr(vlPoint3ui(0,0,0));

    if(nozero)
    { // Inside final map (non-zero voxels)
		// Map 1 & 2 average & dot product
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          weight= *(container_weight+offset); // weight vol voxel
	          voxel *= weight; // weighting "vol"
	          if( voxel != 0.0 )
	          {
	          	average1 += voxel;
	          	average2 += voxel2;
	          	corrTotal += voxel * voxel2;
	          	cont++;
	          }
	        }
	    average1 /= cont;
	    average2 /= cont;

	    // Map 1 & 2 sigma
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          weight= *(container_weight+offset); // weight vol voxel
	          voxel *= weight; // weighting "vol"
	          if( voxel != 0.0 )
	          {
	          	sigma1 += pow(voxel-average1,2);
	          	sigma2 += pow(voxel2-average2,2);
	          }
	        }
	     sigma1 = sqrt(sigma1/cont);
	     sigma2 = sqrt(sigma2/cont);
    }
    else
    { // Inside Cubic Frame
		// Map 1 & 2 average & dot product
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          weight= *(container_weight+offset); // weight vol voxel
	          voxel *= weight; // weighting
	       	  average1 += voxel;
	       	  average2 += voxel2;
	          corrTotal += voxel * voxel2;
	          cont++;
	        }
	    average1 /= cont;
	    average2 /= cont;

	    // Map 1 & 2 sigma
	    for(i=border;i<dim.x()-border;i++)
	      for(j=border;j<dim.y()-border;j++)
	        for(z=border;z<dim.z()-border;z++)
	        {
	          offset = (i*stepx)+(j*stepy)+(z*stepz);
	          voxel= *(container+offset); // 1st vol voxel
	          voxel2= *(container2+offset); // 2nd vol voxel
	          weight= *(container_weight+offset); // weight vol voxel
	          voxel *= weight; // weighting
	          sigma1 += pow(voxel-average1,2);
	          sigma2 += pow(voxel2-average2,2);
	        }
	     sigma1 = sqrt(sigma1/cont);
	     sigma2 = sqrt(sigma2/cont);
    }

	// Normalized cross-correlation computation
    corrTotal = (corrTotal-(cont*average1*average2))/(cont*sigma1*sigma2);
    return corrTotal;
}


//Convolucion de un kernel sobre un mapa
// vlm= Mapa sobre el que se aplica el kernel
// kernel= Kernel a aplicar
// size= Dimension del Kernel
vlVolume* convoluteK(vlVolume *vlm,float *kernel, int size, bool self)
  {
    int i,j,z;
    int kx,ky,kz;
    int cont;
    float total;
    float *voxel,*voxel2;
    int stepx,stepy,stepz;
    int stepx2,stepy2,stepz2;

//    int numVoxels=vlm->voxelCount();
    vlStep step=vlm->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();

    vlDim dim=vlm->dim();
    int border= (size-1)/2;

    float *container= (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));


    vlVolume *n= new vlVolume(vlm,vlDim(border,border,border));
    step=n->stepping();
    stepx2=step.x();
    stepy2=step.y();
    stepz2=step.z();


    float *container2= (float*)n->getVoxelVoidPtr(vlPoint3ui(border,border,border));
    float voxval;

    for(i=0;i<dim.x();i++)
      for(j=0;j<dim.y();j++)
        for(z=0;z<dim.z();z++)
        {
          voxval= *(container+(i*stepx)+(j*stepy)+(z*stepz));
          if (voxval!=0)  {
          voxel2=container2+(i*stepx2)+(j*stepy2)+(z*stepz2);
          total=0;
          cont=0;
          for(kx=-border;kx<=border;kx++)
            for(ky=-border;ky<=border;ky++)
              for(kz=-border;kz<=border;kz++)
              {
                *(voxel2+(kx*stepx2)+(ky*stepy2)+(kz*stepz2))+=(float)kernel[cont]* voxval;;
                cont++;
              }
          }
        }

    if(self)
    {
      delete vlm;
    }
    return n;
  }

// Mon made (29/5/2008)
// Convolucion de un kernel sobre un mapa
// vlm= Mapa sobre el que se aplica el kernel
// kernel= Kernel a aplicar
// size= Dimension del Kernel
// Volume must be big enought to allocate the kernel in the border (pad it before!)
// vlVolume "n" should be already allocated!
vlVolume* convoluteK_nopad(vlVolume *vlm, vlVolume *n, float *kernel, int size, bool self)
  {
    int i,j,z,kx,ky,kz,cont;
    float *voxel,*voxel2;
    int stepx,stepy,stepz;

	// clearing dummy volume
    n->clear(0);

    vlStep step=vlm->stepping();
    stepx=step.x(); // Both maps should have equal "step" (I would rather say "offset"...)
    stepy=step.y();
    stepz=step.z();

    vlDim dim=vlm->dim();
    int border= (size-1)/2;

    float *container= (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container2= (float*)n->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float voxval;

    int xf = dim.x()-border;
    int yf = dim.y()-border;
    int zf = dim.z()-border;

    for(i=border;i<xf;i++)
      for(j=border;j<yf;j++)
        for(z=border;z<zf;z++)
        {
          voxval= *(container+(i*stepx)+(j*stepy)+(z*stepz));
          if (voxval!=0)
          {
              voxel2=container2+(i*stepx)+(j*stepy)+(z*stepz);
        	  cont=0;
        	  for(kx=-border;kx<=border;kx++)
        		  for(ky=-border;ky<=border;ky++)
        			  for(kz=-border;kz<=border;kz++)
        			  {
        				  *(voxel2+(kx*stepx)+(ky*stepy)+(kz*stepz))+=(float)kernel[cont]* voxval;;
        				  cont++;
        			  }
          }
        }

    if(self) // false by default
    {
      delete vlm;
    }
    return n;
  }


#ifdef USE_TBB // Enables Intel's Threading Building Blocks parallel routines

// Convolutes just one "i" slice by a kernel
void ConvSlice(int i, convoluteK_data const *data)
{
	int j,z,kx,ky,kz,cont;
	float voxval;
	float *voxel2;
	int border = ((convoluteK_data *) data)->border;
	int yf = ((convoluteK_data *) data)->j_final;
	int zf = ((convoluteK_data *) data)->z_final;
	int stepx = ((convoluteK_data *) data)->stepx;
	int stepy = ((convoluteK_data *) data)->stepy;
	int stepz = ((convoluteK_data *) data)->stepz;
	float *container = ((convoluteK_data *) data)->container;
	float *container2 = ((convoluteK_data *) data)->container2;
	float *kernel = ((convoluteK_data *) data)->kernel;

	for(j=border;j<yf;j++)
		for(z=border;z<zf;z++)
		{
			voxval= *(container+(i*stepx)+(j*stepy)+(z*stepz));
			if (voxval!=0)
			{
				voxel2=container2+(i*stepx)+(j*stepy)+(z*stepz);
				cont=0;
				for(kx=-border;kx<=border;kx++)
					for(ky=-border;ky<=border;ky++)
						for(kz=-border;kz<=border;kz++)
						{
							*(voxel2+(kx*stepx)+(ky*stepy)+(kz*stepz))+=(float)kernel[cont]* voxval;
							cont++;
						}
			}
		}
}

void ConvSliceX(int i, convoluteK_data const *data)
{
	int x,j,z,kx,ky,kz,cont,end;
	float voxval;
	float *voxel2;
	int chunk = ((convoluteK_data *) data)->chunk;
	int border = ((convoluteK_data *) data)->border;
	int xf = ((convoluteK_data *) data)->i_final;
	int yf = ((convoluteK_data *) data)->j_final;
	int zf = ((convoluteK_data *) data)->z_final;
	int stepx = ((convoluteK_data *) data)->stepx;
	int stepy = ((convoluteK_data *) data)->stepy;
	int stepz = ((convoluteK_data *) data)->stepz;
	float *container = ((convoluteK_data *) data)->container;
	float *container2 = ((convoluteK_data *) data)->container2;
	float *kernel = ((convoluteK_data *) data)->kernel;

	if( (((xf-border)/chunk)-1) == i ) // if last chunk
		end = xf; // the thread working in last chunk will do some extra slices...
	else
		end = chunk*(i+1) + border;

//	fprintf(stderr,"working from slice %d to %d\n",border+chunk*i,end);

	for(x = border + chunk*i; x < end; x++)
	{
		for(j=border;j<yf;j++)
			for(z=border;z<zf;z++)
			{
				voxval= *(container+(x*stepx)+(j*stepy)+(z*stepz));
				if (voxval!=0)
				{
					voxel2=container2+(x*stepx)+(j*stepy)+(z*stepz);
					cont=0;
					for(kx=-border;kx<=border;kx++)
						for(ky=-border;ky<=border;ky++)
							for(kz=-border;kz<=border;kz++)
							{
								*(voxel2+(kx*stepx)+(ky*stepy)+(kz*stepz))+=(float)kernel[cont]* voxval;
								cont++;
							}
				}
			}
	}
}


//#include "tbb/tbb.h"
//using namespace tbb;
class ApplyConvSlice {
	convoluteK_data const *my_a;
public:
	ApplyConvSlice (convoluteK_data const *data) // constructor form "convoluteK_data"?
	{
		my_a = data;
	}
	void operator()( const tbb::blocked_range<int>& r ) const
	{
		convoluteK_data const *data = my_a;
		for( int i=r.begin(); i!=r.end(); ++i )
			ConvSlice(i, data);
	}
};

class ApplyConvSliceX {
	convoluteK_data const *my_a;
public:
	ApplyConvSliceX (convoluteK_data const *data) // constructor from "convoluteK_data"?
	{
		my_a = data;
	}
	void operator()( const tbb::blocked_range<int>& r ) const
	{
		convoluteK_data const *data = my_a;
		for( int i=r.begin(); i!=r.end(); ++i )
			ConvSliceX(i, data);
	}
};


// Mon made (3/6/2012)
// Convolucion de un kernel sobre un mapa using Intel's TBB "parallel_for"
// vlm= Mapa sobre el que se aplica el kernel
// kernel= Kernel a aplicar
// size= Dimension del Kernel
// Volume must be big enought to allocate the kernel in the border (pad it before!)
// vlVolume "n" should be already allocated!
vlVolume* convoluteK_nopad_tbb(vlVolume *vlm, vlVolume *n, float *kernel, int size, int nthreads, int granularity, bool self)
{
	int i,j,z,kx,ky,kz,cont;
	float *voxel,*voxel2;
	int stepx,stepy,stepz;

	// clearing dummy volume
	n->clear(0);

	vlStep step=vlm->stepping();
	stepx=step.x(); // Both maps should have equal "step" (I would rather say "offset"...)
	stepy=step.y();
	stepz=step.z();

	vlDim dim=vlm->dim();
	int border= (size-1)/2;

	float *container= (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container2= (float*)n->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float voxval;

	int xf = dim.x()-border;
	int yf = dim.y()-border;
	int zf = dim.z()-border;

	convoluteK_data data;
	data.border = border;
	data.j_final = yf;
	data.z_final = zf;
	data.stepx = stepx;
	data.stepy = stepy;
	data.stepz = stepz;
	data.container = container;
	data.container2 = container2;
	data.kernel = kernel;
	convoluteK_data *pdata=&data;

//	// One slice per thread...
//	if(granularity > 0)
//		tbb::parallel_for(tbb::blocked_range<int>(border,xf,granularity), ApplyConvSlice(pdata));
//	else
//		tbb::parallel_for(tbb::blocked_range<int>(border,xf), ApplyConvSlice(pdata));

	// One chunk per thread... (better map reliability 100%)
	data.chunk = (xf-border)/nthreads;
	if(granularity > 0) // (always even number of chunks (= nthreads*2)
		tbb::parallel_for(tbb::blocked_range<int>(0,nthreads,granularity), ApplyConvSliceX(pdata)); // do even index slices
	else
		tbb::parallel_for(tbb::blocked_range<int>(0,nthreads), ApplyConvSliceX(pdata)); // do even index slices

	if(self) // false by default
	{
		delete vlm;
	}
	return n;
}

#endif


#ifdef USE_PTHREAD // Enables PThread parallel routines
// Mon made (22/5/2013)
// Multi-thread routine to convolute some cubic region of a map with a given kernel.
void *convoluteK_nopad_thread(void *threadarg)
{
	int border,i,j,z,xi,xf,yi,yf,zi,zf,cont,kx,ky,kz,stepx,stepy,stepz;
	float *container,*container2,*kernel,*voxel2;
    float voxval;

    while(true) // loops for ever...
    {
        pthread_mutex_lock(((convoluteK_data *) threadarg)->p_mutex_begin);
        while( ((convoluteK_data *) threadarg)->begin ) // Checking whether the thread would begin ("re-signal" for "wait" safety)
    		pthread_cond_wait(((convoluteK_data *) threadarg)->p_cond_begin, ((convoluteK_data *) threadarg)->p_mutex_begin);
        ((convoluteK_data *) threadarg)->begin = true; // for next iteration... (enables "beginability")
        pthread_mutex_unlock(((convoluteK_data *) threadarg)->p_mutex_begin);

    	// Loading data from "thread argument"
    	border = ((convoluteK_data *) threadarg)->border;
    	xi = ((convoluteK_data *) threadarg)->i_initial;
    	xf = ((convoluteK_data *) threadarg)->i_final;
    	yi = ((convoluteK_data *) threadarg)->j_initial;
    	yf = ((convoluteK_data *) threadarg)->j_final;
    	zi = ((convoluteK_data *) threadarg)->z_initial;
    	zf = ((convoluteK_data *) threadarg)->z_final;
    	stepx = ((convoluteK_data *) threadarg)->stepx;
    	stepy = ((convoluteK_data *) threadarg)->stepy;
    	stepz = ((convoluteK_data *) threadarg)->stepz;
    	container = ((convoluteK_data *) threadarg)->container;
    	container2 = ((convoluteK_data *) threadarg)->container2;
    	kernel = ((convoluteK_data *) threadarg)->kernel;

//    	for(i = ((convoluteK_data *) threadarg)->i_initial; i < ((convoluteK_data *) threadarg)->i_final; i++)
//    		for(j = ((convoluteK_data *) threadarg)->j_initial; j < ((convoluteK_data *) threadarg)->j_final; j++)
//    			for(z = ((convoluteK_data *) threadarg)->z_initial; z < ((convoluteK_data *) threadarg)->z_final; z++)
//    			{
//    				voxval= *(((convoluteK_data *) threadarg)->container+(i*((convoluteK_data *) threadarg)->stepx)+(j*((convoluteK_data *) threadarg)->stepy)+(z*((convoluteK_data *) threadarg)->stepz));
//    				if (voxval!=0)
//    				{
//    					voxel2 = ((convoluteK_data *) threadarg)->container2+(i*((convoluteK_data *) threadarg)->stepx)+(j*((convoluteK_data *) threadarg)->stepy)+(z*((convoluteK_data *) threadarg)->stepz);
//    					cont=0;
//    					for(kx = -((convoluteK_data *) threadarg)->border; kx <= ((convoluteK_data *) threadarg)->border; kx++)
//    						for(ky = -((convoluteK_data *) threadarg)->border; ky <= ((convoluteK_data *) threadarg)->border; ky++)
//    							for(kz = -((convoluteK_data *) threadarg)->border; kz <= ((convoluteK_data *) threadarg)->border; kz++)
//    							{
//    								*(voxel2+(kx*((convoluteK_data *) threadarg)->stepx)+(ky*((convoluteK_data *) threadarg)->stepy)+(kz*((convoluteK_data *) threadarg)->stepz)) += (float)((convoluteK_data *) threadarg)->kernel[cont]* voxval;
//    								cont++;
//    							}
//    				}
//    			}

    	for(i = xi; i < xf; i++)
    		for(j = yi; j < yf; j++)
    			for(z = zi; z < zf; z++)
    			{
    				voxval= *(container+(i*stepx)+(j*stepy)+(z*stepz));
    				if (voxval!=0)
    				{
    					voxel2 = container2+(i*stepx)+(j*stepy)+(z*stepz);
    					cont=0;
    					for(kx = -border; kx <= border; kx++)
    						for(ky = -border; ky <= border; ky++)
    							for(kz = -border; kz <= border; kz++)
    							{
    								*(voxel2+(kx*stepx)+(ky*stepy)+(kz*stepz)) += (float)kernel[cont]* voxval;
    								cont++;
    							}
    				}
    			}

    	// Here, sending the end signal to main thread...
    	pthread_mutex_lock(((convoluteK_data *) threadarg)->p_mutex_end);
    	(*(((convoluteK_data *) threadarg)->p_nended))++; // increasing "nended" counter... for "wait" safety...
    	pthread_cond_signal(((convoluteK_data *) threadarg)->p_cond_end); // Sending end signal
        pthread_mutex_unlock(((convoluteK_data *) threadarg)->p_mutex_end);
    }
    pthread_exit(NULL);
}

// Mon made (27/5/2013)
// Initialization routine for "convoluteK_nopad_par"
void convoluteK_nopad_par_init(int nthreads, convoluteK_data **p_threads_data, pthread_t **p_threads, float *kernel, int kern_size)
{
	bool debug = false;
	int i,j,rc;
	convoluteK_data *threads_data;
	pthread_t *threads;
	pthread_mutex_t *p_mutex_begin; // begin mutex pointer
	pthread_cond_t *p_cond_begin; // begin condition pointer
	pthread_mutex_t *p_mutex_end; // end mutex pointer
	pthread_cond_t *p_cond_end; // end condition pointer
	int *p_nended; // nended pointer

	// Mutexes and condition creation
	if( !(p_mutex_begin = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)) ) ||
			!(p_mutex_end = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)) ) ||
			!(p_cond_begin = (pthread_cond_t *) malloc(sizeof(pthread_cond_t)) ) ||
			!(p_cond_end = (pthread_cond_t *) malloc(sizeof(pthread_cond_t)) ) ||
			!(p_nended = (int *) malloc(sizeof(int)) ) )
	{
		fprintf(stderr,"Msg(convoluteK_nopad_par_init): I'm sorry, mutex/condition/nended memory allocation failed!\nForcing exit!\n");
		exit(1);
	}

	// Initializing "nended"
	*p_nended = 0;

	// Allocating threads data
	if( !(threads_data = (convoluteK_data *) malloc(sizeof(convoluteK_data) * nthreads)) )
	{
		fprintf(stderr,"Msg(convoluteK_nopad_par_init): I'm sorry, thread memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(threads = (pthread_t *) malloc(sizeof(pthread_t) * nthreads)) )
	{
		fprintf(stderr,"Msg(convoluteK_nopad_par_init): I'm sorry, thread allocation failed!\nForcing exit!\n");
		exit(1);
	}
	*p_threads_data = threads_data;
	*p_threads = threads;

	// Initializing mutexes and conditions (only the first time)
	pthread_mutex_init(p_mutex_begin, NULL); // mutex initialization
	pthread_cond_init(p_cond_begin, NULL); // condition initialization
	pthread_mutex_init(p_mutex_end, NULL); // mutex initialization
	pthread_cond_init(p_cond_end, NULL); // condition initialization

	// kern_size = kern_size * kern_size * kern_size; // Number of voxels of the Convolution Kernel
//	float *kernel2; // Allocate memory for the kernel replica

	for(i = 0; i < nthreads; i++)
	{
		if( !(threads_data[i].kernel = (float *) malloc( sizeof(float) * kern_size * kern_size * kern_size )) )
		{
			fprintf(stderr,"Msg(convoluteK_nopad_par_init): I'm sorry, convolution kernel2 memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		for(j = 0; j < kern_size * kern_size * kern_size; j++)
			(threads_data[i].kernel)[j] = kernel[j]; // Convolution kernel replication (faster memory access?)

		threads_data[i].p_mutex_begin = p_mutex_begin; // passing by reference the begin mutex
		threads_data[i].p_cond_begin = p_cond_begin; // passing by reference the begin condition
		threads_data[i].p_mutex_end = p_mutex_end; // passing by reference the end mutex
		threads_data[i].p_cond_end = p_cond_end; // passing by reference the end condition
		threads_data[i].p_nended = p_nended; // passing by reference the nended variable
		threads_data[i].begin = true; // by default, thread is able to begin...
	}

	// Creating threads
	for(i = 0; i < nthreads; i++)
	{
		if(debug)
			fprintf(stderr,"Creating thread (convoluteK_nopad_par_init): %d\n", i);
		rc = pthread_create(&threads[i], NULL, FOPS::convoluteK_nopad_thread, (void *) &threads_data[i]);
		if (rc)
		{
			fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
}


// Mon made (22/5/2013)
// Parallel Kernel Convolution ("convoluteK_nopad_par_init" should have been called before this one!)
// vlm= Mapa sobre el que se aplica el kernel
// kernel= Kernel a aplicar
// size= Dimension del Kernel
// Volume must be big enought to allocate the kernel in the border (pad it before!)
// vlVolume "n" should be already allocated!
vlVolume *convoluteK_nopad_par(vlVolume *vlm, vlVolume *n, float *kernel, int size, int nthreads, convoluteK_data *threads_data, bool self)
  {
    bool debug = false;
    int i,j,z;
    int kx,ky,kz;
    int cont;
    float *voxel,*voxel2;
    int stepx,stepy,stepz;

	// clearing dummy volume
    n->clear(0);
    vlStep step=vlm->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();

    vlDim dim=vlm->dim();
    int border= (size-1)/2;
    float *container= (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container2= (float*)n->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Loading threads input data
	float nslices = ((float)dim.x()-2*(float)border) / (float)nthreads; // Number of slices per thread (to balance computational burden)
	if(debug)
		fprintf(stderr,"Loading thread data:  tot_slices= %d   nslices= %f  border= %d\n",(dim.x()-2*border),nslices,border);
	for(i = 0; i < nthreads; i++)
	{
		threads_data[i].border = border;
		threads_data[i].i_initial = border + (int)((float)nslices*i); // Initial i-slice index
		threads_data[i].i_final = border + (int)((float)nslices*(i+1)); // Final i-slice index
		threads_data[i].j_initial = border; // Initial j-slice index
		threads_data[i].j_final = dim.y()-border; // Final j-slice index
		threads_data[i].z_initial = border; // Initial z-slice index
		threads_data[i].z_final = dim.z()-border; // Final z-slice index
		threads_data[i].stepx = stepx;
		threads_data[i].stepy = stepy;
		threads_data[i].stepz = stepz;
		threads_data[i].container = container; // Volume pointer
		threads_data[i].container2 = container2; // Dummy volume pointer
		if(debug && i != nthreads-1)
			fprintf(stderr,"Loaded data for thread %d:  i_initial= %d  i_final= %d\n",i,threads_data[i].i_initial,threads_data[i].i_final);
	}
	threads_data[nthreads-1].i_final = dim.x()-border; // this should ensure that all voxels will be convoluted
	if(debug)
		fprintf(stderr,"Loaded data for thread %d:  i_initial= %d  i_final= %d\n",i-1,threads_data[nthreads-1].i_initial,threads_data[nthreads-1].i_final);

	// Awaking all threads!
	if(debug)
		fprintf(stderr, "Awaking all threads! (convoluteK_nopad_par)\n");
    pthread_mutex_lock(threads_data[0].p_mutex_end); // This prevents any end signal miss by controller thread.
    pthread_mutex_lock(threads_data[0].p_mutex_begin);
    for(int n=0; n<nthreads; n++)
    	threads_data[n].begin = false; // "re-signal" beginning...
	pthread_cond_broadcast(threads_data[0].p_cond_begin); // Send beginning  signal...
    pthread_mutex_unlock(threads_data[0].p_mutex_begin);
    // Waiting till all working threads finish.
	while( *(threads_data[0].p_nended) < nthreads )
	{
    	pthread_cond_wait(threads_data[0].p_cond_end, threads_data[0].p_mutex_end); // gets END signal from threads
    	if(debug)
    		fprintf(stderr,"Some thread has finished:  nended= %d\n",*(threads_data[0].p_nended));
	}
    pthread_mutex_unlock(threads_data[0].p_mutex_end);
    (*(threads_data[0].p_nended)) = 0; // reset "ended treads" counter...

    if(self)
    {
      delete vlm;
    }
    return n;
  }

#endif

void compute_kernel_Gaussian(float **kernel,int *dim_vox,float unit_vox,float resolution,float sigma_factor)
{
  int exth,nvox;
  float bvalue,cvalue;
  float mscale;
  int i,j,k;
  int dim_vox2;
  float dsqu;
  //
  // convention res_situs=sqrt(ln2)xres_eman  0.83255461
  //            res_eman=1.2011224*res_situs
  //
  float sig=  0.83255461*resolution/2.0;
  float kmsd=sig*sig/(unit_vox*unit_vox);
  float sigma1d=sqrt(kmsd/3.0);
  //float sigma1d=resolution/(2.0*(float)vol->units().x()*sqrt(3.0));

  exth=(int)ceil(sigma_factor*sigma1d);
  *dim_vox=2*exth-1;
  dim_vox2=*dim_vox * *dim_vox;
  nvox=*dim_vox * dim_vox2;
  *kernel= (float*)malloc(sizeof(float)*nvox);
  memset(*kernel,0,sizeof(float)*nvox);

  bvalue=-1/(2.0*sigma1d*sigma1d);
  cvalue= sigma_factor*sigma_factor*sigma1d*sigma1d;

  mscale=0;
  for(i=0;i<*dim_vox;i++)
    for(j=0;j<*dim_vox;j++)
      for(k=0;k<*dim_vox;k++)
      {
        dsqu=(i-exth+1)*(i-exth+1)+
              (j-exth+1)*(j-exth+1)+
              (k-exth+1)*(k-exth+1);
        if(dsqu <= cvalue)
          (*kernel)[i*dim_vox2 + j*(*dim_vox)+ k]= exp(dsqu*bvalue);
        mscale+=(*kernel)[i*dim_vox2 + j*(*dim_vox)+ k];
      }
  for(i=0;i<nvox;i++)
    (*kernel)[i]/=mscale;






}

void compute_kernel_Lap_Gaussian(float **kernel,int *dim_vox,float unit_vox,float resolution,float sigma_factor)
{
  int exth,nvox;
  float bvalue,cvalue;
  float mscale;
  int i,j,k,kx,ky,kz;
  int dim_vox2;
  float dsqu, *gauss;

  double lap[3][3][3] = {
    { {0.0, 0.0, 0.0}, {0.0, 1. / 12.0, 0.0}, {0.0, 0.0, 0.0} },
    { {0.0, 1.0 / 12.0, 0.0}, {1.0 / 12.0, -6.0 / 12.0, 1.0 / 12.0}, {0.0, 1.0 / 12.0, 0.0} },
    { {0.0, 0.0, 0.0}, {0.0, 1.0 / 12.0, 0.0}, {0.0, 0.0, 0.0} }
  };
  //
  // convention res_situs=sqrt(ln2)xres_eman  0.83255461
  //            res_eman=1.2011224*res_situs
  //
  float sig=  0.83255461*resolution/2.0;
  float kmsd=sig*sig/(unit_vox*unit_vox);
  float sigma1d=sqrt(kmsd/3.0);
  //float sigma1d=resolution/(2.0*(float)vol->units().x()*sqrt(3.0));

  exth=(int)ceil(sigma_factor*sigma1d);

  // gauss
  int dim_voxg,dim_vox2g, nvoxg;
  dim_voxg=2*exth-1;
  dim_vox2g=dim_voxg * dim_voxg;
  nvoxg=dim_voxg * dim_vox2g;
  gauss= (float*)malloc(sizeof(float)*nvoxg);
  memset(gauss,0,sizeof(float)*nvoxg);

  // lap of gauss
  *dim_vox=2*exth+1;
  dim_vox2=*dim_vox * *dim_vox;
  nvox=*dim_vox * dim_vox2;
  *kernel= (float*)malloc(sizeof(float)*nvox);
  memset(*kernel,0,sizeof(float)*nvox);



  bvalue=-1/(2.0*sigma1d*sigma1d);
  cvalue= sigma_factor*sigma_factor*sigma1d*sigma1d;

  mscale=0;
  for(i=0;i<dim_voxg;i++)
    for(j=0;j<dim_voxg;j++)
      for(k=0;k<dim_voxg;k++)
      {
        dsqu=(i-exth+1)*(i-exth+1)+
             (j-exth+1)*(j-exth+1)+
             (k-exth+1)*(k-exth+1);
        if(dsqu <= cvalue)
          gauss[i*dim_vox2g + j*(dim_voxg)+ k]= exp(dsqu*bvalue);
        mscale+=gauss[i*dim_vox2g + j*(dim_voxg)+ k];
      }

  for(i=0;i<nvoxg;i++)
	  gauss[i]/=mscale;

  float voxval;
  for(i=1;i<*dim_vox-1;i++)
      for(j=1;j<*dim_vox-1;j++)
        for(k=1;k<*dim_vox-1;k++)
        {
        	voxval= gauss[(i-1)*dim_vox2g + (j-1)*(dim_voxg)+ (k-1)];
        	//(*kernel)[(i)*dim_vox2 + (j)*(*dim_vox)+ k]=voxval;
        	for(kx=-1;kx<=1;kx++)
        	   for(ky=-1;ky<=1;ky++)
        	         for(kz=-1;kz<=1;kz++)  {
        	        	 (*kernel)[(i+kx)*dim_vox2 + (j+ky)*(*dim_vox)+ (k+kz)]
								   +=lap[kx+1][ky+1][kz+1]* voxval;;
        	         }



        }
  mscale=0;
  for(i=0;i<nvox;i++)
	  mscale+=(*kernel)[i];

  for(i=0;i<nvox;i++)
	  (*kernel)[i]/=mscale;

}

// Mon made (5/2/2009)
// Dotproduct inside a cubic-frame (avoiding padding...)
// size = kernel size (avoids non-valid computations)
// vol = Target (Reference) volume
float dotprod_frame(vlVolume *vol2, vlVolume *vol, vlDim size)
{
	if(vol->voxelCount() != vol2->voxelCount())
	{
		fprintf(stderr,"Correlation Error: Dimension mismatch\n");
		exit(1);
	}

	float aux,corrTotal=0.0;
	float average1=0,average2=0,sigma1=0,sigma2=0;
	int i,j,z;
	int cont=0;
	float *p_voxel,*p_voxel2;
	float voxel,voxel2;
	int stepx,stepy,stepz;
	int offset;

	vlStep step=vol->stepping();
	stepx=step.x();
	stepy=step.y();
	stepz=step.z();

	vlDim dim=vol->dim();
	//    int border= (size-1)/2; // respect to the kernell size

	float *container= (float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container2= (float*)vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Inside CUBIC FRAME
	for(i=size.x();i<dim.x()-size.x();i++)
		for(j=size.y();j<dim.y()-size.y();j++)
			for(z=size.z();z<dim.z()-size.z();z++)
			{
				offset = (i*stepx)+(j*stepy)+(z*stepz);
				voxel= *(container+offset); // 1st vol voxel
				voxel2= *(container2+offset); // 2nd vol voxel
				corrTotal += voxel * voxel2;
				//    				cont++;
			}
//    printf("   corr= %f\n",corrTotal);

//	return 1e5/corrTotal;
	return 1e6-corrTotal;
}


/**
* Returns a mask of indices (array of integers) for fast cross-correlation computations. Mon made (24/8/2018)
* IMPORTANT: A negative integer indicates array end!
* vol = Input volume (map)
* thr = Threshold to consider "in-mask" voxels (values above "thr" will be indexed)
*/
int *indices_mask(vlVolume *vol, float thr)
{
	// Pointer to access voxel values
	float *container = (float*) vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Allocate indices array memory
	int *indices = (int *) malloc( sizeof(int) * vol->voxelCount()); // allocate the maximum

	// Screen all voxels
	int n = 0;
	for(int i = 0; i < vol->voxelCount(); i++)
	{
		if( (*( container + i )) > thr ) // if voxel value greater than "thr"
		{
			indices[ n ] = i;
			n++; // Count the total number of voxels avove "thr"
		}
	}

	indices[ n ] = -1; // A negative integer indicates the end of the indices array (the mask)
	n++;

	indices = (int *) realloc(indices, sizeof(int) * n); // Use just the minimum required memory

	return indices;
}

/**
* Returns a mask of indices (array of integers) with two maps. Values above or below the respective thresholds can be selected.
* (Mon made 27/8/2018)
* IMPORTANT: A negative integer indicates array end!
* vol   = Input mask volume (map)
* thr   = Threshold to consider "in-mask" voxels
* over  = If TRUE, values ABOVE "thr" will be indexed, otherwise the values BELOW "thr" will be indexed.
* vol2  = Input mask volume 2 (map)
* thr2  = Threshold 2 to consider "in-mask" voxels
* over2 = If TRUE, values ABOVE "thr2" will be indexed, otherwise the values BELOW "thr2" will be indexed.
*/
int *indices_mask(vlVolume *vol, float thr, bool over, vlVolume *vol2, float thr2, bool over2)
{
	if( vol->voxelCount() != vol2->voxelCount() )
	{
		fprintf(stderr,"indices_mask> Error: Voxel Count mismatch! vol= %d and vol2= %d voxels. Forcing exit!\n", vol->voxelCount(), vol2->voxelCount());
		exit(1);
	}

	// Pointer to access voxel values
	float *container = (float*) vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container2 = (float*) vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Allocate indices array memory
	int *indices = (int *) malloc( sizeof(int) * vol->voxelCount()); // allocate the maximum

	// 4 different combinations of "over" and "over2" are possible:
	int n = 0;
	if(over && over2)
		for(int i = 0; i < vol->voxelCount(); i++) // Screen all voxels
		{
			if( (*( container + i )) > thr && (*( container2 + i )) > thr2 ) // if both voxel values are greater than the respective thresholds
			{
				indices[ n ] = i;
				n++; // Count the total number of voxels avove "thr"
			}
		}
	if(over && !over2)
		for(int i = 0; i < vol->voxelCount(); i++) // Screen all voxels
		{
			if( (*( container + i )) > thr && (*( container2 + i )) < thr2 ) // if vol voxel value is greater than "thr" and vol2 voxel value is smaller than "thr2"
			{
				indices[ n ] = i;
				n++; // Count the total number of voxels avove "thr"
			}
		}
	if(!over && !over2)
		for(int i = 0; i < vol->voxelCount(); i++) // Screen all voxels
		{
			if( (*( container + i )) < thr && (*( container2 + i )) < thr2 ) // if both voxel values are smaller than the respective thresholds
			{
				indices[ n ] = i;
				n++; // Count the total number of voxels avove "thr"
			}
		}
	if(!over && over2)
		for(int i = 0; i < vol->voxelCount(); i++) // Screen all voxels
		{
			if( (*( container + i )) < thr && (*( container2 + i )) > thr2 ) // if vol voxel value is smaller than "thr" and vol2 voxel value is greater than "thr2"
			{
				indices[ n ] = i;
				n++; // Count the total number of voxels avove "thr"
			}
		}

	indices[ n ] = -1; // A negative integer indicates the end of the indices array (the mask)
	n++;

	indices = (int *) realloc(indices, sizeof(int) * n); // Use just the minimum required memory

	return indices;
}

/**
* Create a masked output map from the input "vol" map and an integers-array mask
* (Mon made 27/8/2018)
* IMPORTANT: A negative integer indicates array end!
* vol = Input volume (map)
* mask = Mask of indices (array of integers), e.g. for fast cross-correlation. WARNING: A negative integer must indicate the array end!
*/
vlVolume *apply_mask(vlVolume *vol, int *mask)
{
	// Copy input "vol" map
	vlVolume *out = new vlVolume( vol );
	out->clear(0.0); // MON: this is not efficient... Is there any "empty map constructor"???

	// Pointer to access voxel values
	float *container_vol = (float*) vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container_out = (float*) out->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Screen all voxels
	int n = 0;
	int x = 0;
	for(int i = 0; i < out->voxelCount(); i++) // Screen the whole map
	{
		if( (x = mask[n]) >= 0 ) // Only copy "in-mask" values from "vol" to "out"
		{
			*( container_out + x ) = *( container_vol + x );
			n++; // counter to screen mask
		}
		else
			break; // exit (mask end reached)
	}

	return out; // Dump output map
}


/**
* Returns the dot product between two maps (faster than cross-correlation). Mon made (24/8/2018)
* IMPORTANT: A negative integer indicates array end!
* vol  = Input volume 1
* vol2 = Input volume 2
* mask = Mask of indices (array of integers) for fast cross-correlation. WARNING: A negative integer must indicate the array end!
*/
float dotprod_mask(vlVolume *vol, vlVolume *vol2, int *mask)
{
	if( vol->voxelCount() != vol2->voxelCount() )
	{
		fprintf(stderr,"dotprod_mask> Error: Voxel Count mismatch! vol= %d and vol2= %d voxels. Forcing exit!\n", vol->voxelCount(), vol2->voxelCount());
		exit(1);
	}

	// Pointer to access voxel values
	float *container  = (float*) vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container2 = (float*) vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Screen all voxels
	int n = 0;
	int x;
	double dot = 0.0;
	for(unsigned int i = 0; i < vol->voxelCount(); i++)
	{
		x = mask[n]; // get index
		if( x >= 0 ) // only non-negative indices are valid
		{
			dot += *( container + x ) * ( *( container2 + x ) ); // Accumulate dot product
			n++; // Count the total number of voxels
			// fprintf(stderr,"dot= %f   x= %d   n= %d   i= %d   vol2= %f\n",dot,x,n,i,*( container2 + x ));
		}
		else
			break; // break if negative (mask end reached)
	}
	dot /= n; // Average dot

	return dot;
}


//	score = 1 - FOPS::correlation_frame(volT,volR,pad,corr_nozero,false,avgR,sigR);

/**
* Returns the normalized cross-correlation between two maps. Mon made (30/8/2018)
* IMPORTANT: A negative integer indicates array end!
* vol  = Input volume 1
* vol2 = Input volume 2
* mask = Mask of indices (array of integers) for fast cross-correlation. WARNING: A negative integer must indicate the array end!
*/
float correlation_mask(vlVolume *vol, vlVolume *vol2, int *mask)
{
	float value, value2; // Value of the voxels
	int nvox; // Total number of voxels
	int n = 0; // Index of Mask array
	int x; // Current map index
	double dot = 0.0; // Dot product
	double avg = 0.0, avg2 = 0.0; // Map averages
	float sig = 0.0, sig2 = 0.0; // Map sigmas

	// Pointer to access voxel values
	float *container  = (float*) vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float *container2 = (float*) vol2->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Some initial checking
	if( vol->voxelCount() != vol2->voxelCount() )
	{
		fprintf(stderr,"correlation_mask> Error: Voxel Count mismatch! vol= %d and vol2= %d voxels. Forcing exit!\n", vol->voxelCount(), vol2->voxelCount());
		exit(1);
	}

	// Get total number of voxels
	nvox = vol->voxelCount();

	// Compute the accumulated dot-product and the averages
	for(int i = 0; i < nvox; i++) // Screen all voxels
	{
		x = mask[n]; // get index
		if( x >= 0 ) // only non-negative indices are valid
		{
			// Get values from maps
			value  = *( container  + x );
			value2 = *( container2 + x );

			// Accumulate dot product
			dot  += value * value2;

			// Accumulate averages
			avg  += value;
			avg2 += value2;

			n++; // Count the total number of voxels
			// fprintf(stderr,"dot= %f   x= %d   n= %d   i= %d   vol2= %f\n",dot,x,n,i,*( container2 + x ));
		}
		else
			break; // break if negative (mask end reached)
	}
	avg  /= nvox; // Average dot
	avg2 /= nvox; // Average dot

	// Compute the sigmas
	n = 0;
	for(int i = 0; i < nvox; i++) // Screen all voxels
	{
		x = mask[n]; // get index
		if( x >= 0 ) // only non-negative indices are valid
		{
			// Accumulate de difference of squares
		   	sig  += powf(*( container  + x ) - avg,  2);
		   	sig2 += powf(*( container2 + x ) - avg2, 2);

		   	n++; // Count the total number of voxels
		}
		else
			break; // break if negative (mask end reached)
	}
    sig  = sqrt(sig / nvox);
    sig2 = sqrt(sig2 / nvox);


 	// Normalized cross-correlation computation
	dot = ( dot - (nvox * avg * avg2) ) / (nvox * sig * sig2);
    return dot;
}



}
