/***************************************************************************
                          fmod.cpp  -  description
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

#define swapF(a,b,t) ((t)=(a),(a)=(b),(b)=(t))


#include "floatops.h"


 namespace FOPS
 {

 //***************MODIFICACION DEL VOLUMEN********************************************
  void threshold(vlVolume *vol, float limit)
  {
    vlVolIter<float, vlLayout::Linear>  iter(vol);
    iter.begin();
    do{
      if(iter.get()<limit)
        iter.set(0.0);
    }while(iter.next());
  }

	 // Quick-sort recursive routine to sort everything in between `low' <-> `high'
	 //  arr   --> Input array of values (float) [v0,v1,v2,...vN-1]
	 //  index --> Input array of indices (int) [0,1,2,...,N-1]
	 //  low   --> Low pivot point (for recursion)
	 //  high  --> High pivot point (for recursion)
 void quicksort(float *arr, int *index, int low, int high)
	 {
	 	int y; // dummy variable for swap ints
	 	int i = low; // left
	 	int j = high; // right
	 	float z = arr[ index[(low + high) / 2] ]; // pivot

	 	do {
	 		while(arr[index[i]] < z) // find member above
	 			i++;

	 		while(arr[index[j]] > z) // // find member below
	 			j--;

	 		if(i <= j) // Swap two elements
	 		{
	            swapF(index[i],index[j],y);
	 			i++;
	 			j--;
	 		}
	 	} while(i <= j);

	 	// recurse
	 	if(low < j)
	 		quicksort(arr, index, low, j);

	 	if(i < high)
	 		quicksort(arr, index, i, high);
	 };

 void norm_threshold(vlVolume *vol, float limit, float dens)
  {

    vlVolIter<float, vlLayout::Linear>  iter(vol);
    double sum, sum0;

    float *arrayD,swaP, value;
    int  dim,i;

    // Allocate memory for the SORTed indices array for RMSDs (sortr)
    dim=vol->dim().x()*vol->dim().y()*vol->dim().z();
	arrayD=(float*)malloc(sizeof(float)*dim);
    int *sorted = (int *) malloc(sizeof(int) * dim);

    double Vfactor=vol->units().x()*vol->units().y()*vol->units().z();

    printf("voxel %f %f \n", Vfactor, vol->units().x());


    iter.begin();
    dim=0;
    do{
      value=iter.get();
      if (value>1E-15) {
      arrayD[dim]=value/Vfactor; dim++;
      }
    } while(iter.next());


    sum0=0;
     for(i=0;i<dim ;i++) {
        sum0+=arrayD[i];
        sorted[i] = i; // Load with the sequential indices
     }
    printf("dim %d sum0 %f dens %f   %f %f %f \n",dim, sum0, dens, arrayD[0],arrayD[1],arrayD[2] );

    if (dens>sum) dens=sum;

    quicksort(arrayD, sorted, 0, dim-1); // Quick sort algorithm

    sum=0;
    for(i=0;i<dim ;i++) {
    	sum+=arrayD[sorted[dim-1-i]];
       // printf("%d %f %f\n", i, sum, arrayD[sorted[dim-1-i]]);
    	if (sum>=dens) break;
    }

    printf("%d %f %f %f\n", i, sum, arrayD[sorted[dim-1-i]], arrayD[sorted[dim-1-i]]*Vfactor );

    // normalize
    // mul(vol,limit/arrayD[sorted[dim-1-i]]);

    free(arrayD);

  }


	void threshold_zero(vlVolume *vol, float limit)
  {
    float limit_neg=limit*-1.0;
    vlVolIter<float, vlLayout::Linear>  iter(vol);
    iter.begin();
    do{
      if(iter.get()<limit && iter.get()>limit_neg)
        iter.set(0.0);
    }while(iter.next());
 }

	void thresholdDown(vlVolume *vol, float limit)
  {
    vlVolIter<float, vlLayout::Linear>  iter(vol);
    iter.begin();
    do{
      if(iter.get()<limit)
        iter.set(limit);
    }while(iter.next());
  }

  void thresholdUp(vlVolume *vol, float limit)
  {
    vlVolIter<float, vlLayout::Linear>  iter(vol);
    iter.begin();
    do{
      if(iter.get()>limit)
        iter.set(limit);
    }while(iter.next());
  }

  void normalize(vlVolume *vol, float factor)
  {
    mul(vol,1.0/factor);
  }



  void changeDensity(vlVolume* vol,float max, float min)
  {
    vlVolIter<float, vlLayout::Linear> iter(vol);
    float value,oldmax,oldmin;

    oldmin=min_mass(vol);
    oldmax=max_mass(vol);

    iter.begin();
    do{
      value=iter.get();
      value= min+ ((max-min)/(oldmax-oldmin)) * (value-oldmin);
      iter.set(value);
    }while(iter.next());
  }

  void changeDensity2(vlVolume* vol,float mean, float sigm)
  {
    vlVolIter<float, vlLayout::Linear> iter(vol);
    float value,oldmean,oldsigma;

    oldmean=calc_average(vol);
    oldsigma=sigma(vol);

    iter.begin();
    do{
      value=iter.get();
      value= ((value-oldmean)*(sigm/oldsigma)) +mean;
      iter.set(value);
    }while(iter.next());
  }

void surface(vlVolume* vol, float limit, std::vector<vlPoint3ui> & boundaryVoxels)
{

  threshold(vol,limit);


  vlNeighborhood neigh;
  bool borde;


  neigh.add(vlOffset(0,0,-1));
  neigh.add(vlOffset(0,0,1));
  neigh.add(vlOffset(0,1,0));
  neigh.add(vlOffset(0,-1,0));
  neigh.add(vlOffset(-1,0,0));
  neigh.add(vlOffset(1,0,0));

  vlVolIter<float, vlLayout::Linear>  iter(vol);
  iter.setNeighborhood(neigh);
  iter.begin();
  do{
    if(iter.get()!=0.0)
    {

      //std::cout<<"hola1"<<std::endl;
      iter.firstNeighbor();
      borde=false;
      do
      {
        //std::cout<<iter.getNeighbor()<<std::endl;
        if(iter.getNeighbor()==0.0)
        {
          //std::cout<<"hola3"<<std::endl;

          borde=true;
          boundaryVoxels.push_back(iter.pos());
        }
      }while(iter.nextNeighbor() && !borde);


    }
  }while(iter.next());

}


void erode(vlVolume* vol, const std::vector<vlPoint3ui> & boundaryVoxels
    ,const uint32 iter, std::vector<vlPoint3ui> & newBoundaryVoxels)
{
  std::vector<vlPoint3ui> currBoundaryVoxels(boundaryVoxels);
  std::vector<vlPoint3ui> nextLayerBoundaryVoxels;

  vlPoint3ui cVoxel; // current voxel
  uint32 growCount(iter);

  vlVolIter<float, vlLayout::Linear> cVoxIter(vol); // iterator at the current voxel

  vlTriple<int32> delta;
  bool first=true;

  // indicate starting
  //std::cout << "-" << std::flush;

  vlNeighborhood neigh;
  neigh.add(vlOffset(0,0,-1));
  neigh.add(vlOffset(0,0,1));
  neigh.add(vlOffset(0,-1,0));
  neigh.add(vlOffset(0,1,0));
  neigh.add(vlOffset(-1,0,0));
  neigh.add(vlOffset(1,0,0));

  cVoxIter.setNeighborhood(neigh);

  // keep growing the specified number of times
  while(growCount>0) {

  if(first)
  {
    while(!currBoundaryVoxels.empty()) {
      cVoxel = currBoundaryVoxels.back();
      cVoxIter.moveTo(cVoxel);
      if(cVoxIter.get()!=0.0)
          cVoxIter.set(0.0);
      currBoundaryVoxels.pop_back();
      nextLayerBoundaryVoxels.push_back(cVoxIter.pos());
    }
    currBoundaryVoxels = nextLayerBoundaryVoxels;
    first=false;
    --growCount;
  }
  else
  {

    // loop over each voxel in the array
    while(!currBoundaryVoxels.empty()) {
      // get the first boundary voxel from the back (since stl vectors only support
      // popping from the back end)
      cVoxel = currBoundaryVoxels.back();
      cVoxIter.moveTo(cVoxel);
      cVoxIter.firstNeighbor();

      do{

        if ( cVoxIter.getNeighbor() != 0.0 ) {
          cVoxIter.setNeighbor(0.0);
          // add the voxel to the next level of boundary voxels
          nextLayerBoundaryVoxels.push_back(cVoxIter.pos()+cVoxIter.getNeighborOffset());
        }
      }while(cVoxIter.nextNeighbor());

      // remove the current voxel from the array
      currBoundaryVoxels.pop_back();
      } // end while(currBoundaryVoxel.size())

      // indicate growing
      //std::cout << "*" << std::flush;

      // decrement the growCount
      --growCount;

      // now copy the new layer voxels back to boundaryVoxels for the next loop
      // do this only when growCount is non-zero..
      if(growCount>0){
        // copy the voxel arrays
        currBoundaryVoxels = nextLayerBoundaryVoxels;

      // empty the next layer voxel array
      nextLayerBoundaryVoxels.clear();
      } else {
      newBoundaryVoxels.clear();
      newBoundaryVoxels = nextLayerBoundaryVoxels;
    }

 } // end-while(growCount)
 }
 }



 void dilate(vlVolume* vol, const std::vector<vlPoint3ui> & boundaryVoxels
      ,const uint32 iter, std::vector<vlPoint3ui> & newBoundaryVoxels, float newTag)
  {
    std::vector<vlPoint3ui> currBoundaryVoxels(boundaryVoxels);
    std::vector<vlPoint3ui> nextLayerBoundaryVoxels;

    vlPoint3ui cVoxel; // current voxel
    uint32 growCount(iter);
    vlVolIter<float, vlLayout::Linear> cVoxIter(vol); // iterator at the current voxel

    bool valid;
    // indicate starting
    //std::cout << "-" << std::flush;

    vlNeighborhood neigh;
    neigh.add(vlOffset(0,0,-1));
    neigh.add(vlOffset(0,0,1));
    neigh.add(vlOffset(0,-1,0));
    neigh.add(vlOffset(0,1,0));
    neigh.add(vlOffset(-1,0,0));
    neigh.add(vlOffset(1,0,0));

    cVoxIter.setNeighborhood(neigh);


    // keep growing the specified number of times
    while(growCount>0) {
    // loop over each voxel in the array
    while(!currBoundaryVoxels.empty()) {
      // get the first boundary voxel from the back (since stl vectors only support
      // popping from the back end)
      cVoxel = currBoundaryVoxels.back();
      cVoxIter.moveTo(cVoxel);

      cVoxIter.firstNeighbor();
      do{

         if ( cVoxIter.getNeighbor() == 0.0 ) {
           valid=cVoxIter.setNeighbor(newTag);
           // add the voxel to the next level of boundary voxels
           if(valid)
           {
            nextLayerBoundaryVoxels.push_back(cVoxIter.pos()+cVoxIter.getNeighborOffset());
            //std::cout<<"voxel: "<< cVoxIter.pos()+cVoxIter.getNeighborOffset() <<std::endl;
           }
         }
      }while(cVoxIter.nextNeighbor());

      // remove the current voxel from the array
      currBoundaryVoxels.pop_back();
    } // end while(currBoundaryVoxel.size())

    // indicate growing
    //std::cout << "*" << std::flush;

    // decrement the growCount
    --growCount;

    // now copy the new layer voxels back to boundaryVoxels for the next loop
    // do this only when growCount is non-zero..
    if(growCount>0){
      // copy the voxel arrays
      currBoundaryVoxels = nextLayerBoundaryVoxels;

      // empty the next layer voxel array
      nextLayerBoundaryVoxels.clear();
    } else {
      newBoundaryVoxels.clear();
      newBoundaryVoxels = nextLayerBoundaryVoxels;
      nextLayerBoundaryVoxels.clear();
    }

   } // end-while(growCount)


  }
//***************FIN MODIFICACION DEL VOLUMEN********************************************



 }
