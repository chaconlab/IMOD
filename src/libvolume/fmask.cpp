#include "floatops.h"
#include "vlipolator_trilinear.h"

 namespace FOPS
 {


 // ***********OPERACIONES DE CREACION DE MASCARA********************

  /*Creacion de mascara de un volumen
    1: voxels con densidad mayor que cutoff
    0: voxels con densidad menor que cutoff
  */
 vlVolume* createMask(vlVolume *vol, float cutoff)
 {
    float value;

    vlVolume *mask= new vlVolume(vol,vlDim(0,0,0));
    vlVolIter<float, vlLayout::Linear> iter(vol);
    vlVolIter<float, vlLayout::Linear> iter2(mask);

    iter.begin();
    iter2.begin();

    while(!iter.end())
    {
      value=iter.get();
      if(value>=cutoff)
        iter2.set(1.0);
      else
        iter2.set(0.0);

      iter.next();
      iter2.next();
    }
    return mask;
 }

  int count_points(vlVolume *vol, float cutoff)
 {
    float value;
    int count;

    vlVolIter<float, vlLayout::Linear> iter(vol);

    iter.begin();
    count=0;
    while(!iter.end())
    {
      value=iter.get();
      if(value>=cutoff)
        count++;
      iter.next();
    }
    return count;
 }


bool inMask(vlVolume *mask,const vlPoint3i & position)
  {
    float aux;
    vlPoint3ui position2;


    if(position.x()<0 || position.x()>=mask->dim().x())
      return false;
    if(position.y()<0 || position.y()>=mask->dim().y())
      return false;
    if(position.z()<0 || position.z()>=mask->dim().z())
      return false;

    position2.x((unsigned int)position.x());
    position2.y((unsigned int)position.y());
    position2.z((unsigned int)position.z());
    mask->getVoxel(position2, aux);

    if(aux==1.0)
      return true;
    else
      return false;
  }


 void eatMask(vlVolume *mask, uint32  dist)
  {
     uint32  layers;
     layers=dist;

     std::vector<vlPoint3ui> boundaryVoxels,newBoundaryVoxels;

     surface(mask, 0.5, boundaryVoxels);

     erode(mask, boundaryVoxels,(uint32) layers, newBoundaryVoxels);
  }

  void dilateMask(vlVolume *mask, uint32  dist)
   {
      uint32  layers;
      layers=dist;

      std::vector<vlPoint3ui> boundaryVoxels,newBoundaryVoxels;

      surface(mask, 0.5, boundaryVoxels);

      dilate(mask, boundaryVoxels,(uint32) layers, newBoundaryVoxels,1.0);

   }


  void beatMask(vlVolume *mask, uint32 layers)
  {
   std::vector<vlPoint3ui> boundaryVoxels,newBoundaryVoxels,newBoundaryVoxels2;

   surface(mask, 0.5, boundaryVoxels);

   dilate(mask,boundaryVoxels, layers, newBoundaryVoxels,1.0);
   boundaryVoxels.clear();
   surface(mask, 0.5, boundaryVoxels);
   erode(mask, boundaryVoxels, layers, newBoundaryVoxels2);
  }

// Fill holes of size fill

void fill_mask(vlVolume *mask, int fill)
  {
    int i,j,z,k;
    float value;
    int stepx,stepy,stepz,last_one;
    bool hole,ones;

    fill+=2; // next pos and the last one

    vlDim dim=mask->dim();
    float *container= (float*)mask->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    vlStep step=mask->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();


    for(i=0;i<dim.x();i++)
      for(j=0;j<dim.y();j++) {
        for(z=0;z<dim.z()-fill;z++)
        {
          value= *(container+(i*stepx)+(j*stepy)+(z*stepz));
          last_one=0;
          if(value==1) {  // check is fill
            hole=false; ones=false;
            for(k=1;k<fill;k++)  // check if is holes nearby
             if ( *(container+(i*stepx)+(j*stepy)+((z+k)*stepz))==0)
             { hole=true;}
             else { if (hole) ones=true; last_one=k;}
            if ((hole)&&(ones))
            for(k=1;k<last_one;k++)
               *(container+(i*stepx)+(j*stepy)+((z+k)*stepz))=1;
            //z+=last_one-1;
          }
        }
        //printf("\ne");
        //for(z=0;z<dim.z();z++)
        //printf("%1.0f",*(container+(i*stepx)+(j*stepy)+(z*stepz)));
  }


  for(z=0;z<dim.z();z++)
   for(j=0;j<dim.y();j++) {
     for(i=0;i<dim.x()-fill;i++)
     {
       value= *(container+(i*stepx)+(j*stepy)+(z*stepz));
       last_one=0;
       if(value==1) {  // check is fill
         hole=false; ones=false;
         for(k=1;k<fill;k++)  // check if is holes nearby
          if ( *(container+((i+k)*stepx)+(j*stepy)+(stepz))==0)
          { hole=true;}
          else { if (hole) ones=true; last_one=k;}
         if ((hole)&&(ones))
         for(k=1;k<last_one;k++)
            *(container+((i+k)*stepx)+(j*stepy)+(stepz))=1;
         //i+=last_one-1;
       }
     }
     //printf("\ne");
     //for(z=0;z<dim.z();z++)
     //printf("%1.0f",*(container+(i*stepx)+(j*stepy)+(z*stepz)));
}




  }






  void invert(vlVolume *mask)
  {
    mul(mask,-1.0);
    add(mask,1.0);
  }

  void fillcenter_onmask(vlVolume *mask, float cx, float cy, float cz, float sampling)
 {
   int kx,ky,kz,border;
   vlPoint3ui pos;
   float aux;

   border=(int) ceil(sampling)-1;
   aux=1.0;

   if (border<=0) {
     pos.x((int)floor(cx));
     pos.y((int)floor(cy));
     pos.z((int)floor(cz));
     mask->setVoxel(pos,aux);
   }
   else
   for(kx=-border;kx<=border;kx++)
      for(ky=-border;ky<=border;ky++)
        for(kz=-border;kz<=border;kz++) {

   pos.x((int)floor(cx+kx));
   pos.y((int)floor(cy+ky));
   pos.z((int)floor(cz+kz));
   mask->setVoxel(pos,aux);

        }


 }

 /***************AUXILIAR FUNCTIONS AND TYPES***************************/
  typedef struct listentry *LLISTPTR; /* Auxiliary pointer to list for non-recursive flood fill */
  /**
   * Auxilary structure to explore voxels in a region
   */
  typedef struct listentry {
    int      x;
    int      y;
    int      z;
    LLISTPTR next;
  } LLIST;
  /**
   * Auxiliary function to introduce a voxel in a collection for creating a region
   *
   *@param vol: Initial map
   *@param mask: result map.
   *@param tail: pointer to the end of the collection
   *@param point: voxel to include in the collection
   *@param mark: value to identify the voxels in the region in the result map.
   *@param cutoff: connected voxels with values over this cutoff will be considered in the region
   *@param cont: keep the number of voxels in the region
   *
   */
  void introduce_voxel(vlVolume *vol,vlVolume *mask,LLIST **tail, vlPoint3ui point, float cutoff,int mark,int *cont);

/********************************************************************/

 vlVolume * regions(vlVolume *vol, float cutoff, int **list, int *num_max,int min_voxels, bool erase,int connectivity )
 {
	 vlVolume *mask;
	 vlVolIter<float, vlLayout::Linear> iter(vol);
	 float value,value_mask;
	 int mark=1;
	 int cont=0;
	 int cont2=0;

	 //Creation of the mask
	 mask=new vlVolume(vol);
	 mask->clear(0.0);
	 vlVolIter<float, vlLayout::Linear> iter2(mask);
	 *list =(int*)malloc(0);
	 *num_max=0;

	 //For each voxel in original map
	 iter.begin();iter2.begin();
	 while(!iter.end())
	 {
		 value=iter.get();
		 value_mask=iter2.get();

		 //If voxel is valid and is not introduced in a region, a region is created
		 if(value>cutoff && value_mask==0)
		 {
			 //Region creation function
			 cont2=flood(vol, mask, iter.pos(), mark, cutoff, min_voxels, erase, connectivity);
			 cont++;
			 mark++;
			 //fprintf(stderr,"Nº of voxel in region %d: %d\n",cont,cont2);
			 //Set new number of voxels for region in array
			 *num_max=*num_max+1;
			 (*list)=(int*)realloc((*list),(*num_max)*sizeof(int));
			 (*list)[(*num_max)-1]=cont2;
		 }
		 iter.next();
		 iter2.next();
	 }

	 //
	 FOPS::thresholdDown(mask, 0);






	 //ORDENAR por volumen
	 	int pos_maximo,olpm, olpt;
	    int maximo, *olp, *olp2;
	 	int tem;
	 	int *list2;

	 	list2=(int*)malloc(sizeof(int)* *num_max);
	 	olp=(int*)malloc(sizeof(int)* *num_max);
	 	olp2=(int*)malloc(sizeof(int)* *num_max);
	 	memcpy(list2,(*list),sizeof(int)* *num_max);

	 	for (int i=0; i<*num_max; i++) olp[i]=i;

	 	for (int i=0; i<*num_max; i++)
	 	{
	 		pos_maximo=-1;
	 		maximo=0;
	 		for (int k=i; k<*num_max; k++)
	 		{
	 			if (list2[k]>maximo)
	 			{
	 				pos_maximo=k;
	 				maximo=list2[k];
	 				olpm=olp[k];
	 			}
	 				//printf("región %d  k %d número de voxel %d %d %d\n", i, k, pos_maximo , maximo, olpm);
	 		}
	 		tem=list2[i];
	 		olpt=olp[i];

	 		list2[i]=maximo;
	 		olp[i]=olpm;

	 		list2[pos_maximo]=tem;
	 		olp[pos_maximo]=olpt;
	 		//printf("max %d  número de voxel %d %d\n", i, pos_maximo , list2[i]);
	 	}


	 	for (int k=0; k<*num_max; k++) {
	 		// printf("%d --> %d\n", k,  olp[k]);
	 		olp2[olp[k]]=k;


	 	}
//	 	for (int k=0; k<*num_max; k++)
//	 		 printf("región %4d  número de voxel %4d %4d  %4d %4d\n", k , (*list)[k], list2[k], olp2[k], olp[k]);
//
//	 	for (int k=0; k<*num_max; k++)
//	 		 		printf("%d --> %d %d --> %d\n", k,  olp[k], k+1, olp2[k]+1  );


        int val;
	 	iter2.begin();
	 	while(!iter2.end())
	 		 {
	 			 value_mask=iter2.get();
	 			 if (value_mask>0) {
	 			 val = (int) value_mask-1.0;
	 			 iter2.set(olp2[val]+1);
	 			// printf("%d --> %d\n", val+1, olp2[val]+1);

	 			 }
	 			 iter2.next();
	 		 }

	 	for (int k=0; k<*num_max; k++) (*list)[k]=list2[k];

//		for (int k=0; k<*num_max; k++)
//		 		 		 			{
//		 		 		 			printf("region %d  número de voxel %d %d\n", k , (*list)[k], list2[k]);
//
//
//
//		 		 		 			}
//		getchar();
//
//	 	 //For each voxel in original map
//	 		 iter.begin();
//	 		 iter2.begin();
//	 		*num_max=0; cont=0; mark=1;
//	 		mask->clear(0.0);
//	 		 while(!iter.end())
//	 		 {
//	 			 value=iter.get();
//	 			 value_mask=iter2.get();
//
//	 			// fprintf(stderr,"%f %f\n",value,value_mask);
//	 			 //If voxel is valid and is not introduced in a region, a region is created
//	 			 if(value>cutoff && value_mask==0)
//	 			 {
//	 				 //Region creation function
//	 				 cont2=flood(vol, mask, iter.pos(), mark, cutoff, min_voxels, erase, connectivity);
//	 				 cont++;
//	 				 mark++;
//	 				 fprintf(stderr,"Nº of voxel in region %d: %d\n",cont,cont2);
//	 				 //Set new number of voxels for region in array
//	 				 *num_max=*num_max+1;
//	 				 // (*list)=(int*)realloc((*list),(*num_max)*sizeof(int));
//	 				 (*list)[(*num_max)-1]=cont2;
//	 			 }
//	 			 iter.next();
//	 			 iter2.next();
//	 		 }
//
//	 		for (int k=0; k<*num_max; k++)
//	 		 		 			{
//	 		 		 			printf("región-----> %d  número de voxel %d\n", k , (*list)[k]);
//	 		 		 			}


	 return mask;
 }


 int flood(vlVolume *vol, vlVolume *mask, vlPoint3ui initial_point, int mark, float cutoff, int min_voxels, bool erase, int connectivity )
 {

	 float value;
	 float value_mask;
	 int cont=1;
	 //Init collection of voxel in region
	 LLIST *list_init,*list_head, *list_tail, *list_current;
	 list_head       = (LLIST *) malloc(sizeof(LLIST));
	 list_head->x    = initial_point.x();
	 list_head->y    = initial_point.y();
	 list_head->z    = initial_point.z();
	 list_head->next = NULL;
	 list_tail       = list_head;
	 list_init		 = list_head;

	 mask->setVoxel(initial_point,(float)mark);
	 //For each voxel in the region, search all the neighbors in function of the connectivity
	 while(list_head != NULL)
	 {
			 if(list_head->x>0)
			 {
				 //Introduce neighbor in collection
				 introduce_voxel(vol,mask,&list_tail,
				vlPoint3ui(list_head->x-1,list_head->y,list_head->z),
				cutoff,mark,&cont);

				if(list_head->y>0 && connectivity>0)
				 {
					 introduce_voxel(vol,mask,&list_tail,
					  vlPoint3ui(list_head->x-1,list_head->y-1,list_head->z),
					  cutoff,mark,&cont);

					 if(list_head->z>0 && connectivity>1)
						 introduce_voxel(vol,mask,&list_tail,
						  vlPoint3ui(list_head->x-1,list_head->y-1,list_head->z-1),
						  cutoff,mark,&cont);

					 if(list_head->z<vol->dim().z()-1 && connectivity>0)
						introduce_voxel(vol,mask,&list_tail,
						 vlPoint3ui(list_head->x-1,list_head->y-1,list_head->z+1),
						 cutoff,mark,&cont);
				 }
				 if(list_head->y<vol->dim().y()-1 && connectivity>0)
				 {
					 introduce_voxel(vol,mask,&list_tail,
					  vlPoint3ui(list_head->x-1,list_head->y+1,list_head->z),
					  cutoff,mark,&cont);

					 if(list_head->z>0 && connectivity>1)
						introduce_voxel(vol,mask,&list_tail,
						 vlPoint3ui(list_head->x-1,list_head->y+1,list_head->z-1),
						 cutoff,mark,&cont);

					 if(list_head->z<vol->dim().z()-1 && connectivity>1)
						introduce_voxel(vol,mask,&list_tail,
						 vlPoint3ui(list_head->x-1,list_head->y+1,list_head->z+1),
						 cutoff,mark,&cont);
				}

			 }



			 if(list_head->x<vol->dim().x()-1)
			 {

				 introduce_voxel(vol,mask,&list_tail,
				vlPoint3ui(list_head->x+1,list_head->y,list_head->z),
				cutoff,mark,&cont);

				if(list_head->y>0 && connectivity>0)
				 {
					 introduce_voxel(vol,mask,&list_tail,
					  vlPoint3ui(list_head->x+1,list_head->y-1,list_head->z),
					  cutoff,mark,&cont);

					 if(list_head->z>0 && connectivity>1)
						 introduce_voxel(vol,mask,&list_tail,
						  vlPoint3ui(list_head->x+1,list_head->y-1,list_head->z-1),
						  cutoff,mark,&cont);

					 if(list_head->z<vol->dim().z()-1 && connectivity>1)
						introduce_voxel(vol,mask,&list_tail,
						 vlPoint3ui(list_head->x+1,list_head->y-1,list_head->z+1),
						 cutoff,mark,&cont);
				 }
				 if(list_head->y<vol->dim().y()-1 && connectivity>0)
				 {
					 introduce_voxel(vol,mask,&list_tail,
					  vlPoint3ui(list_head->x+1,list_head->y+1,list_head->z),
					  cutoff,mark,&cont);

					 if(list_head->z>0 && connectivity>1)
						introduce_voxel(vol,mask,&list_tail,
						 vlPoint3ui(list_head->x+1,list_head->y+1,list_head->z-1),
						 cutoff,mark,&cont);

					 if(list_head->z<vol->dim().z()-1 && connectivity>1)
						introduce_voxel(vol,mask,&list_tail,
						 vlPoint3ui(list_head->x+1,list_head->y+1,list_head->z+1),
						 cutoff,mark,&cont);
				}

			 }


			 if(list_head->y>0)
			 {
				 introduce_voxel(vol,mask,&list_tail,
				  vlPoint3ui(list_head->x,list_head->y-1,list_head->z),
				  cutoff,mark,&cont);

				 if(list_head->z>0 && connectivity>0)
					 introduce_voxel(vol,mask,&list_tail,
					  vlPoint3ui(list_head->x,list_head->y-1,list_head->z-1),
					  cutoff,mark,&cont);

				 if(list_head->z<vol->dim().z()-1 && connectivity>0)
					introduce_voxel(vol,mask,&list_tail,
					 vlPoint3ui(list_head->x,list_head->y-1,list_head->z+1),
					 cutoff,mark,&cont);
			 }



			 if(list_head->y<vol->dim().y()-1)
			 {
				 introduce_voxel(vol,mask,&list_tail,
				  vlPoint3ui(list_head->x,list_head->y+1,list_head->z),
				  cutoff,mark,&cont);

				 if(list_head->z>0 && connectivity>0)
					introduce_voxel(vol,mask,&list_tail,
					 vlPoint3ui(list_head->x,list_head->y+1,list_head->z-1),
					 cutoff,mark,&cont);

				 if(list_head->z<vol->dim().z()-1 && connectivity>0)
					introduce_voxel(vol,mask,&list_tail,
					 vlPoint3ui(list_head->x,list_head->y+1,list_head->z+1),
					 cutoff,mark,&cont);
			}

			 if(list_head->z>0)
				introduce_voxel(vol,mask,&list_tail,
				 vlPoint3ui(list_head->x,list_head->y,list_head->z-1),
				 cutoff,mark,&cont);

			 if(list_head->z<vol->dim().z()-1)
				introduce_voxel(vol,mask,&list_tail,
				 vlPoint3ui(list_head->x,list_head->y,list_head->z+1),
				 cutoff,mark,&cont);

			 //Once all the neighbors are explored, the voxel is deleted of the collection
			 list_current = list_head;
			 list_head    = list_head->next;
			 if(!erase)
				 free(list_current);
	 }

	 //If delete small regions is set, set to -1 small regions
	 if(erase)
	 {
		 value_mask=-1;
		 while(list_init!=NULL)
		 {
			 if(cont<min_voxels)
				 mask->setVoxel(vlPoint3ui(list_init->x,list_init->y,list_init->z),value_mask);

			 list_current = list_init;
			 list_init    = list_init->next;
			 free(list_current);
		 }
	 }

	 return cont;
 }


 void introduce_voxel(vlVolume *vol,vlVolume *mask,LLIST **tail, vlPoint3ui point, float cutoff,int mark,int *cont)
 {
	  // fprintf(stderr,"point %d %d %d\n",point.x(),point.y(),point.z());

	 /* adds new voxels which need to be visited to the linked list */
	   LLIST *current;
	   float value, value_mask=0;
	   vol->getVoxel(point,value);
	   mask->getVoxel(point,value_mask);

	   //fprintf(stderr,"point %d %d %d %f %f\n",point.x(),point.y(),point.z(),value,value_mask );
	   //If voxels is valid and it is not previously introduced in a region, it is introduced at the en of the collection
	   if(value >= cutoff && value_mask==0)
	   {
		  mask->setVoxel(point,(float)mark);

		   current = (LLIST *) malloc(sizeof(LLIST));
	   	   current->next = NULL;
	   	   current->x    = point.x();
	   	   current->y    = point.y();
	   	   current->z    = point.z();

	   	   (*tail)->next = current;
	   	   (*tail)       = current;
	   	   *cont=*cont+1;
	   }
	   else
	   {
		   //fprintf(stderr,"Saltado %f %f %f %d %d %d\n",value,value_mask,cutoff,point.x(),point.y(),point.z());
	   }
	   return;
 }

 }//FOPS


