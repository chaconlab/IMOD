#include "floatops.h"


void signDistanceGridMol(S_Grid g, float voxelSize,vlVolume *vol)
 {

         vlVolIter<float, vlLayout::Linear> iter(vol);
        do
         {

                if (iter.get()>0)
               {
               	   g->matrix[iter.pos().x()][iter.pos().y()][iter.pos().z()].phi = 0;
               }
                iter.next();
          } while(!iter.end());


 }

void signDistanceGridMol2(S_Grid g, float voxelSize,vlVolume *vol)
 {
         vlVolIter<float, vlLayout::Linear> iter(vol);
        do
         {

                if (iter.get()==0)
                  g->matrix[iter.pos().x()][iter.pos().y()][iter.pos().z()].phi = 0;
                iter.next();
          } while(!iter.end());

 }

namespace FOPS
{
  vlVolume *projectSurface( vlVolume *vol,float PR, bool inner )
  {
    vlVolume *vol2,*msk,*volAux;
    vlDim mask_pad;
    int dim_max;
    int n,m,k;
    int cont;

    volAux=new vlVolume(vol);

    if(volAux->dim().x()%2==1)
      mask_pad.x(1);
    else
      mask_pad.x(0);

    if(volAux->dim().y()%2==1)
      mask_pad.y(1);
    else
      mask_pad.y(0);

    if(volAux->dim().z()%2==1)
      mask_pad.z(1);
    else
      mask_pad.z(0);

    vol2=FOPS::padVolume(volAux, mask_pad,1);
    delete volAux;

    if(vol2->dim().x()>=vol2->dim().y() && vol2->dim().x()>=vol2->dim().z() )
      dim_max=vol2->dim().x();
    else
      if(vol2->dim().y()>=vol2->dim().x() && vol2->dim().y()>=vol2->dim().z() )
        dim_max=vol2->dim().y();
      else
        dim_max=vol2->dim().z();
    //dim_max*=2;

    mask_pad.x( (dim_max-vol2->dim().x())/2 );
    mask_pad.y( (dim_max-vol2->dim().y())/2 );
    mask_pad.z( (dim_max-vol2->dim().z())/2 );
    vol2=FOPS::padVolume(vol2, mask_pad,true);
    msk=new vlVolume(vol2);
    msk->clear();

    if(vol2->dim().x()!=vol2->dim().y() || vol2->dim().x()!=vol2->dim().z())
    {
      fprintf(stderr,"Error: Incorrect dimensions: %d %d %d - %d\n",
              vol2->dim().x(),vol2->dim().y(),vol2->dim().z(),dim_max);
      exit(1);
    }
    //fprintf(stderr,"createGrid dims: %d dim_max=%d\n",vol2->dim().x(),dim_max);
    //fprintf(stderr,"crear Grid dim_max=%d\n",dim_max);
    //fprintf(stderr,"PR=%f\n",PR);
    //getchar();
    S_Grid grid=  createGrid(dim_max);

    /*cont=0;
    for(m=0;m<dim_max;m++)
      for (n=0;n<dim_max;n++)
        for (k=0;k<dim_max;k++)
        {
          if(grid->matrix[m][n][k].phi==1)
          {
            cont++;
          }
        }
    fprintf(stderr,"creacion Grid=%d\n",cont);*/

    //fprintf(stderr,"Grid creado\n");
    //getchar();

    //fprintf(stderr,"signDistanceGridMol\n");
    signDistanceGridMol(grid, vol2->units().x(),vol2);

	/*cont=0;
    for(m=0;m<dim_max;m++)
      for (n=0;n<dim_max;n++)
        for (k=0;k<dim_max;k++)
        {
          if(grid->matrix[m][n][k].phi==0)
          {
            cont++;
          }
        }*/
   // fprintf(stderr,"signDistanceGridMol=%d\n",cont);


    delete vol2;

    //fprintf(stderr,"expand\n");
    expand(grid,PR);

	cont=0;
    for(m=0;m<dim_max;m++)
      for (n=0;n<dim_max;n++)
        for (k=0;k<dim_max;k++)
        {
          if(grid->matrix[m][n][k].phi==0)
          {
            cont++;
          }
        }

    fastMarching(grid, inner);


    for(m=0;m<dim_max;m++)
      for (n=0;n<dim_max;n++)
        for (k=0;k<dim_max;k++)
        {
          if(grid->matrix[m][n][k].phi==0)
          {
            msk->setVoxel(vlPoint3ui(m,n,k),(float)1.0);
          }
        }
    destroyGrid(grid);
    return(msk);

  }

S_Grid surfaceMask_grid( vlVolume *vol,float PR, bool inner )
{
  vlVolume *msk;
  vlDim mask_pad;
  int dim_max;

  msk=new vlVolume(vol);
  //msk= new vlVolume( vlDim( vol->dim().x()-20, vol->dim().y()-20, vol->dim().z()-20 ), Float, vlUnit( unit, unit, unit ) );
  msk->clear();

  if(vol->dim().x()>=vol->dim().y() && vol->dim().x()>=vol->dim().z() )
     dim_max=vol->dim().x();
   else
     if(vol->dim().y()>=vol->dim().x() && vol->dim().y()>=vol->dim().z() )
       dim_max=vol->dim().y();
     else
       dim_max=vol->dim().z();


  S_Grid grid= createGrid(dim_max);
  signDistanceGridMol(grid, vol->units().x(),vol);

  expand(grid,PR);
  shrink(grid,PR);
  fastMarching(grid, inner);
  delete msk;
  return(grid);

}

vlVolume *erode_grid( vlVolume *vol,S_Grid grid2,float PR,int offset)
{
  vlVolume *msk;
  int n,m,k;
  msk=new vlVolume(vol);
  msk->clear();
  S_Grid grid=copyGrid(grid2);

  shrink(grid,PR);

  for(m=offset;m<vol->dim().x()+offset;m++)
    for (n=offset;n<vol->dim().y()+offset;n++)
      for (k=offset;k<vol->dim().z()+offset;k++)
      {
        if((grid)->matrix[m][n][k].phi==0)
          msk->setVoxel(vlPoint3ui(m-offset,n-offset,k-offset),(float)1.0);
      }
  destroyGrid(grid);

  return(msk);

}

}




