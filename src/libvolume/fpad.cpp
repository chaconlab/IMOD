/***************************************************************************
                          fpad.cpp  -  description
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
 #include "vlipolator_trilinear.h"



  //***************FUNCIONES AUXILIARES****************************************
 void translate(vlPoint3ui source,vlDim margin, vlPoint3ui *destiny)
 {
    destiny->x(source.x()+margin.x());
    destiny->y(source.y()+margin.y());
    destiny->z(source.z()+margin.z());

 };

 float conversion(int pos, float unit, float unit2)
 {
  float realPos= unit * (float)pos;

  return realPos/unit2;
 }




 namespace FOPS
 {

  //***************REDIMENSION DEL VOLUMEN****************************************
/*Introduce un relleno por los bordes a un volumen, devuelve puntero
a un nuevo volumen rellenado
  old: puntero al volumen a rellenar
  margin: triplete que indica el numero de voxels de los margenes en x,y y z. estos
  margenes son introducidos a ambos lados del volumen*/
  vlVolume * padVolume(vlVolume *old, vlDim margin, bool self)
  {
    vlPoint3ui newPos;
    vlPoint3f oldPosition, newPosition;

    //Creacion de volumen donde se almacenara el volumen padeado
    vlVolume *padded= new vlVolume(old,margin);
    //Inicialmente vacio
    padded->clear(0);

    //Creacion de un iterador sobre el volumen a padear
    vlVolIter<float, vlLayout::Linear> iterador(old);

    //Se recorre todo el volumen
    iterador.begin();
    while(!iterador.end())
    {
      //En funcion de la posicion de un voxel en el volumen original se calcula
      //la posicion en el volumen padeado

      translate(iterador.pos(),margin,&newPos);
      //Introducion del valor del voxel(en el original) en el voxel del volumen padeado
      padded->setVoxel(newPos,iterador.get());

       //std::cout << " A la posicion: "<<iterador.pos()<<"con valor:"<<iterador.get()<<" le corresponde la posicion: " << newPos<< std::endl;
      iterador.next();
    }
    translate(iterador.pos(),margin,&newPos);
    padded->setVoxel(newPos,iterador.get());

    // Recalculacion de la posicion del origen de coordenadas, ya que este se ha desplazado
    // al introducir los margenes
    old->getPosition(&oldPosition);
    newPosition.x(oldPosition.x()-(float)(margin.x()* old->units().x() ));
    newPosition.y(oldPosition.y()-(float)(margin.y()* old->units().y() ));
    newPosition.z(oldPosition.z()-(float)(margin.z()* old->units().z() ));
    padded->setPosition(newPosition);

    if(self)
    {
      delete old;
    }


    return padded;
  }


/*Introduce un relleno por los bordes finales de un volumen, devuelve puntero
a un nuevo volumen rellenado
  old: puntero al volumen a rellenar
  margin: triplete que indica el numero de voxels de los margenes en x,y y z. estos
  margenes son introducidos en un solo lado del volumen*/
  vlVolume * padVolume(vlVolume *old, vlDim margin, int dummy,bool self)
  {
    vlPoint3ui newPos;


    //Creacion de volumen donde se almacenara el volumen padeado
    vlVolume *padded= new vlVolume(old,margin,1);
    //Inicialmente vacio
    padded->clear(0);

    //Creacion de un iterador sobre el volumen a padear
    vlVolIter<float, vlLayout::Linear> iterador(old);

    //Se recorre todo el volumen
    iterador.begin();
    do
    {
      //Introducion del valor del voxel(en el original) en el voxel del volumen padeado
      padded->setVoxel(iterador.pos(),iterador.get());
//      std::cout << " A la posicion: "<<iterador.pos()<<" le corresponde la posicion: " << newPos<< std::endl;
    }while(iterador.next());

    if(self)
    {
      delete old;
    }
    return padded;

   }


  /*Calculo de la distancia maxima (en voxels)entre dos voxels con masa de un volumen
  vol: puntero al volumen*/
  float maxWidth(vlVolume *vol)
  {
    float max=0,length;
    //creacion de dos iteradores diferentes sobre el mismo volumen
    vlVolIter<float, vlLayout::Linear> iterador(vol);
    vlVolIter<float, vlLayout::Linear> iterador2(vol);

    bool ok=true,ok2=true;

    iterador.begin();
    while(ok)
    {

      if(iterador.get()>0.0)
      {

        ok2=true;
        iterador2.begin();
        while(ok2)
        {

          if(iterador2.get()>0.0)
          {
            //Si para dos voxels dados sus masas son mayores que 0 se calcula la
            //distancia en voxels entre ellos
            length= pow((float)iterador.pos().x() - iterador2.pos().x(),2);
            length+= pow((float)iterador.pos().y() - iterador2.pos().y(),2);
            length+= pow((float)iterador.pos().z() - iterador2.pos().z(),2);
            length=sqrt(length);
            // Si la distancia entre los dos puntos es mayor que la maxima encontrada
            // la distancia es la nueva maxima
            if(max<length)
              max=length;
          }

          if(iterador2.end())
            ok2=false;
          else
            iterador2.next();
        }
      }

      if(iterador.end())
        ok=false;
      else
        iterador.next();
    }


    return max;
  }

vlVolume *interpolate_map(vlVolume *old, vlUnit newUnit, bool self)
  {
    vlPoint3f dimension, pos,aux;
    vlDim newDim;

    // Mon: BUG 24/08/2018 "The "-1" seems a chapa...
//    dimension.x( (old->dim().x()-1)*old->units().x());
//    dimension.y( (old->dim().y()-1)*old->units().y());
//    dimension.z( (old->dim().z()-1)*old->units().z());
    dimension.x( old->dim().x()*old->units().x());
    dimension.y( old->dim().y()*old->units().y());
    dimension.z( old->dim().z()*old->units().z());


    // Mon: BUG 24/08/2018 "The "+1" seems a chapa...
//    newDim.x( (int)floor(dimension.x()/newUnit.x())+1);
//    newDim.y( (int)floor(dimension.y()/newUnit.y())+1);
//    newDim.z( (int)floor(dimension.z()/newUnit.z())+1);
    newDim.x( (int)ceil(dimension.x()/newUnit.x()) );
    newDim.y( (int)ceil(dimension.y()/newUnit.y()) );
    newDim.z( (int)ceil(dimension.z()/newUnit.z()) );
    //CHECK THIS!!!!!

    if(newDim.x()<2 || newDim.y()<2 || newDim.z()<2)
    {
    	fprintf(stderr,"ERROR: New units too large\n");
      return NULL;
    }

    vlVolume *nuevo= new vlVolume(newDim,Float,newUnit);

    vlVolIter<float, vlLayout::Linear> itern(nuevo);
    vlVolIter<float, vlLayout::Linear> itero(old);

    vlInterpolatorTriLinear<float,vlLayout::Linear> interpolador;
    itero.setInterpolation(&interpolador);


    itern.begin();
    do{
      pos.x(conversion(itern.pos().x(), nuevo->units().x(), old->units().x())) ;
      pos.y(conversion(itern.pos().y(), nuevo->units().y(), old->units().y())) ;
      pos.z(conversion(itern.pos().z(), nuevo->units().z(), old->units().z())) ;

 //     std::cout <<" CONVERSION DE: " << itern.pos() << " a "<< pos<< std::endl;
      itern.set(itero.getValueAt(pos));
    } while(itern.next());



    old->getPosition(&aux);
    nuevo->setPosition(aux);


    if(self)
      delete old;
    return nuevo;
  }




vlVolume * resize(vlVolume *old, vlDim size, bool self)
{
  vlPoint3i minCorner((old->dim().x()-size.x())/2 ,(old->dim().y()-size.y())/2,(old->dim().z()-size.z())/2);
  float aux;

  vlVolume *vOut=new vlVolume(size,Float,old->units());
  vlVolIter<float, vlLayout::Linear> iter2(vOut);
  vOut->clear();


  iter2.begin();
  do{
    vlPoint3ui olpP(iter2.pos().x()+minCorner.x(),
                  iter2.pos().y()+minCorner.y(),
                  iter2.pos().z()+minCorner.z());

    old->getVoxel(olpP,aux);
    iter2.set(aux);
  }while(iter2.next());

  vlPoint3f oldPosition,newPosition;
  old->getPosition(&oldPosition);
  newPosition.x(oldPosition.x()+(float)(minCorner.x()* old->units().x() ));
  newPosition.y(oldPosition.y()+(float)(minCorner.y()* old->units().y() ));
  newPosition.z(oldPosition.z()+(float)(minCorner.z()* old->units().z() ));

  vOut->setPosition(newPosition);

  if(self)
  {
      delete old;
  }

  return vOut;

}

// DO IT WITHOUT INTERPOLATION!!! (WHENEVER...)
/* Re-sizes a map defined by its origin and max corner (real space units). (Mon 3/3/2009)
 * The final map will be defined by its "pmin" and "pmax" corners,
 * and will be placed in the same framework as "old" one.
 *  */
vlVolume *resize(vlVolume *old, vlPoint3f pmin, vlPoint3f pmax, bool self)
{
	  vlUnit size = old->units();
	  vlDim dimnew;
	  dimnew.x( ceil( ( pmax.x() - pmin.x() ) / size.x() ) );
	  dimnew.y( ceil( ( pmax.y() - pmin.y() ) / size.y() ) );
	  dimnew.z( ceil( ( pmax.z() - pmin.z() ) / size.z() ) );
	  vlVolume *map = new vlVolume( dimnew, Float, size );
	  map->clear(0);
	  map->setPosition( pmin );
	  projectVolume(old, map);
	  if(self)
	  {
	      delete old;
	  }
	  return map;
}

//***************FIN REDIMENSION DEL VOLUMEN****************************************


	void projectVolume(vlVolume *orig, vlVolume *dest)
	{
		vlPoint3f pos_orig, pos_dest,real_pos,relative_pos;
		int x,y,z;
		float value;

		orig->getPosition(&pos_orig);
		dest->getPosition(&pos_dest);

		dest->clear(0.0);

		for(x=0;x<orig->dim().x();x++)
		{
			real_pos.x((float)x*orig->units().x() + pos_orig.x());
			relative_pos.x( (real_pos.x()-pos_dest.x()) / dest->units().x() );

			for(y=0;y<orig->dim().y();y++)
			{
				real_pos.y((float)y*orig->units().y() + pos_orig.y());
				relative_pos.y( (real_pos.y()-pos_dest.y()) / dest->units().y() );

				for(z=0;z<orig->dim().z();z++)
				{
					real_pos.z((float)z*orig->units().z() + pos_orig.z());
					relative_pos.z( (real_pos.z()-pos_dest.z()) / dest->units().z() );

					orig->getVoxel(vlPoint3ui(x,y,z),value);
					//fprintf(stderr,"%f %f %f\n",relative_pos.x(),relative_pos.y(),relative_pos.z());
					dest->interpolateVoxel(relative_pos, value);
				}
			}
		}


	}

 }
