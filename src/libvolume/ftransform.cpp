/***************************************************************************
                          ftransform.cpp  -  description
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



//****************FUNCION AUXILIAR***************************************************
//Calculo de la posicion real de un voxel en un volumen obtenido de una transformada inversa de fourier
//Para volumenes de dimension X par
vlPoint3ui realPosition(vlPoint3ui center, vlPoint3ui pos,vlPoint3ui odd)
{

	vlPoint3ui newPos;

	if(pos.x()<=center.x())
		newPos.x(pos.x()+center.x());
	else
		if(odd.x()==1)
			newPos.x(pos.x()-center.x()-1);
		else
			newPos.x(pos.x()-center.x());

	if(pos.y()<=center.y())
		newPos.y(pos.y()+center.y());
	else
		if(odd.y()==1)
			newPos.y(pos.y()-center.y()-1);
		else
			newPos.y(pos.y()-center.y());

	if(pos.z()<=center.z())
		newPos.z(pos.z()+center.z());
	else
		if(odd.z()==1)
			newPos.z(pos.z()-center.z()-1);
		else
			newPos.z(pos.z()-center.z());

	return newPos;

}


// this one seem to work ok....pablo

vlPoint3ui realPosition2(vlPoint3ui dim, vlPoint3ui posold)
{

	vlPoint3ui newPos;
	int pos[8];
	unsigned quad_boundary_x, quad_boundary_y, quad_boundary_z;
	unsigned extx_half,exty_half,extz_half,rix,riy,riz;

	extx_half  =(dim.x()-1) / 2;
	exty_half = (dim.y()-1) / 2;
	extz_half = (dim.z()-1) / 2;

	rix=posold.x();
	riy=posold.y();
	riz=posold.z();

	quad_boundary_x = extx_half+1;
	quad_boundary_y = exty_half+1;
	quad_boundary_z = extz_half+1;


	pos[0] = extx_half + rix ;
	pos[1] = rix - quad_boundary_x   ;
	pos[2] = exty_half + riy ;
	pos[3] = riy - quad_boundary_y  ;
	pos[4] = extz_half + riz  ;
	pos[5] = riz - quad_boundary_z  ;



	if (rix < quad_boundary_x) {
		if (riy < quad_boundary_y) {
			if (riz < quad_boundary_z) {
				newPos.x(pos[0]);
				newPos.y(pos[2]);
				newPos.z(pos[4]);


			} else {
				newPos.x(pos[0]);
				newPos.y(pos[2]);
				newPos.z(pos[5]);


			}
		} else {
			if (riz < quad_boundary_z) {
				newPos.x(pos[0]);
				newPos.y(pos[3]);
				newPos.z(pos[4]);

			} else {
				newPos.x(pos[0]);
				newPos.y(pos[3]);
				newPos.z(pos[5]);


			}
		}
	} else {
		if (riy < quad_boundary_y) {
			if (riz < quad_boundary_z) {
				newPos.x(pos[1]);
				newPos.y(pos[2]);
				newPos.z(pos[4]);


			} else {
				newPos.x(pos[1]);
				newPos.y(pos[2]);
				newPos.z(pos[5]);


			}
		} else {
			if (riz < quad_boundary_z) {
				newPos.x(pos[1]);
				newPos.y(pos[3]);
				newPos.z(pos[4]);

			} else {
				newPos.x(pos[1]);
				newPos.y(pos[3]);
				newPos.z(pos[5]);

			}
		}
	}


	return newPos;

}





namespace FOPS
{

//****************TRANSFORMACIONES***************************************************




vlVolume *Fourier(vlVolume * vol)
{
	float *pointer;
	int dimx,dimy,dimz;
	vlPoint3f aux;

	fftwf_complex *fout;
	fftwf_plan p;
	vlVolume *fourier;


	// Obtencion del puntero al array de datos del volumen
	pointer=(float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	//Dimensiones del array de datos
	dimx=vol->dim().x();
	dimy=vol->dim().y();
	dimz=vol->dim().z();

	//Creacion del volumen donde se alamacenara la transformada de Fourier
	fourier= new vlVolume(vlDim((dimx/2)+1,dimy,dimz),Complexf,vol->units());

	//Obtencion del puntero al array de datos de la transformada de Fourier
	fout=(fftwf_complex*)(fourier->getVoxelVoidPtr(vlPoint3ui(0,0,0)));

	//Creacion del plan para realizar la transformada. Se le pasa como argumentos
	// las dimensiones del array de entrada asi como las direcciones de memoria de los
	// array de entrada y salida
	p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
			pointer,fout,FFTW_ESTIMATE);

	// Calculo de la transformada de Fourier
	fftwf_execute(p);
	fftwf_destroy_plan(p);



	//En fourier esta el volumen en el espacio de fourier
	vol->getPosition(&aux);
	fourier->setPosition(aux);
	return fourier;
}


void Fourier(vlVolume * vol,vlVolume *fourier)
{
	float *pointer;
	int dimx,dimy,dimz;
	vlPoint3f aux;

	fftwf_complex *fout;
	fftwf_plan p;


	// Obtencion del puntero al array de datos del volumen
	pointer=(float*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	//Dimensiones del array de datos
	dimx=vol->dim().x();
	dimy=vol->dim().y();
	dimz=vol->dim().z();


	//Obtencion del puntero al array de datos de la transformada de Fourier
	fout=(fftwf_complex*)(fourier->getVoxelVoidPtr(vlPoint3ui(0,0,0)));

	//Creacion del plan para realizar la transformada. Se le pasa como argumentos
	// las dimensiones del array de entrada asi como las direcciones de memoria de los
	// array de entrada y salida
	p=fftwf_plan_dft_r2c_3d(dimz,dimy,dimx,
			pointer,fout,FFTW_ESTIMATE);

	// Calculo de la transformada de Fourier
	fftwf_execute(p);
	fftwf_destroy_plan(p);



	//En fourier esta el volumen en el espacio de fourier
	vol->getPosition(&aux);
	fourier->setPosition(aux);

}


vlVolume *iFourier(vlVolume * vol,int odd)
{
	float *fout;
	int dimx,dimy,dimz;
	vlPoint3f aux;

	fftwf_complex *pointer;
	fftwf_plan p;
	vlVolume *real;


	// Obtencion del puntero al array de datos del volumen
	pointer=(fftwf_complex*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	//Dimensiones del array de datos
	dimx=vol->dim().x();
	dimy=vol->dim().y();
	dimz=vol->dim().z();

	//Creacion del volumen donde se alamacenara la transformada inversade Fourier
	real= new vlVolume(vlDim((dimx-1)*2+odd,dimy,dimz),Float,vol->units());
	//Obtencion del puntero al array de datos de la transformada de Fourier
	fout=(float*)(real->getVoxelVoidPtr(vlPoint3ui(0,0,0)));


	p=fftwf_plan_dft_c2r_3d(dimz,dimy,(dimx-1)*2+odd,
			pointer,fout,FFTW_ESTIMATE);

	fftwf_execute(p);
	fftwf_destroy_plan(p);

	mul(real,1.0/( ((dimx-1)*2+odd)*dimy*dimz));

	vol->getPosition(&aux);
	real->setPosition(aux);
	return real;
}

void iFourier(vlVolume * vol,int odd,vlVolume *real)
{
	float *fout;
	int dimx,dimy,dimz;
	vlPoint3f aux;

	fftwf_complex *pointer;
	fftwf_plan p;



	// Obtencion del puntero al array de datos del volumen
	pointer=(fftwf_complex*)vol->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	//Dimensiones del array de datos
	dimx=vol->dim().x();
	dimy=vol->dim().y();
	dimz=vol->dim().z();

	//Obtencion del puntero al array de datos de la transformada de Fourier
	fout=(float*)(real->getVoxelVoidPtr(vlPoint3ui(0,0,0)));


	p=fftwf_plan_dft_c2r_3d(dimz,dimy,(dimx-1)*2+odd,
			pointer,fout,FFTW_ESTIMATE);

	fftwf_execute(p);
	fftwf_destroy_plan(p);

	mul(real,1.0/( ((dimx-1)*2+odd)*dimy*dimz));

	vol->getPosition(&aux);
	real->setPosition(aux);
}




vlVolume *crop(vlVolume *vol, float threshold, bool self)
{
	vlPoint3i maxCorner(-1000,-1000,-1000), minCorner(1000,1000,1000);
	vlVolume *vOut;

	int aux;
	vlVolIter<float, vlLayout::Linear> iter(vol);

	iter.begin();
	do{
		if(iter.get()>threshold)
		{
			//if(iter.get()<0.000000001)
			//std::cout<<std::endl<<"ATENCION: "<<iter.get()<<" "<<iter.pos()<<std::endl;
			aux=iter.pos().x();
			if(aux>maxCorner.x())
				maxCorner.x(aux);

			if(aux<minCorner.x())
				minCorner.x(aux);

			aux=iter.pos().y();
			if(aux>maxCorner.y())
				maxCorner.y(aux);

			if(aux<minCorner.y())
				minCorner.y(aux);

			aux=iter.pos().z();
			if(aux>maxCorner.z())
				maxCorner.z(aux);
			if(aux<minCorner.z())
				minCorner.z(aux);

		}
	}while(iter.next());

	float aux2;
	int limitx, limity,limitz;

	// force to odd limits
	if ( (maxCorner.x()-minCorner.x()) % 2 != 0)
		limitx=maxCorner.x()-minCorner.x();
	else limitx=maxCorner.x()-minCorner.x()+1;


	if ( (maxCorner.y()-minCorner.y()) % 2 != 0)
		limity=maxCorner.y()-minCorner.y();
	else limity=maxCorner.y()-minCorner.y()+1;

	if ( (maxCorner.z()-minCorner.z()) % 2 != 0)
		limitz=maxCorner.z()-minCorner.z();
	else limitz=maxCorner.z()-minCorner.z()+1;

	vlDim newD(limitx,limity, limitz);

	//std::cout << "Esquina max: "<<maxCorner<<" Esquina minima: "<<minCorner<<" New Dim: "<< newD<<std::endl;
	vOut=new vlVolume(newD,Float,vol->units());
	//std::cout << "Esquina max: "<<maxCorner<<" Esquina minima: "<<minCorner<<" New Dim: "<< newD<<std::endl;
	vlVolIter<float, vlLayout::Linear> iter2(vOut);
	vOut->clear();


	iter2.begin();
	do{
		vlPoint3ui olpP(iter2.pos().x()+minCorner.x(),
				iter2.pos().y()+minCorner.y(),
				iter2.pos().z()+minCorner.z());

		vol->getVoxel(olpP,aux2);
		iter2.set(aux2);
	}while(iter2.next());

	vlPoint3f oldPosition,newPosition;
	vol->getPosition(&oldPosition);
	newPosition.x(oldPosition.x()+(float)(minCorner.x()* vol->units().x() ));
	newPosition.y(oldPosition.y()+(float)(minCorner.y()* vol->units().y() ));
	newPosition.z(oldPosition.z()+(float)(minCorner.z()* vol->units().z() ));

	vOut->setPosition(newPosition);

	if(self)
	{
		delete vol;
	}
	return vOut;
}


vlVolume *cropad(vlVolume *vol, float threshold, int pad, bool self)
{
	vlPoint3i maxCorner(-1000,-1000,-1000), minCorner(1000,1000,1000);
	vlVolume *vOut;

	int aux;
	vlVolIter<float, vlLayout::Linear> iter(vol);

	iter.begin();
	do{
		if(iter.get()>threshold)
		{
			//if(iter.get()<0.000000001)
			//std::cout<<std::endl<<"ATENCION: "<<iter.get()<<" "<<iter.pos()<<std::endl;
			aux=iter.pos().x();
			if(aux>maxCorner.x())
				maxCorner.x(aux);

			if(aux<minCorner.x())
				minCorner.x(aux);

			aux=iter.pos().y();
			if(aux>maxCorner.y())
				maxCorner.y(aux);

			if(aux<minCorner.y())
				minCorner.y(aux);

			aux=iter.pos().z();
			if(aux>maxCorner.z())
				maxCorner.z(aux);
			if(aux<minCorner.z())
				minCorner.z(aux);

		}
	}while(iter.next());

	if(pad > 0) // Padding some voxel every dimension (+ and -)
	{
		maxCorner.x(maxCorner.x()+pad);
		minCorner.x(minCorner.x()-pad);
		maxCorner.y(maxCorner.y()+pad);
		minCorner.y(minCorner.y()-pad);
		maxCorner.z(maxCorner.z()+pad);
		minCorner.z(minCorner.z()-pad);

	}

	float aux2;
	int limitx, limity,limitz;

	// force to odd limits
	if ( (maxCorner.x()-minCorner.x()) % 2 != 0)
		limitx=maxCorner.x()-minCorner.x();
	else limitx=maxCorner.x()-minCorner.x()+1;


	if ( (maxCorner.y()-minCorner.y()) % 2 != 0)
		limity=maxCorner.y()-minCorner.y();
	else limity=maxCorner.y()-minCorner.y()+1;

	if ( (maxCorner.z()-minCorner.z()) % 2 != 0)
		limitz=maxCorner.z()-minCorner.z();
	else limitz=maxCorner.z()-minCorner.z()+1;

	vlDim newD(limitx,limity, limitz);

	//std::cout << "Esquina max: "<<maxCorner<<" Esquina minima: "<<minCorner<<" New Dim: "<< newD<<std::endl;
	vOut=new vlVolume(newD,Float,vol->units());
	//std::cout << "Esquina max: "<<maxCorner<<" Esquina minima: "<<minCorner<<" New Dim: "<< newD<<std::endl;
	vlVolIter<float, vlLayout::Linear> iter2(vOut);
	vOut->clear();


	iter2.begin();
	do{
		vlPoint3ui olpP(iter2.pos().x()+minCorner.x(),
				iter2.pos().y()+minCorner.y(),
				iter2.pos().z()+minCorner.z());

		vol->getVoxel(olpP,aux2);
		iter2.set(aux2);
	}while(iter2.next());

	vlPoint3f oldPosition,newPosition;
	vol->getPosition(&oldPosition);
	newPosition.x(oldPosition.x()+(float)(minCorner.x()* vol->units().x() ));
	newPosition.y(oldPosition.y()+(float)(minCorner.y()* vol->units().y() ));
	newPosition.z(oldPosition.z()+(float)(minCorner.z()* vol->units().z() ));

	vOut->setPosition(newPosition);

	if(self)
	{
		delete vol;
	}
	return vOut;
}


vlVolume *crop_square(vlVolume *vol, float threshold, bool self)
{
	vlPoint3i maxCorner(-1000,-1000,-1000), minCorner(1000,1000,1000);
	vlVolume *vOut;
	float threshold_neg=threshold*(-1.0);

	int aux;
	vlVolIter<float, vlLayout::Linear> iter(vol);

	iter.begin();
	do{
		if(iter.get()>threshold || iter.get()<threshold_neg)
		{
			//if(iter.get()<0.000000001)
			//std::cout<<std::endl<<"ATENCION: "<<iter.get()<<" "<<iter.pos()<<std::endl;
			aux=iter.pos().x();
			if(aux>maxCorner.x())
				maxCorner.x(aux);

			if(aux<minCorner.x())
				minCorner.x(aux);

			aux=iter.pos().y();
			if(aux>maxCorner.y())
				maxCorner.y(aux);

			if(aux<minCorner.y())
				minCorner.y(aux);

			aux=iter.pos().z();
			if(aux>maxCorner.z())
				maxCorner.z(aux);
			if(aux<minCorner.z())
				minCorner.z(aux);

		}
	}while(iter.next());

	if( ((vol->dim().x()-1) - maxCorner.x()) <= minCorner.x())
	{
		minCorner.x((vol->dim().x()-1) - (maxCorner.x()));
		maxCorner.x(((vol->dim().x()-1) - (maxCorner.x()))*2);
	}
	else
		maxCorner.x(minCorner.x()*2);

	if( ((vol->dim().y()-1) - maxCorner.y()) <= minCorner.y())
	{
		minCorner.y((vol->dim().y()-1) - (maxCorner.y()));
		maxCorner.y(((vol->dim().y()-1) - (maxCorner.y()))*2);
	}
	else
		maxCorner.y(minCorner.y()*2);

	if( ((vol->dim().z()-1) - maxCorner.z()) <= minCorner.z())
	{
		minCorner.z((vol->dim().z()-1) - (maxCorner.z()));
		maxCorner.z(((vol->dim().z()-1) - maxCorner.z())*2);
	}
	else
		maxCorner.z(minCorner.z()*2);

	float aux2;
	vlDim newD(vol->dim().x()-maxCorner.x(),
			vol->dim().y()-maxCorner.y(),
			vol->dim().z()-maxCorner.z());

	vOut=new vlVolume(newD,Float,vol->units());
	//std::cout << "Esquina max: "<<maxCorner<<" Esquina minima: "<<minCorner<<" New Dim: "<< newD<<std::endl;
	vlVolIter<float, vlLayout::Linear> iter2(vOut);
	vOut->clear();


	iter2.begin();
	do{
		vlPoint3ui olpP(iter2.pos().x()+minCorner.x(),
				iter2.pos().y()+minCorner.y(),
				iter2.pos().z()+minCorner.z());

		vol->getVoxel(olpP,aux2);
		iter2.set(aux2);
	}while(iter2.next());

	vlPoint3f oldPosition,newPosition;
	vol->getPosition(&oldPosition);
	newPosition.x(oldPosition.x()+(float)(minCorner.x()* vol->units().x() ));
	newPosition.y(oldPosition.y()+(float)(minCorner.y()* vol->units().y() ));
	newPosition.z(oldPosition.z()+(float)(minCorner.z()* vol->units().z() ));

	vOut->setPosition(newPosition);

	if(self)
	{
		delete vol;
	}
	return vOut;
}



void moveCuadrants(vlVolume **volume)
{
	vlVolume *vol=*volume;
	vlPoint3ui newP;
	float aux;
	vlPoint3f aux_pos;

	vlVolIter<float, vlLayout::Linear> iter(vol);

	vlVolume * volAux=new vlVolume(vol->dim(),Float,vol->units());

	vol->getPosition(&aux_pos);
	volAux->setPosition(aux_pos);
	iter.begin();
	do{

		aux=iter.get();
		newP=realPosition2(vol->dim(),iter.pos());
		//newP=realPosition(center,iter.pos(),odd);

		if(!volAux->setVoxel(newP,aux))
			fprintf(stderr,"Error setting a voxel\n");

	}while(iter.next());


	// flip xyz
	unsigned int tx,ty,tz;
	float *g_phi_rot;
	vlStep step;

	step=volAux->stepping();

	g_phi_rot = (float *) volAux->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	for ( tx = 0; tx < volAux->dim().x()/2.0; tx++)
		for ( ty = 0; ty < volAux->dim().y(); ty++ )
			for ( tz = 0; tz < volAux->dim().z(); tz++ )
			{

				aux=*(g_phi_rot+ step.x()*tx+ step.y()*ty+step.z()*tz);
				*(g_phi_rot+ step.x()*tx+ step.y()*ty+step.z()*tz)=

						*(g_phi_rot+
								step.x()*(volAux->dim().x()-tx-1)+
								step.y()*ty+
								step.z()*tz);

				*(g_phi_rot+ step.x()*(volAux->dim().x()-tx-1)+
						step.y()*ty+
						step.z()*tz) = aux;
			}

	for ( tx = 0; tx < volAux->dim().x(); tx++)
		for ( ty = 0; ty < volAux->dim().y()/2.0; ty++ )
			for ( tz = 0; tz < volAux->dim().z(); tz++ )
			{

				aux=*(g_phi_rot+ step.x()*tx+ step.y()*ty+step.z()*tz);
				*(g_phi_rot+ step.x()*tx+ step.y()*ty+step.z()*tz)=

						*(g_phi_rot+
								step.x()*tx+
								step.y()*(volAux->dim().y()-ty-1)+
								step.z()*tz);

				*(g_phi_rot+ step.x()*tx+
						step.y()*(volAux->dim().y()-ty-1)+
						step.z()*tz) = aux;
			}

	for ( tx = 0; tx < volAux->dim().x(); tx++)
		for ( ty = 0; ty < volAux->dim().y(); ty++ )
			for ( tz = 0; tz < volAux->dim().z()/2.0; tz++ )
			{

				aux=*(g_phi_rot+ step.x()*tx+ step.y()*ty+step.z()*tz);
				*(g_phi_rot+ step.x()*tx+ step.y()*ty+step.z()*tz)=

						*(g_phi_rot+
								step.x()*tx+
								step.y()*ty+
								step.z()*(volAux->dim().z()-tz-1));

				*(g_phi_rot+ step.x()*tx+
						step.y()*ty+
						step.z()*(volAux->dim().z()-tz-1)) = aux;
			}



	vol->~vlVolume();
	*volume=volAux;
}





}

///Namespace for operations over integer data volumes
namespace IOPS
{
void moveCuadrants(vlVolume **volume)
{
	vlVolume *vol=*volume;

	vlPoint3ui newP;
	vlPoint3ui odd;
	vlPoint3f aux_pos;
	int aux;


	vlPoint3ui center(vol->dim().x()/2+vol->dim().x()%2,
			vol->dim().y()/2+vol->dim().y()%2,
			vol->dim().z()/2+vol->dim().z()%2);



	vlVolIter<int, vlLayout::Linear> iter(vol);

	if(vol->dim().x()%2==0)
		odd.x(0);
	else
		odd.x(1);

	if(vol->dim().y()%2==0)
		odd.y(0);
	else
		odd.y(1);

	if(vol->dim().z()%2==0)
		odd.z(0);
	else
		odd.z(1);




	vlVolume * volAux=new vlVolume(vol->dim(),SignedInt32,vol->units());
	vol->getPosition(&aux_pos);
	volAux->setPosition(aux_pos);
	iter.begin();
	do{

		aux=iter.get();
		newP=realPosition(center,iter.pos(),odd);

		if(!volAux->setVoxel(newP,aux))
			fprintf(stderr,"ERROR. Wrong setting of a voxel\n");

	}while(iter.next());

	vol->~vlVolume();
	*volume=volAux;
}
}
