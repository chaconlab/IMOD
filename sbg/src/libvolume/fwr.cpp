/***************************************************************************
                          fwr.cpp  -  description
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

#ifndef HAVE_STRCASECMP

//***********************FUNCION AUXILIAR ******************************

int fwr_strcasecmp(const char *s1, const char *s2)
{
	while ((*s1 != '\0')
			&& (tolower(*(unsigned char *) s1) ==
					tolower(*(unsigned char *) s2))) {
		s1++;
		s2++;
	}

	return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
}

#endif




static void swap8_aligned(void *v, long ndata) {
	int *data = (int *) v;
	long i;
	int *N;

	for (i=0; i<ndata; i++) {
		N = data + i;
		*N=(((*N >>24) & 0x000000ff) | ((*N >> 8) & 0x0000ff00) | ((*N << 8) & 0x00ff0000) | ((*N <<24) & 0xff000000) );

	}
}


int test_registration(float origx, float origy, float origz, float width)
{
	float xreg, xreg1, yreg, yreg1, zreg, zreg1;

	xreg = fabs(fmod(origx + 0.00001 * width, width));
	xreg1 = fabs(fmod(origx - 0.00001 * width, width));
	yreg = fabs(fmod(origy + 0.00001 * width, width));
	yreg1 = fabs(fmod(origy - 0.00001 * width, width));
	zreg = fabs(fmod(origz + 0.00001 * width, width));
	zreg1 = fabs(fmod(origz - 0.00001 * width, width));
	if (xreg1 < xreg) xreg = xreg1;
	if (yreg1 < yreg) yreg = yreg1;
	if (zreg1 < zreg) zreg = zreg1;
	if (xreg + yreg + zreg > 0.0001 * width) return 0;
	else return 1;
}


//(((type >>24) & 0x000000ff) | ((type >> 8) & 0x0000ff00) | ((type << 8) & 0x00ff0000) | ((type <<24) & 0xff000000) )


namespace FOPS
{

// **********************ESCRITURA*****************************************************
void write(vlVolume *vol,char *fileName)
{
	vlPoint3f aux;
	vlVolIter<float, vlLayout::Linear> iter(vol);
	vol->getPosition(&aux);

	// iter.write(fileName,&aux);

	FILE *fout;
	fout=fopen(fileName,"w");
	if(fout==NULL)
	{
		fprintf(stderr,"Error: Can't open file\n");
		return;
		// exit();
	}

	fprintf(fout," %f %f %f %f %d %d %d\n", vol->units().x(), aux.x(), aux.y(), aux.z(), vol->dim().x(), vol->dim().y(), vol->dim().z());
	fprintf(fout,"\n");

	int count=0;
	iter.begin();
	do{
		if((count+1)%10==0)
			fprintf(fout," %10.6f \n",iter.get());
		else
			fprintf(fout," %10.6f ",iter.get());
	} while(iter.next() );

	fclose(fout);
}


void writeBRIX(char *filename,vlVolume *vol,bool flag)
{
	FILE *f;
	//  bool headerr=false;


	//  int ibuff[3];
	//  float fbuff[3];
	//  long lbuff;
	float prod;
	float plus;

	//buffer de cabecera
	char header[512];
	//Buffer de datos
	char brick[512];
	float value;
	int xbrix,ybrix,zbrix;
	float realx,realy,realz;
	//min y max valores del volumen
	float min=FOPS::min_mass(vol);
	float max=FOPS::max_mass(vol);
	vlPoint3f aux;

	//Creacion de fichero
	f=fopen(filename,"wb");
	if(f==NULL)
	{
		fprintf(stderr,"writeBrix ERROR: %s Error opening file\n",filename);
		return;
	}

	memset(header,' ',512);
	plus=min;
	prod=255.0/(max-min);

	/*
      Header para brix con origenes enteros:
      sprintf(header, ":-) Origin%5.0f%5.0f%5.0f Extent%5ld%5ld%5ld Grid%5d%5d%5d "
            "Cell %10.3f%10.3f%10.3f%10.3f%10.3f%10.3f "
            "Prod%12.5f Plus%8d Sigma %12.5f",
	 */
	vol->getPosition(&aux);
	sprintf(header, ":-) Origin %12.5f %12.5f %12.5f Extent%5d%5d%5d Grid%5d%5d%5d Cell %10.3f%10.3f%10.3f%10.3f%10.3f%10.3f Prod %12.5f Plus %12.5f Sigma %12.5f",
			(aux.x()/(float)vol->units().x()),//oringinX
			(aux.y()/(float)vol->units().y()),//oringinY
			(aux.z()/(float)vol->units().z()),//oringinZ
			vol->dim().x(),vol->dim().y(),vol->dim().z(),//extend
			vol->dim().x(),vol->dim().y(),vol->dim().z(),//Grid
			vol->dim().x()*vol->units().x(),//CellX
			vol->dim().y()*vol->units().y(),//CellY
			vol->dim().z()*vol->units().z(),//CellZ
			90.0,90.0,90.0,//Cell angles
			prod,//prod
			plus,//plus
			FOPS::sigma(vol)//sigma
	);

	fwrite( header,512,sizeof(unsigned char), f);
	xbrix= (int) ceil ((float) vol->dim().x() /8.0);
	ybrix= (int) ceil ((float) vol->dim().y() /8.0);
	zbrix= (int) ceil ((float) vol->dim().z() /8.0);


	for(int b3=0;b3<zbrix;b3++)
		for(int b2=0;b2<ybrix;b2++)
			for(int b1=0;b1<xbrix;b1++)
			{

				for(int k=0;k<8;k++)
					for(int j=0;j<8;j++)
						for(int i=0;i<8;i++)
						{
							realz=b3*8+k;
							realy=b2*8+j;
							realx=b1*8+i;

							//Escritura de los valores en su posicion correspondiente del volumen
							if(realz<vol->dim().z() && realy<vol->dim().y() && realx<vol->dim().x())
							{
								vol->getVoxel(vlPoint3ui(realx,realy,realz),value);
								brick[i+j*8+k*64]=(unsigned char)((value-plus)*prod);
							}
							else
								brick[i+j*8+k*64]=(unsigned char)0;
						}
				//Se escriben "ladrillos de 8X8X8
				fwrite(brick,sizeof(char),512,f);


			}
	fclose(f);
}



void writeMRC2014(char *filename,vlVolume *vol,bool changendian,bool flag)
{
	FILE *f;
	bool headerr=false;



	int ibuff[3];
	float fbuff[3];
	float new_offset[3];
	long lbuff;
	char label[4];
	vlPoint3f aux;

	//Creacion de fichero
	f=fopen(filename,"wb");
	if(f==NULL)
	{
		fprintf(stderr,"readVol ERROR: %s opening file\n\n",filename);
		return;
	}

	//Escritura de los parametros de la cabecera
	if(flag)
		fprintf(stderr,"writeMRC> Volume dimensions: %d %d %d\n",vol->dim().x(),vol->dim().y(),vol->dim().z());
	ibuff[0]=vol->dim().x();
	ibuff[1]=vol->dim().y();
	ibuff[2]=vol->dim().z();
	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;

	ibuff[0]=2;
	if(changendian)
	{
		swap8_aligned(ibuff, 1);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume type: %d\n",ibuff[0]);


	vol->getPosition(&aux);
	/* if grid is not in register with origin, we set the CCP4 start indices to zero */
	if (test_registration(aux.x(), aux.y(), aux.z(), vol->units().x()) == 0) {
		ibuff[0]  = 0;
		ibuff[0]  = 0;
		ibuff[0]  = 0;
	} else {
		ibuff[0] = (int)floor((aux.x() / vol->units().x()+0.5));
		ibuff[1] = (int)floor((aux.y() / vol->units().x()+0.5));
		ibuff[2] = (int)floor((aux.z() / vol->units().x()+0.5));
	}

	new_offset[0] =aux.x() ;
	new_offset[1] =aux.y() ;
	new_offset[2] =aux.z() ;


	//	fbuff[0]=( aux.x() / vol->units().x());
	//	fbuff[1]=( aux.y() / vol->units().x());
	//	fbuff[2]=( aux.z() / vol->units().x());
	//
	//
	//	if ((ibuff[0]+0.00001>fbuff[0])||(ibuff[0]+0.00001<fbuff[0]))
	//	{
	//		new_offset[0]=(ibuff[0]-fbuff[0])*vol->units().x();
	//		if(flag) fprintf(stderr,"Warning: volume moved %f units on X-axis \n",new_offset[0]);
	//	}
	//	if ((ibuff[1]+0.00001>fbuff[1])||(ibuff[1]+0.00001<fbuff[1]))
	//	{
	//		new_offset[1]=(ibuff[1]-fbuff[1])*vol->units().y();
	//		if(flag) fprintf(stderr,"Warning: volume moved %f units on X-axis \n",new_offset[1]);
	//	}
	//	if ((ibuff[2]+0.00001>fbuff[2])||(ibuff[2]+0.00001<fbuff[2]))
	//	{
	//		new_offset[2]=(ibuff[2]-fbuff[2])*vol->units().z();
	//		if(flag) fprintf(stderr,"Warning: volume moved %f units on X-axis \n",new_offset[2]);
	//	}
	//
	//
	//	new_offset[0] =aux.x() + new_offset[1];
	//	new_offset[1] =aux.y() + new_offset[1];
	//	new_offset[2] =aux.z() + new_offset[2];


	//ibuff[0]=ibuff[1]=ibuff[2]=0;

	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}

	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;

	if(flag)
		fprintf(stderr,"writeMRC> Start indices: %d %d %d\n",ibuff[0],ibuff[1],ibuff[2]);

	ibuff[0]=vol->dim().x();
	ibuff[1]=vol->dim().y();
	ibuff[2]=vol->dim().z();
	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume gridsize: %d %d %d\n",ibuff[0],ibuff[1],ibuff[2]);

	fbuff[0]=vol->dim().x()*vol->units().x();
	fbuff[1]=vol->dim().y()*vol->units().y();
	fbuff[2]=vol->dim().z()*vol->units().z();
	if(changendian)
	{
		swap8_aligned(fbuff, 3);
	}
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume width: %f %f %f\n",fbuff[0],fbuff[1],fbuff[2]);



	fbuff[0]=90.0;
	fbuff[1]=90.0;
	fbuff[2]=90.0;
	if(changendian)
	{
		swap8_aligned(fbuff, 3);
	}
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume angles: %f %f %f\n",fbuff[0],fbuff[1],fbuff[2]);


	ibuff[0]=1;
	ibuff[1]=2;
	ibuff[2]=3;
	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume orientation: %d %d %d\n",ibuff[0],ibuff[1],ibuff[2]);

	fbuff[0]=min_mass(vol);
	fbuff[1]=max_mass(vol);
	fbuff[2]=calc_average(vol);
	if(changendian)
	{
		swap8_aligned(fbuff, 3);
	}
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume statistical data: %f %f %f\n",fbuff[0],fbuff[1],fbuff[2]);

	lbuff=0;
	if(changendian)
	{
		swap8_aligned(&lbuff, 1);
	}
	if(fwrite(&lbuff,4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume Space Group: %ld\n",lbuff);

	lbuff=0;
	if(changendian)
	{
		swap8_aligned(&lbuff, 1);
	}
	for(int i=0;i<2;i++)
		if(fwrite(&lbuff,4,1,f)!=1) headerr=true;

	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Number of symmetric volumes %ld\n",lbuff);

	fbuff[0]=0;
	if(changendian)
	{
		swap8_aligned(fbuff, 1);
	}
	for(int i=0;i<12;i++)
		if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Skew transformation: %f\n",fbuff[0]);

	lbuff=0;
	for(int i=38;i<50;i++)
		if(fwrite(&lbuff,4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Empty space\n");


	//vol->getPosition(&aux);
	//new_offset[0]=aux.x();
	//new_offset[1]=aux.y();
	//new_offset[2]=aux.z();

	if(flag) {
		printf("writeMRC> Volume ORIGIN %f %f %f\n", aux.x(), aux.y(), aux.z());
		printf("writeMRC> Volume ORIGIN %f %f %f\n", new_offset[0], new_offset[1], new_offset[2]);
	}

	if(changendian)
	{
		swap8_aligned(new_offset, 3);
	}

	for(int i=0;i<3;i++)
		if(fwrite(&new_offset[i],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume float point offset\n");


	label[0]='M';
	label[1]='A';
	label[2]='P';
	label[3]=' ';
	if(fwrite(label,sizeof(char),4,f)!=4) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> MAP mark\n");


	char endchar[4] = {1, 1, 0, 0};

	if (*((int *)endchar) > 65536) {   /* big endian */
		label[0]  = 17;
		label[1]  = 17;
	} else { /* little endian */
		label[0]  = 68;
		label[1]  = 65;
	}
	label[2]='\0';
	label[3]='\0';


	if(fwrite(label,sizeof(char),4,f)!=4) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> MACHST mark\n");


	fbuff[0]=FOPS::sigma(vol);
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Volume variance: %f\n",fbuff[0]);

	lbuff=0;
	if(fwrite(&lbuff,4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Number of labels employed: %ld\n",lbuff);

	label[0]='\0';
	for(int i=0;i<800;i++)
		if(fwrite(&label[0],sizeof(char),1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeMRC> Label space\n");


	//Escritura de los datos de densidad del volumen
	vlVolIter<float, vlLayout::Linear> iter(vol);
	iter.begin();
	do{
		fbuff[0]=iter.get();
		if(changendian)
		{
			swap8_aligned(fbuff, 1);
		}
		if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	}
	while(iter.next());
	if(flag)
		fprintf(stderr,"writeMRC> Volume data\n");

	fclose(f);
	if(headerr==true)
	{
		if(flag) fprintf(stdout,"writeMRC> Error writing MRC file\n");

	}
	else
	{
		if(flag) fprintf(stdout,"writeMRC> MRC file created Successfully \n");

	}
	return ;
}


void writeCCP4(char *filename,vlVolume *vol,bool changendian,bool flag)
{
	FILE *f;
	bool headerr=false;


	int ibuff[3];
	float fbuff[3];
	float new_offset[3];
	long lbuff;
	char label[4];
	vlPoint3f aux;

	//Creacion de fichero
	f=fopen(filename,"wb");
	if(f==NULL)
	{
		fprintf(stderr,"readVol ERROR: %s opening file\n\n",filename);
		return;
	}

	//Escritura de los parametros de la cabecera
	if(flag)
		fprintf(stderr,"writeCCP4> Volume dimensions: %d %d %d\n",vol->dim().x(),vol->dim().y(),vol->dim().z());
	ibuff[0]=vol->dim().x();
	ibuff[1]=vol->dim().y();
	ibuff[2]=vol->dim().z();
	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;

	ibuff[0]=2;
	if(changendian)
	{
		swap8_aligned(ibuff, 1);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume type: %d\n",ibuff[0]);

	vol->getPosition(&aux);
	fbuff[0]=( aux.x()   /vol->units().x());
	fbuff[1]=( aux.y()  /vol->units().x());
	fbuff[2]=( aux.z()  /vol->units().x());

	new_offset[0]=new_offset[1]=new_offset[2]=0;

	ibuff[0]=(int) fbuff[0];
	ibuff[1]=(int )fbuff[1];
	ibuff[2]=(int) fbuff[2];

	if(fabs(fbuff[0]-(float)ibuff[0])>0.5)
		if(fbuff[0]>0) ibuff[0]++;
		else ibuff[0]--;

	if(fabs(fbuff[1]-(float)ibuff[1])>0.5)
		if(fbuff[1]>0) ibuff[1]++;
		else ibuff[1]--;

	if(fabs(fbuff[2]-(float)ibuff[2])>0.5)
		if(fbuff[2]>0) ibuff[2]++;
		else ibuff[2]--;


	if(ibuff[0]!=fbuff[0])
	{
		new_offset[0]=((float)ibuff[0]-fbuff[0])*vol->units().x();
		if(flag) fprintf(stderr,"Warning: volume moved %f units on X-axis \n",new_offset[0]);
	}
	if(ibuff[1]!=fbuff[1])
	{
		new_offset[1]=((float)ibuff[1]-fbuff[1])*vol->units().y();
		if(flag) fprintf(stderr,"Warning: volume moved %f units on X-axis \n",new_offset[1]);
	}
	if(ibuff[2]!=fbuff[2])
	{
		new_offset[2]=((float)ibuff[2]-fbuff[2])*vol->units().z();
		if(flag) fprintf(stderr,"Warning: volume moved %f units on X-axis \n",new_offset[2]);
	}
	new_offset[0]*=-1.0;
	new_offset[1]*=-1.0;
	new_offset[2]*=-1.0;

	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume origin: %d %d %d\n",ibuff[0],ibuff[1],ibuff[2]);

	ibuff[0]=vol->dim().x();
	ibuff[1]=vol->dim().y();
	ibuff[2]=vol->dim().z();
	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume gridsize: %d %d %d\n",ibuff[0],ibuff[1],ibuff[2]);

	fbuff[0]=vol->dim().x()*vol->units().x();
	fbuff[1]=vol->dim().y()*vol->units().y();
	fbuff[2]=vol->dim().z()*vol->units().z();
	if(changendian)
	{
		swap8_aligned(fbuff, 3);
	}
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume width: %f %f %f\n",fbuff[0],fbuff[1],fbuff[2]);



	fbuff[0]=90.0;
	fbuff[1]=90.0;
	fbuff[2]=90.0;
	if(changendian)
	{
		swap8_aligned(fbuff, 3);
	}
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume angles: %f %f %f\n",fbuff[0],fbuff[1],fbuff[2]);


	ibuff[0]=1;
	ibuff[1]=2;
	ibuff[2]=3;
	if(changendian)
	{
		swap8_aligned(ibuff, 3);
	}
	if(fwrite(&ibuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&ibuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume orientation: %d %d %d\n",ibuff[0],ibuff[1],ibuff[2]);

	fbuff[0]=min_mass(vol);
	fbuff[1]=max_mass(vol);
	fbuff[2]=calc_average(vol);
	if(changendian)
	{
		swap8_aligned(fbuff, 3);
	}
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[1],4,1,f)!=1) headerr=true;
	if(fwrite(&fbuff[2],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume stadistical data: %f %f %f\n",fbuff[0],fbuff[1],fbuff[2]);

	lbuff=0;
	if(changendian)
	{
		swap8_aligned(&lbuff, 1);
	}
	if(fwrite(&lbuff,4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume Space Group: %ld\n",lbuff);

	lbuff=0;
	if(changendian)
	{
		swap8_aligned(&lbuff, 1);
	}
	for(int i=0;i<2;i++)
		if(fwrite(&lbuff,4,1,f)!=1) headerr=true;

	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Number of symetric volumes %ld\n",lbuff);

	fbuff[0]=0;
	if(changendian)
	{
		swap8_aligned(fbuff, 1);
	}
	for(int i=0;i<12;i++)
		if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Skew transformation: %f\n",fbuff[0]);

	lbuff=0;
	for(int i=0;i<11;i++)
		if(fwrite(&lbuff,4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Empty space\n");

	label[0]='A';
	label[1]='D';
	label[2]='P';
	label[3]='\0';
	if(fwrite(label,sizeof(char),4,f)!=4) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> ADP mark\n");

	if(changendian)
	{
		swap8_aligned(new_offset, 3);
	}

	for(int i=0;i<3;i++)
		if(fwrite(&new_offset[i],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume float point offset\n");


	label[0]='M';
	label[1]='A';
	label[2]='P';
	label[3]=' ';
	if(fwrite(label,sizeof(char),4,f)!=4) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> MAP mark\n");

	label[0]='\0';
	label[1]='\0';
	label[2]='\0';
	label[3]='\0';
	if(fwrite(label,sizeof(char),4,f)!=4) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> MACHST mark\n");


	fbuff[0]=FOPS::sigma(vol);
	if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Volume variance: %f\n",fbuff[0]);

	lbuff=0;
	if(fwrite(&lbuff,4,1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Number of labels employed: %ld\n",lbuff);

	label[0]='\0';
	for(int i=0;i<800;i++)
		if(fwrite(&label[0],sizeof(char),1,f)!=1) headerr=true;
	if(headerr)
		return;
	if(flag)
		fprintf(stderr,"writeCCP4> Label space\n");


	//Escritura de los datos de densidad del volumen
	vlVolIter<float, vlLayout::Linear> iter(vol);
	iter.begin();
	do{
		fbuff[0]=iter.get();
		if(changendian)
		{
			swap8_aligned(fbuff, 1);
		}
		if(fwrite(&fbuff[0],4,1,f)!=1) headerr=true;
	}
	while(iter.next());
	if(flag)
		fprintf(stderr,"writeCCP4> Volume data\n");

	fclose(f);
	if(headerr==true)
	{
		if(flag) fprintf(stdout,"writeCCP4> Error writing CCP4 file\n");

	}
	else
	{
		if(flag) fprintf(stdout,"writeCCP4> CCP4 file created Successfully \n");

	}
	return ;
}


void writeMAT5(char *filename,vlVolume *vol, char *name,int indian)
{
	FILE *f;
	int i,j,k;
	char cab[128];
	char cabElement[8];
	char arrayFlags[16];
	char dimensions[24];
	char arrayName[100];
	char cabData[8];
	float *buf_aux;
	vlPoint3ui pos;

	int textlenght=strlen(name);


	int dimx=vol->dim().x(),dimy=vol->dim().y(),dimz=vol->dim().z();


	f=fopen(filename,"wb");
	for(i=0;i<128;i++)
		cab[i]=' ';
	strcpy(cab,"MATLAB 5.0 MAT-file from Situs");

	//construccion cabecera
	if(indian==0)
	{
		cab[125]=0x01;
		cab[124]=0x00;
		cab[127]='M';
		cab[126]='I';
	}
	else
	{
		cab[124]=0x01;
		cab[125]=0x00;
		cab[126]='M';
		cab[127]='I';
	}


	//construccion array flags
	((int*)arrayFlags)[0]=6;
	((int*)arrayFlags)[1]=8;
	arrayFlags[1+8]=0x0;
	arrayFlags[8]=0x7;
	for(i=0;i<7;i++)
		arrayFlags[9+i]=0x0;

	//construccion dimensiones
	((int*)dimensions)[0]=5; //5
	((int*)dimensions)[1]=12; //12
	((int*)dimensions)[2]=dimy;
	((int*)dimensions)[3]=dimx;
	((int*)dimensions)[4]=dimz;
	((int*)dimensions)[5]=0;


	//Nombre de la matriz
	((int*)arrayName)[0]=1;
	((int*)arrayName)[1]=textlenght;
	strcpy(&arrayName[8],name);
	if(textlenght%8!=0)
		textlenght+=(8-textlenght%8);


	//Construccion de cabecera de datos
	((int*)cabData)[0]=7;
	((int*)cabData)[1]=dimx*dimy*dimz*sizeof(float);


	//construccion cabecera elemento
	((int*)cabElement)[0]=14;
	((int*)cabElement)[1]=dimx*dimy*dimz*sizeof(float)+56+textlenght;



	//Escritura en fichero
	fwrite(cab,sizeof(char),128,f);
	fwrite(cabElement,sizeof(char),8,f);
	fwrite(arrayFlags,1,16,f);
	fwrite(dimensions,1,24,f);
	fwrite(arrayName,1,8+textlenght,f);
	fwrite(cabData,1,8,f);
	pos.x(0); pos.y(0);pos.z(0);

	buf_aux=(float*)malloc(sizeof(float)*dimx*dimy*dimz);
	for(i=0;i<dimx;i++)
		for(j=0;j<dimy;j++)
			for(k=0;k<dimz;k++)
			{
				//pos.x(j),pos.y(i),pos.z(k);
				//vol->getVoxel(pos, *(buf_aux+(dimx-1-i)+((dimy-1-j)*dimx)+(k*dimx*dimy)) );

				pos.x(i),pos.y(j),pos.z(k);
				vol->getVoxel(pos, *(buf_aux+(j)+((i)*dimy)+((k)*dimx*dimy)) );
			}
	//fwrite( vol->getVoxelVoidPtr(pos),sizeof(float),dimx*dimy*dimz,f);
	fwrite( buf_aux,sizeof(float),dimx*dimy*dimz,f);

	free(buf_aux);
	fclose(f);

}

bool writeFile(vlVolume *vol,char *filename,bool ccp4Endian)
{
	int i=0;
	char *point,extension[10],nombre[100];

	point=strrchr(filename,'.');

	if(point==NULL)
	{
		fprintf(stderr,"writeFile> ERROR. Not recognized File extension %s\n",filename);
		return false;
	}

	while(point[i]!='\0')
	{
		extension[i]=point[i];
		i++;
	}
	extension[i]='\0';

	if (fwr_strcasecmp(extension,".mrc")==0)
	{
		writeMRC2014(filename,vol,ccp4Endian,false);
		return true;
	}

	if (fwr_strcasecmp(extension,".ccp4")==0)
	{
		writeCCP4(filename,vol,ccp4Endian,false);
		return true;
	}
	if(fwr_strcasecmp(extension,".brix")==0)
	{
		writeBRIX(filename,vol,true);
		return true;
	}
	if(fwr_strcasecmp(extension,".sit")==0)
	{
		write(vol,filename);
		return true;
	}
	if(fwr_strcasecmp(extension,".mat")==0)
	{

		/*nombre[0]='V';nombre[1]='o';nombre[2]='l';nombre[3]='_';
      i=0;
      while(filename[i]!='.')
      {
        nombre[i+4]=filename[i];
        i++;
      }
      nombre[i+4]='\0';*/
		i=0;
		while(filename[i]!='.')
		{
			nombre[i]=filename[i];
			i++;
		}
		nombre[i]='\0';
		if(nombre[0]>=48 &&nombre[0]<=57)
			nombre[0]='N';


		writeMAT5(filename,vol,nombre,0);
		return true;
	}

	return false;
}

// **********************FIN DE ESCRITURA*****************************************************





// **********************LECTURA*****************************************************
vlVolume* readVol(char *filename)
{
	FILE *f;

	vlVolume *vol;
	float units,posx,posy,posz;
	int dimx,dimy,dimz;
	float buffer;

	f=fopen(filename,"rt");
	if(f==NULL) {
		fprintf(stderr,"\n  Can't open %s for reading\n\n",filename);
		exit(1);
	}



	fscanf(f," %f %f %f %f %d %d %d\n", &units, &posx,&posy,&posz,&dimx,&dimy,&dimz);


	vol=new vlVolume(vlDim(dimx,dimy,dimz),Float,vlUnit(units,units,units));
	vlVolIter<float, vlLayout::Linear> iter(vol);


	iter.begin();
	do{
		fscanf(f,"%f",&buffer);
		iter.set(buffer);
	}while(iter.next());

	fclose(f);
	vol->setPosition(vlPoint3f(posx,posy,posz) );
	return vol;
}


vlVolume* readBrix(char *filename, bool flag)
{
	FILE *f;
	char keyWord[81];
	int  ext_x, ext_y, ext_z;
	float org_x, org_y, org_z;
	float grid_x, grid_y, grid_z, cell_x, cell_y, cell_z,
	cell_alpha, cell_beta, cell_gamma, prod, plus, sigma;
	int xbrix,ybrix,zbrix,realx,realy,realz;
	unsigned char brick[512];

	//Apertura de fichero
	f=fopen(filename,"rb");
	if(f==NULL) {
		fprintf(stderr,"\n  Can't open %s for reading: No such brix file\n\n",filename);
		exit(1);
	}


	// read in BRIX file header information -- stored as a 512-element char array
	fscanf(f, "%s", keyWord);
	if (fwr_strcasecmp(keyWord, ":-)") != 0) {
		fprintf(stderr, "Improperly formatted header.\n");
		return NULL;
	}

	// "Origin: the origin of the map in grid units, along X, Y, Z."
	fscanf(f, " %s %f %f %f", keyWord, &org_x, &org_y, &org_z);
	if (fwr_strcasecmp(keyWord, "origin") != 0) {
		fprintf(stderr, "Error reading origin.\n");
		return NULL;
	}

	// "Extent: the extent (size) of the map, along X, Y, Z, in grid units"
	fscanf(f, " %s %d %d %d", keyWord, &ext_x, &ext_y, &ext_z);
	if (fwr_strcasecmp(keyWord, "extent") != 0) {
		fprintf(stderr, "Error reading extent.\n");
		return NULL;
	}

	// "Grid: number of grid points along the whole unit cell, along X, Y, Z."
	fscanf(f, " %s %f %f %f", keyWord, &grid_x, &grid_y, &grid_z);
	if (fwr_strcasecmp(keyWord, "grid") != 0) {
		fprintf(stderr, "Error reading grid.\n");
		return NULL;
	}

	// "Cell: Unit cell parameters"
	// cell x, y, and z are the length of the unit cell along each axis,
	// cell alpha, beta, and cell_gamma are the angles between axes.
	fscanf(f, " %s %f %f %f %f %f %f", keyWord,
			&cell_x, &cell_y, &cell_z, &cell_alpha, &cell_beta, &cell_gamma);
	if (fwr_strcasecmp(keyWord, "cell") != 0) {
		fprintf(stderr, "Error reading cell.\n");
		return NULL;
	}

	// Convert angles to radians
	cell_alpha *= M_PI / 180.0;
	cell_beta *= M_PI / 180.0;
	cell_gamma *= M_PI / 180.0;

	// "Prod, plus: Constants that bring the electron density from byte to
	// normal scale."
	fscanf(f, " %s %f", keyWord, &prod);
	if (fwr_strcasecmp(keyWord, "prod") != 0) {
		fprintf(stderr, "Error reading prod.\n");
		return NULL;
	}
	fscanf(f, " %s %f", keyWord, &plus);
	if (fwr_strcasecmp(keyWord, "plus") != 0) {
		fprintf(stderr, "Error reading plus.\n");
		return NULL;
	}

	// "Sigma: Rms deviation of electron density map."
	fscanf(f, " %s %f", keyWord, &sigma);
	if (fwr_strcasecmp(keyWord, "sigma") != 0) {
		fprintf(stderr, "Error reading sigma.\n");
		return NULL;
	}

	if(flag)
	{
		fprintf(stderr,"readBrix> Origin: %f %f %f\n",org_x,org_y,org_z);
		fprintf(stderr,"readBrix> Extent: %d %d %d\n",ext_x,ext_y,ext_z);
		fprintf(stderr,"readBrix> Grid: %f %f %f\n",grid_x,grid_y,grid_z);
		fprintf(stderr,"readBrix> Cell: %f %f %f\n",cell_x,cell_y,cell_z);
		fprintf(stderr,"readBrix> Pros: %f plus: %f sigma: %f\n",prod,plus,sigma);
		fprintf(stderr,"readBrix> cell_alpha: %f cell_beta: %f cell_gamma: %f\n",cell_alpha,cell_beta,cell_gamma);
	}

	//Creacion del volumen de salida donde se almacenara los datos leidos
	vlVolume *vol= new vlVolume(vlDim(ext_x,ext_y,ext_z),Float,vlUnit(cell_z/ext_z,cell_y/ext_y,cell_x/ext_x));


	xbrix= (int) ceil ((float) ext_x /8.0);
	ybrix= (int) ceil ((float) ext_y /8.0);
	zbrix= (int) ceil ((float) ext_z /8.0);

	//Lectura de los datos del volumen
	fseek(f,512,SEEK_SET);
	for(int b3=0;b3<zbrix;b3++)
		for(int b2=0;b2<ybrix;b2++)
			for(int b1=0;b1<xbrix;b1++)
			{
				//Se leen "ladrillos de 8X8X8
				fread(brick,sizeof(char),512,f);

				for(int k=0;k<8;k++)
					for(int j=0;j<8;j++)
						for(int i=0;i<8;i++)
						{
							realz=b3*8+k;
							realy=b2*8+j;
							realx=b1*8+i;

							//Escritura de los valores en su posicion correspondiente del volumen
							if(realz<vol->dim().z() && realy<vol->dim().y() && realx<vol->dim().x())
								vol->setVoxel(vlPoint3ui(realx,realy,realz),(((float)brick[i+j*8+k*64])/prod)+plus);
						}
			}

	vol->setPosition(vlPoint3f(org_x*vol->units().x(),org_y*vol->units().y(),org_z*vol->units().z()) );
	fclose(f);
	return vol;
}

vlVolume* readCCP4(char *filename,bool flag)
{
	FILE *f=NULL;

	//Apertura del fichero
	f=fopen(filename,"rb");
	if(f==NULL) {
		fprintf(stderr,"\n  Can't open %s for reading: No such ccp4 file\n\n",filename);
		exit(1);
	}


	char mapString[4]="\0", symData[81]="\0";
	int origin[3], extent[3], grid[3], crs2xyz[3], mode, symBytes;
	int swap, i;
	long dataOffset;
	float cellDimensions[3], cellAngles[3], new_origin[3];
	float  xScale, yScale, zScale;
	int coord[3],x,y,z;
	float *rowdata;
	vlPoint3f position;
	bool marca_adp=false;

	// Lectura de la informacion mas importante de la cabecera
	if ( (fread(extent, sizeof(int), 3, f) != 3) ||
			(fread(&mode, sizeof(int), 1, f) != 1) ||
			(fread(origin, sizeof(int), 3, f) != 3) ||
			(fread(grid, sizeof(int), 3, f) != 3) ||
			(fread(cellDimensions, sizeof(float), 3, f) != 3) ||
			(fread(cellAngles, sizeof(float), 3, f) != 3) ||
			(fread(crs2xyz, sizeof(int), 3, f) != 3) ) {
		fprintf(stderr, "Improperly formatted line.\n");
		return NULL;
	}

	// Check the number of bytes used for storing symmetry operators
	//en principio deberia de ser 0
	fseek(f, 96, SEEK_SET);
	if ( (fread(&symBytes, sizeof(int), 1, f) != 1) ) {
		fprintf(stderr, "Problem reading the file.\n");
		return NULL;
	}

	fseek(f,192, SEEK_SET);
	if ( (fgets(mapString, 4, f) == NULL) ||
			(fwr_strcasecmp(mapString, "ADP") == 0) ) {
		marca_adp=true;
	}

	fseek(f,196, SEEK_SET);
	if((fread(new_origin,sizeof(float),3,f)!=3) ){
		fprintf(stderr, "Problem reading offset.\n");
		return NULL;
	}



	// Check for the string "MAP" at byte 208, indicating a CCP4 file.
	fseek(f, 208, SEEK_SET);
	if ( (fgets(mapString, 4, f) == NULL) ||
			(fwr_strcasecmp(mapString, "MAP") != 0) ) {
		fprintf(stderr, "File not in CCP4 format.\n");
		return NULL;
	}

	swap = 0;
	// Check the data type of the file.
	//Comprobacion de que estructura de almacenamiento se esta utilizando y adaptacion a esta
	// Si el mode es igual a 2, el almacenamiento se realiza de la forma adecuada
	if (mode != 2) {
		// fprintf(stderr,"mode=%d\n",mode);
		//Si no se ha leido el mode correctamente se comproueba si este esta invertido
		swap8_aligned(&mode, 1);
		//fprintf(stderr,"mode=%d\n",mode);

		if (mode != 2) {
			// Si una vez invertido se sigue sin leer el tipo adecuado, el fichero
			// es de un tipo no adecuado
			fprintf(stderr, "Non-real data types are currently not supported.\n");
			//return NULL;
		}
		else
		{  swap = 1;
		// fprintf(stderr,"swap done\n");
		}
	}

	// Swap all the information obtained from the header
	if (swap == 1) {
		swap8_aligned(extent, 3);
		swap8_aligned(origin, 3);
		swap8_aligned(grid, 3);
		swap8_aligned(cellDimensions, 3);
		swap8_aligned(cellAngles, 3);
		swap8_aligned(crs2xyz, 3);
		swap8_aligned(&symBytes, 1);
		swap8_aligned(new_origin, 3);

	}

	//flag=1;
	/// Presentacion de informacion leida
	if(flag)
	{
		fprintf(stderr,"readCCP4> read extent: %d %d %d\n",extent[0],extent[1],extent[2]);
		fprintf(stderr,"readCCP4> read mode: %d\n",mode);
		fprintf(stderr,"readCCP4> read origin: %d %d %d\n",origin[0],origin[1],origin[2]);
		fprintf(stderr,"readCCP4> read grid: %d %d %d\n",grid[0],grid[1],grid[2]);
		fprintf(stderr,"readCCP4> read cellDimensions: %f %f %f\n",cellDimensions[0],cellDimensions[1],cellDimensions[2]);
		fprintf(stderr,"readCCP4> read cellAngels: %f %f %f\n",cellAngles[0],cellAngles[1],cellAngles[2]);
		fprintf(stderr,"readCCP4> read crs2xyz: %d %d %d\n",crs2xyz[0],crs2xyz[1],crs2xyz[2]);
		fprintf(stderr,"readCCP4> read symBytes: %d\n",symBytes);
		if(marca_adp) fprintf(stderr,"readCCP4> read shift origin: %f %f %f\n",new_origin[0],new_origin[1],new_origin[2]);

	}
	// Check the dataOffset: this fixes the problem caused by files claiming
	// to have symmetry records when they do not.
	fseek(f, 0, SEEK_END);
	dataOffset = ftell(f) - 4*(extent[0]*extent[1]*extent[2]);

	if (fabs(origin[0]) < 0.0001 && fabs(new_origin[1]) < 0.0001 && fabs(origin[2]) < 0.0001) { /* seem to have crystallographic origin */
		fprintf(stderr,"readCCP4>  Using crystallographic (CCP4) style origin defined by unit cell start indices.\n");
		origin[0] = extent[0] * grid[0];
		origin[1] = extent[1] * grid[1];
		origin[2] = extent[2] * grid[2];
	}


	if (dataOffset != (1024 + symBytes)) {
		if (dataOffset == 1024) {
			// Bogus symmetry record information
			fprintf(stdout, "Warning: file contains bogus symmetry record.\n");
			symBytes = 0;
		}
		else {
			fprintf(stderr, "File size does not match header.\n");
			return NULL;
		}
	}


	// Read symmetry records -- organized as 80-byte lines of text.
	if (symBytes != 0) {
		fprintf(stdout, "readCCP4> Symmetry records found:\n");
		fseek(f, 1024, SEEK_SET);
		for (i = 0; i < symBytes/80; i++) {
			fgets(symData, 81, f);
			fprintf(stdout, "readCCP4> %s\n", symData);
		}
	}

	xScale = cellDimensions[0] / grid[0];
	yScale = cellDimensions[1] / grid[1];
	zScale = cellDimensions[2] / grid[2];
	if(new_origin[0]>xScale || new_origin[1]>yScale || new_origin[2]>zScale )
	{
		fprintf(stderr,"Incorrect origin shift: %f %f %f\n",
				new_origin[0],new_origin[0],new_origin[0]);
		return NULL;
	}

	//Creacion del objeto volumen que almacenara la informacion obtenida
	vlVolume* vol=new vlVolume(vlDim(extent[crs2xyz[0]-1],extent[crs2xyz[1]-1],extent[crs2xyz[2]-1]),
			Float,vlUnit(xScale,yScale,zScale));

	//Colocacion del origen de coordenadas
	if(marca_adp)
	{
		position.x(origin[0]*xScale+new_origin[0]);
		position.y(origin[1]*yScale+new_origin[1]);
		position.z(origin[2]*zScale+new_origin[2]);
	}
	else
	{
		position.x(origin[0]*xScale);
		position.y(origin[1]*yScale);
		position.z(origin[2]*zScale);
	}

	vol->setPosition(position);

	rowdata= new float[extent[0]];
	fseek(f,dataOffset, SEEK_SET);

	for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
		for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
			// Read an entire row of data from the file, then write it into the
			// datablock with the correct slice ordering.

			if (feof(f)) {
				fprintf(stderr, "Unexpected end-of-file.\n");
				return NULL;
			}

			if (ferror(f)) {
				fprintf(stderr, "Problem reading the file.\n");
				return NULL;
			}

			//Lectura de una fila
			if ( fread(rowdata, sizeof(float), extent[0], f) != (unsigned int)extent[0] ) {
				fprintf(stderr, "Error reading data row.\n");
				return NULL;
			}

			//Modificacion de la estructura de los datos en caso necesario
			if(swap==1)
				swap8_aligned(rowdata,extent[0]);

			//Escritura de los valores del buffer en su correspondiente posicion
			//en el volumen de salida
			for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
				x = coord[crs2xyz[0]-1];
				y = coord[crs2xyz[1]-1];
				z = coord[crs2xyz[2]-1];

				vol->setVoxel(vlDim(x,y,z),rowdata[coord[0]]);
			}

		}
	}

	fclose(f);
	return vol;
}

vlVolume* readMRC2014(char *filename,bool flag)
{
	FILE *f=NULL;

	//flag = 1;

	//Apertura del fichero
	f=fopen(filename,"rb");
	if(f==NULL) {
		fprintf(stderr,"\n  Can't open %s for reading: No such EM file\n\n",filename);
		exit(1);
	}

	//flag=1;

	char mapString[4]="\0", symData[81]="\0";
	int origin[3], extent[3], grid[3], crs2xyz[3], mode, symBytes;
	int swap, i;
	long dataOffset;
	float cellDimensions[3], cellAngles[3], new_origin[3];
	float  xScale, yScale, zScale;
	int coord[3],x,y,z;
	float *rowdata;
	vlPoint3f position;
	bool marca_adp=false;

	// Lectura de la informacion mas importante de la cabecera
	if ( (fread(extent, sizeof(int), 3, f) != 3) ||
			(fread(&mode, sizeof(int), 1, f) != 1) ||
			(fread(origin, sizeof(int), 3, f) != 3) ||
			(fread(grid, sizeof(int), 3, f) != 3) ||
			(fread(cellDimensions, sizeof(float), 3, f) != 3) ||
			(fread(cellAngles, sizeof(float), 3, f) != 3) ||
			(fread(crs2xyz, sizeof(int), 3, f) != 3) ) {
		fprintf(stderr, "Improperly formatted line.\n");
		return NULL;
	}

	// Check the number of bytes used for storing symmetry operators
	//en principio deberia de ser 0
	fseek(f, 96, SEEK_SET);
	if ( (fread(&symBytes, sizeof(int), 1, f) != 1) ) {
		fprintf(stderr, "Problem reading the file.\n");
		return NULL;
	}

	//    fseek(f,192, SEEK_SET);
	//    if ( (fgets(mapString, 4, f) == NULL) ||
	//         (fwr_strcasecmp(mapString, "ADP") == 0) ) {
	//           marca_adp=true;
	//    }

	fseek(f,196, SEEK_SET);
	if((fread(new_origin,sizeof(float),3,f)!=3) ){
		fprintf(stderr, "Problem reading offset.\n");
		return NULL;
	}


	// Check for the string "MAP" at byte 208, indicating a CCP4 file.
	fseek(f, 208, SEEK_SET);
	if ( (fgets(mapString, 4, f) == NULL) ||
			(fwr_strcasecmp(mapString, "MAP") != 0) ) {
		fprintf(stderr, "File not in MRC format.\n");
		return NULL;
	}

	swap = 0;
	// Check the data type of the file.
	//Comprobacion de que estructura de almacenamiento se esta utilizando y adaptacion a esta
	// Si el mode es igual a 2, el almacenamiento se realiza de la forma adecuada
	if (mode != 2) {
		fprintf(stderr,"mode=%d\n",mode);
		//Si no se ha leido el mode correctamente se comproueba si este esta invertido
		swap8_aligned(&mode, 1);
		fprintf(stderr,"mode=%d\n",mode);

		if (mode != 2) {
			// Si una vez invertido se sigue sin leer el tipo adecuado, el fichero
			// es de un tipo no adecuado
			fprintf(stderr, "Non-real data types are currently not supported.\n");
			//return NULL;
		}
		else
		{  swap = 1;
		fprintf(stderr,"swap done\n");
		}
	}

	// Swap all the information obtained from the header
	if (swap == 1) {
		swap8_aligned(extent, 3);
		swap8_aligned(origin, 3);
		swap8_aligned(grid, 3);
		swap8_aligned(cellDimensions, 3);
		swap8_aligned(cellAngles, 3);
		swap8_aligned(crs2xyz, 3);
		swap8_aligned(&symBytes, 1);
		swap8_aligned(new_origin, 3);

	}


	int nx, ny, nz;
	if (crs2xyz[0] == 1 && crs2xyz[1] == 2 && crs2xyz[2] == 3) {
		nx = 0;
		ny = 1;
		nz = 2;
	}
	else if (crs2xyz[0] == 1 && crs2xyz[1] == 3 && crs2xyz[2] == 2) {
		nx = 0;
		ny = 2;
		nz = 1;
	}
	else if (crs2xyz[0] == 2 && crs2xyz[1] == 1 && crs2xyz[2] == 3) {
		nx = 1;
		ny = 0;
		nz = 2;
	}
	else if (crs2xyz[0] == 2 && crs2xyz[1] == 3 && crs2xyz[2] == 1) {
		nx = 2;
		ny = 0;
		nz = 1;
	}
	else if (crs2xyz[0] == 3 && crs2xyz[1] == 1 && crs2xyz[2] == 2) {
		nx = 1;
		ny = 2;
		nz = 0;

	}
	else if (crs2xyz[0] == 3 && crs2xyz[1] == 2 && crs2xyz[2] == 1) {
		nx = 2;
		ny = 1;
		nz = 0;
	} else {
		fprintf(stderr,"readMRC>   Error Axis Assignment\n");
	}

	/// Presentacion de informacion leida
	if(flag)
	{
		fprintf(stderr,"readMRC> read extent: %d %d %d\n",extent[0],extent[1],extent[2]);
		fprintf(stderr,"readMRC> read mode: %d\n",mode);
		fprintf(stderr,"readMRC> read Intervals: %d %d %d\n",origin[0],origin[1],origin[2]);
		fprintf(stderr,"readMRC> read origin: %f %f %f\n",new_origin[0],new_origin[1],new_origin[2]);
		fprintf(stderr,"readMRC> read grid: %d %d %d\n",grid[0],grid[1],grid[2]);
		fprintf(stderr,"readMRC> read cellDimensions: %f %f %f\n",cellDimensions[0],cellDimensions[1],cellDimensions[2]);
		fprintf(stderr,"readMRC> read cellAngels: %f %f %f\n",cellAngles[0],cellAngles[1],cellAngles[2]);
		fprintf(stderr,"readMRC> read crs2xyz: %d %d %d\n",crs2xyz[0],crs2xyz[1],crs2xyz[2]);
		fprintf(stderr,"readMRC> read symBytes: %d\n",symBytes);

		//    if(marca_adp) fprintf(stderr,"readCCP4> read shift origin: %f %f %f\n",new_origin[0],new_origin[1],new_origin[2]);

	}
	// Check the dataOffset: this fixes the problem caused by files claiming
	// to have symmetry records when they do not.
	fseek(f, 0, SEEK_END);
	dataOffset = ftell(f) - 4*(extent[0]*extent[1]*extent[2]);




	if (dataOffset != (1024 + symBytes)) {
		if (dataOffset == 1024) {
			// Bogus symmetry record information
			fprintf(stdout, "Warning: file contains bogus symmetry record.\n");
			symBytes = 0;
		}
		else {
			fprintf(stderr, "File size does not match header.\n");
			return NULL;
		}
	}


	// Read symmetry records -- organized as 80-byte lines of text.
	if (symBytes != 0) {
		fprintf(stdout, "readCCP4> Symmetry records found:\n");
		fseek(f, 1024, SEEK_SET);
		for (i = 0; i < symBytes/80; i++) {
			fgets(symData, 81, f);
			fprintf(stdout, "readCCP4> %s\n", symData);
		}
	}

	xScale = cellDimensions[0] / grid[0];
	yScale = cellDimensions[1] / grid[1];
	zScale = cellDimensions[2] / grid[2];

//	if(new_origin[0]>xScale || new_origin[1]>yScale || new_origin[2]>zScale )
//	{
//		fprintf(stderr,"readMRC> Warning incorrect origin shift: %f %f %f\n",
//				new_origin[0],new_origin[1],new_origin[2]);
		// return NULL;
//	}

	if (fabs(new_origin[0]) < 0.0001 && fabs(new_origin[1]) < 0.0001 && fabs(new_origin[2]) < 0.0001) { /* seem to have crystallographic origin */
		fprintf(stderr,"readMRC>\nreadMRC>    --> Warning CCP4 style origin defined by unit cell start indices\nreadMRC>\n");

		new_origin[0] = origin[nx] * xScale ;
		new_origin[1] = origin[ny] * yScale ;
		new_origin[2] = origin[nz] * zScale ;
	}


	//Creacion del objeto volumen que almacenara la informacion obtenida
	vlVolume* vol=new vlVolume(vlDim(extent[nx],extent[ny],extent[nz]),
			Float,vlUnit(xScale,yScale,zScale));

	//Colocacion del origen de coordenadas
	// if(marca_adp)
	// {
	//  position.x(origin[0]*xScale+new_origin[0]);
	// position.y(origin[1]*yScale+new_origin[1]);
	// position.z(origin[2]*zScale+new_origin[2]);
	// }
	//else
	//{
	// position.x(origin[0]*xScale);
	// position.y(origin[1]*yScale);
	// position.z(origin[2]*zScale);
	//}

	position.x(new_origin[0]);
	position.y(new_origin[1]);
	position.z(new_origin[2]);
	vol->setPosition(position);

//	fprintf( stdout,"MRC> Original dimensions %d %d %d  Grid size %f \n",vol->dim().x(), vol->dim().y(),vol->dim().z());
//	fprintf(stderr,"readMRC> read extent: %d %d %d\n",extent[nx],extent[ny],extent[nz]);
//	fprintf(stderr,"readMRC> read extent: %d %d %d\n",nx,ny,nz);

	rowdata= new float[extent[0]];
	fseek(f,dataOffset, SEEK_SET);

	for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
		for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
			// Read an entire row of data from the file, then write it into the
			// datablock with the correct slice ordering.

			if (feof(f)) {
				fprintf(stderr, "Unexpected end-of-file.\n");
				return NULL;
			}

			if (ferror(f)) {
				fprintf(stderr, "Problem reading the file.\n");
				return NULL;
			}

			//Lectura de una fila
			if ( fread(rowdata, sizeof(float), extent[0], f) != (unsigned int)extent[0] ) {
				fprintf(stderr, "Error reading data row.\n");
				return NULL;
			}

			//Modificacion de la estructura de los datos en caso necesario
			if(swap==1)
				swap8_aligned(rowdata,extent[0]);

			//Escritura de los valores del buffer en su correspondiente posicion
			//en el volumen de salida
			for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
				x = coord[nx];
				y = coord[ny];
				z = coord[nz];

				vol->setVoxel(vlDim(x,y,z),rowdata[coord[0]]);
			}

		}
	}





	fclose(f);
	return vol;
}


vlVolume *readMAT5(char *filename,float unit,int indian)
{
	FILE *f;
	f=fopen(filename,"rb");
	if(f==NULL) {
		fprintf(stderr,"\n  Can't open %s for reading: No such mat file\n\n",filename);
		exit(1);
	}

	int miMatrix, tamMatrix,miUint32,mi8,miType,
	miInt32,mi12,dimx,dimy,dimz,miInt8,textLenght,
	miData,tamData,iaux;

	void *data;
	char cab[128],text[100];
	int kk;
	int cont=0,i,j,k;
	vlPoint3ui pos;
	int size;
	float aux;

	fread(cab,1,128,f);
	if(cab[125]!=0x01 || cab[124]!=0x00 || cab[127]!='M' || cab[126]!='I' )
	{
		//std::cout<<"readMAT5 ERROR: Bad header"<<std::endl;
		return NULL;
	}

	fread(&miMatrix,4,1,f);
	if(miMatrix!=14)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad element type\n");
		return NULL;
	}
	fread(&tamMatrix,4,1,f);

	fread(&miUint32,4,1,f);
	if(miUint32!=6)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad Array Flags\n");
		return NULL;
	}

	fread(&mi8,4,1,f);
	if(mi8!=8)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad Array Flags size\n");
		return NULL;
	}

	fread(&miType,4,1,f);
	//if(miType!=7)
	//{
	//  std::cout<<"readMAT5 ERROR: Bad Array Flags type"<<std::endl;
	//  return NULL;
	//}
	fread(&kk,4,1,f);

	fread(&miInt32,4,1,f);
	if(miInt32!=5)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad Dimension Array\n");
		return NULL;
	}
	fread(&mi12,4,1,f);
	if(mi12!=12)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad Dimension Array size\n");
		return NULL;
	}

	fread(&dimy,4,1,f);
	fread(&dimx,4,1,f);
	fread(&dimz,4,1,f);
	fread(&kk,4,1,f);

	fread(&miInt8,4,1,f);
	if(miInt8!=1)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad Array name\n");
		return NULL;
	}

	fread(&textLenght,4,1,f);
	if(textLenght%8==0)
		fread(&text,1,textLenght,f);
	else
		fread(&text,1,textLenght+(8-textLenght%8),f);
	text[textLenght]='\0';

	fread(&miData,4,1,f);

	switch(miData)
	{
	case 1:     iaux=sizeof(char);
	break;
	case 3:     iaux=sizeof(short int);
	break;
	case 4:     iaux=sizeof(unsigned short int);
	break;
	case 5:      iaux=sizeof(int);
	break;
	case 6:     iaux=sizeof(unsigned int);
	break;
	case 7:     iaux=sizeof(float);
	break;
	case 9:     iaux=sizeof(double);
	break;
	default:
		fprintf(stderr,"readMAT5 ERROR: Bad Data Array\n");
		return NULL;
		break;
	}


	fread(&tamData,4,1,f);
	if(tamData/iaux!=dimx*dimy*dimz)
	{
		fprintf(stderr,"readMAT5 ERROR: Bad number of Data\n");
		return NULL;
	}

	vlVolume *vol= new vlVolume(vlDim(dimx,dimy,dimz),Float,vlUnit(unit,unit,unit));


	switch(miData)
	{
	case 1:
		data=(char*)malloc(sizeof(char)*dimx*dimy*dimz);
		size=1;
		break;
	case 3:
		data=(short int*)malloc(sizeof(short int)*dimx*dimy*dimz);
		size=2;
		break;
	case 4:
		data=(unsigned short int*)malloc(sizeof(unsigned short int)*dimx*dimy*dimz);
		size=2;
		break;
	case 5:
		data=(int*)malloc(sizeof(int)*dimx*dimy*dimz);
		size=4;
		break;
	case 6:
		data=(unsigned int*)malloc(sizeof(unsigned int)*dimx*dimy*dimz);
		size=4;
		break;
	case 7:
		data=(float*)malloc(sizeof(float)*dimx*dimy*dimz);
		size=sizeof(float);
		break;
	case 9:
		data=(double*)malloc(sizeof(double)*dimx*dimy*dimz);
		size=sizeof(double);
		break;
	default:
		fprintf(stderr,"readMAT5 ERROR: Bad Data Array\n");
		return NULL;
		break;
	}

	/*if(miData==7)
    {
      pos.x(0);
      pos.y(0);
      pos.z(0);
      if(fread(vol->getVoxelVoidPtr(pos),size,dimx*dimy*dimz,f)!=dimx*dimy*dimz)
      {
        std::cout<<"readMAT5 ERROR: not enough Data"<<std::endl;
        return NULL;
      }
    }*/

	if(fread(data,size,dimx*dimy*dimz,f)!=dimx*dimy*dimz)
	{
		fprintf(stderr,"readMAT5 ERROR: not enough Data\n");
		return NULL;
	}


	for(k=0;k<dimz;k++)
	{
		pos.z(k);
		for(i=0;i<dimx;i++)
		{
			pos.x(i);
			for(j=0;j<dimy;j++)
			{
				pos.y(j);
				switch(miData)
				{
				case 1: aux=(float) (((char*)data)[cont]); break;
				case 3: aux=(float) (((short int*)data)[cont]); break;
				case 4: aux=(float) (((unsigned short int*)data)[cont]); break;
				case 5: aux=(float) (((int*)data)[cont]) ; break;
				case 6: aux=(float) (((unsigned short int*)data)[cont]) ; break;
				case 7: aux= ((float*)data)[cont]; break;
				case 9: aux=(float) (((double*)data)[cont]) ; break;
				default: break;
				}

				vol->setVoxel(pos,aux);
				cont++;
			}
		}
	}


	return vol;

	fclose(f);
}

void changeE(void *x)
{
	int i;
	char kk;
	for(i=0;i<2;i++)
	{
		kk=((char *)x)[i];
		((char *)x)[i]=((char *)x)[3-i];
		((char *)x)[3-i]=kk;
	}
}



vlVolume *readGRD(char *filename)
{
	FILE *f;
	char title[136];
	int tam;
	int iii,iii2,nbyte,intdat,i1,i2,i3;
	float a1,a2,a3,x11,x22,y11,y22,z11,z22;
	float xext,yext,zext;
	vlPoint3ui pos;
	int i,x,y,z;
	float *buffer,aux;

	f=fopen(filename,"rb");
	if(f==NULL) {
		fprintf(stderr,"\n  Can't open %s for reading: No such grd file\n\n",filename);
		exit(1);
	}



	fseek(f,4,SEEK_SET);

	/*char kk;
		int j;
		for(i=0;i<2000;i++)
		{
   	fread(&kk,1,1,f);
   	j=i%4;
   	fprintf(stderr,"%s char=(%d) %d %d\n",filename,kk,i,j);
		}
  	return NULL;*/

	fread(title,1,131,f);
	title[132]='\0';
	fseek(f,144,SEEK_SET);


	fprintf(stderr,"title=(%s)\n",title);

	fread(&iii,sizeof(int),1,f);
	changeE(&iii);
	fread(&nbyte,sizeof(int),1,f);
	changeE(&nbyte);
	fread(&intdat,sizeof(int),1,f);
	changeE(&intdat);
	fread(&xext,sizeof(float),1,f);
	changeE(&xext);
	fread(&yext,sizeof(float),1,f);
	changeE(&yext);
	fread(&zext,sizeof(float),1,f);
	changeE(&zext);
	fread(&a1,sizeof(float),1,f);
	changeE(&a1);
	fread(&a2,sizeof(float),1,f);
	changeE(&a2);
	fread(&a3,sizeof(float),1,f);
	changeE(&a3);
	fread(&x11,sizeof(float),1,f);
	changeE(&x11);
	fread(&x22,sizeof(float),1,f);
	changeE(&x22);
	fread(&y11,sizeof(float),1,f);
	changeE(&y11);
	fread(&y22,sizeof(float),1,f);
	changeE(&y22);
	fread(&z11,sizeof(float),1,f);
	changeE(&z11);
	fread(&z22,sizeof(float),1,f);
	changeE(&z22);
	fread(&i1,sizeof(int),1,f);
	changeE(&i1);
	fread(&i2,sizeof(int),1,f);
	changeE(&i2);
	fread(&i3,sizeof(int),1,f);
	changeE(&i3);
	fread(&i,sizeof(int),1,f);
	fread(&i,sizeof(int),1,f);

	i1++;
	i2++;
	i3++;

	fprintf(stderr,"iii=%d nbyte=%d intdat=%d\n",iii,nbyte,intdat);
	fprintf(stderr,"xext=%f yext=%f zext=%f\n",xext,yext,zext);
	fprintf(stderr,"a1=%f a2=%f a3=%f\n",a1,a2,a3);
	fprintf(stderr,"x11=%f x22=%f y11=%f y22=%f z11=%f z22=%f\n",x11,x22,y11,y22,z11,z22);
	fprintf(stderr,"i1=%d i2=%d i3=%d\n",i1,i2,i3);


	vlVolume *vol= new vlVolume(vlDim(i1,i2,i3),Float,vlUnit(1.0,1.0,1.0));
	vol->clear(0.0);
	buffer=(float*)malloc(sizeof(float)*i1);


	i=0;
	for(z=0;z<i3;z++)
		for(y=0;y<i2;y++)
		{
			fread(buffer,sizeof(float),i1,f);
			i=0;
			for(x=0;x<i1;x++)
			{

				aux=buffer[i];
				changeE(&aux);
				pos.x(x);pos.y(y);pos.z(z);
				vol->setVoxel(pos,aux);
				i++;
			}
			fread(&i,sizeof(int),1,f);
			fread(&i,sizeof(int),1,f);
		}


	fclose(f);

	vol->setPosition(vlPoint3f(xext*x11,yext*y11,zext*z11) );
	return vol;
}

vlVolume* readACNT(char *filename)
{
	FILE *f;
	vlVolume *vol;
	char comment[1000];
	float unitsx,unitsy,unitsz,posx,posy,posz;
	int dimx,dimy,dimz;
	float buffer;

	f=fopen(filename,"rt");
	if(f==NULL) {
		fprintf(stderr," Can't open %s for reading: No such acnt file\n\n",filename);
		exit(1);
	}

	fgets(comment,1000,f);
	//printf("%s\n",comment);

	fscanf(f,"%f %f %d\n",&posx,&unitsx,&dimx);
	fscanf(f,"%f %f %d\n",&posy,&unitsy,&dimy);
	fscanf(f,"%f %f %d\n",&posz,&unitsz,&dimz);
	fgets(comment,1000,f);



	vol=new vlVolume(vlDim(dimx,dimy,dimz),Float,vlUnit(unitsx,unitsy,unitsz));
	vlVolIter<float, vlLayout::Linear> iter(vol);


	iter.begin();
	do{
		fscanf(f,"%f\n",&buffer);
		iter.set(buffer);
	}while(iter.next());

	fclose(f);
	vol->setPosition(vlPoint3f(posx,posy,posz) );
	return vol;
}

vlVolume *readFile(char *filename,float unit)
{
	int i=0;
	char *point,extension[10];
	point=strrchr(filename,'.');
	int type=0;

	if(point!=NULL)
	{
		while(point[i]!='\0')
		{
			extension[i]=point[i];
			i++;
		}
		extension[i]='\0';

		type=7;
		if(     (fwr_strcasecmp(extension,".mrc")==0) ||
				(fwr_strcasecmp(extension,".map")==0) )
		{
			type=6;
		}
		if(fwr_strcasecmp(extension,".ccp4")==0)
		{
			type=0;
		}
		if(fwr_strcasecmp(extension,".sit")==0)
		{
			type=1;
		}
		if(fwr_strcasecmp(extension,".situs")==0)
		{
			type=1;
		}
		if(fwr_strcasecmp(extension,".brix")==0)
		{
			type=2;
		}
		if(fwr_strcasecmp(extension,".mat")==0)
		{
			type=3;
		}
		if(fwr_strcasecmp(extension,".grd")==0)
		{
			type=4;
		}
		if(fwr_strcasecmp(extension,".acnt")==0)
		{
			type=5;
		}

	}



	switch(type)
	{
	case 0:
		return readCCP4(filename,false);
	case 1:
		return readVol(filename);
	case 2:
		return readBrix(filename,true);
	case 3:
		return readMAT5(filename,unit,0);
	case 4:
		return readGRD(filename);
	case 5:
		return readACNT(filename);
	case 6:
		return readMRC2014(filename,false);
	default:
		fprintf(stderr,"Unknown volume datatype (%s). Please use .ccp4, .sit or .acnt extension\n",extension);
		exit(1);
		break;
	}



}

// **********************FIN DE LECTURA*****************************************************


}
