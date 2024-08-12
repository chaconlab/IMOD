/*

 <GaussIO.c>
 for Input/Output of Gaussian Mixture Model


>> EXAMPLE OF THE FILE FOR GAUSSIAN MIXTURE MODEL << 
HEADER 3D Gaussian Mixture Model
REMARK COMMAND GaussMix -ip 1finA.pdb -ng 2
REMARK DATE Jun 25,2007 16:49:6
REMARK NGAUSS 2
REMARK logLike_final -26994.083495
HETATM    0  GAU GAU A   0     -19.800 223.788 124.704 0.571 0.571
REMARK GAUSS   0 W 0.5708323788
REMARK GAUSS   0 det 475387.1231151074
REMARK GAUSS   0 M -19.800051 223.787840 124.704350
REMARK GAUSS   0 iCovM 00 0.014196 01 0.001106 02 0.006553
REMARK GAUSS   0 iCovM 11 0.017889 12 0.001085 22 0.011367
HETATM    1  GAU GAU A   1     -20.756 203.517 112.367 0.429 0.429
REMARK GAUSS   1 W 0.4291676212
REMARK GAUSS   1 det 380806.7991351850
REMARK GAUSS   1 M -20.755753 203.517304 112.366505
REMARK GAUSS   1 iCovM 00 0.017786 01 0.004456 02 0.010503
REMARK GAUSS   1 iCovM 11 0.011111 12 0.000952 22 0.021257
TER

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
//#include "globalvar.h"
//#include "gauss.h"
//#include "BasicIO.h"
//#include "TabExpMinX.h"
//#include "Cramer3D.h"
#include "gausscorr.h"

/*** Functions (GLOBAL) ***/


/*****************/
/*** FUNCTIONS ***/
/*****************/

void Write_Gaussian3D_File(char *fname,int Ngauss,GAUSS3D *Garray,char  chain,char *comment)
{
 FILE *fp;
 int i,j,g; 
 fp = fopen(fname,"w");
 if (fp==NULL) {printf("#ERROR:Can't write to gaussfile \"%s\"\n",fname);} 
 printf("#Write_Gaussian3D_File()-->\"%s\"\n",fname);
 fprintf(fp,"HEADER 3D Gaussian Mixture Model\n");
 // Tiene un fallo: un solo "%s" para dos variables "PAR","COMMAND"
 // fprintf(fp,"REMARK COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK DATE %s\n",Get_Date_String());
 fprintf(fp,"REMARK FILENAME %s\n",fname);
 fprintf(fp,"REMARK NGAUSS %d\n",Ngauss);
 if (comment[0]!='\0') fprintf(fp,"REMARK COMMENT %s\n",comment);
 for (g=0;g<Ngauss;++g)
 {
  fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
       g,"GAU","GAU",chain,g, Garray[g].M[0],Garray[g].M[1],Garray[g].M[2],Garray[g].Weight,Garray[g].Weight);
  fprintf(fp,"REMARK GAUSS%4d W %.10lf\n",g,Garray[g].Weight);
  fprintf(fp,"REMARK GAUSS%4d det %.10lf\n",g,Garray[g].det);
  fprintf(fp,"REMARK GAUSS%4d M %lf %lf %lf\n",g,Garray[g].M[0], Garray[g].M[1], Garray[g].M[2]);
  fprintf(fp,"REMARK GAUSS%4d iCovM 00 %lf 01 %lf 02 %lf\n",
        g,Garray[g].iSigM[0][0],Garray[g].iSigM[0][1],Garray[g].iSigM[0][2]); 
  fprintf(fp,"REMARK GAUSS%4d iCovM 11 %lf 12 %lf 22 %lf\n",
        g,Garray[g].iSigM[1][1],Garray[g].iSigM[1][2],Garray[g].iSigM[2][2]); 
  fprintf(fp,"REMARK GAUSS%4d CovM  00 %lf 01 %lf 02 %lf\n",
        g,Garray[g].SigM[0][0],Garray[g].SigM[0][1],Garray[g].SigM[0][2]);
  fprintf(fp,"REMARK GAUSS%4d CovM  11 %lf 12 %lf 22 %lf\n",
        g,Garray[g].SigM[1][1],Garray[g].SigM[1][2],Garray[g].SigM[2][2]);

 } 
 fprintf(fp,"TER\n");
 fclose(fp);
} /* end of Write_Gaussian3D_File() */


char Read_Gaussian3D_File(char *fname,int *Ngauss,GAUSS3D *Garray)
{
	FILE *fp;
	int  g,g0,w,L;
	char line[512],buff[512],word[10][100],chain;
	int  Wsta[100],Wend[100],Nword;
	double det;

	*Ngauss = 0; g0 = -1;  chain = ' ';
	fp = fopen(fname,"r");
	if (fp==NULL)
	{
		printf("#ERROR:Can't open gaussfile \"%s\"\n",fname);
		exit(1);
	}
	printf("#Read_Gaussian3D_File(\"%s\")\n",fname);
	/*
>> EXAMPLE OF THE FILE FOR GAUSSIAN MIXTURE MODEL << 
HEADER 3D Gaussian Mixture Model
REMARK COMMAND GaussMix -ip 1finA.pdb -ng 2
REMARK DATE Jun 25,2007 16:49:6
REMARK NGAUSS 2
REMARK logLike_final -26994.083495
HETATM    0  GAU GAU A   0     -19.800 223.788 124.704 0.571 0.571
REMARK GAUSS   0 W 0.5708323788
REMARK GAUSS   0 det 475387.1231151074
REMARK GAUSS   0 M -19.800051 223.787840 124.704350
REMARK GAUSS   0 iCovM 00 0.014196 01 0.001106 02 0.006553
REMARK GAUSS   0 iCovM 11 0.017889 12 0.001085 22 0.011367
TER
	 */
	while (feof(fp)==0)
	{
		line[0] = '\0';
		fgets(line,511,fp);

		L = strlen(line);
		if ((L>0)&&(line[L-1]=='\n'))
			line[L-1] = '\0';

		if (strncmp(line,"HETATM",6) == 0)
			chain = line[21];

		if (strncmp(line,"REMARK GAUSS",12) == 0)
		{
			Get_Part_Of_Line(buff,line,12,512);
			Split_to_Words(buff,' ',&Nword,Wsta,Wend,100);
			for (w=0;w<Nword;++w)
				Get_Part_Of_Line(word[w],buff,Wsta[w],Wend[w]);
			g = atoi(word[0]);
			if (g0!=g)
				*Ngauss += 1;
			g0 = g;
			if (g!=(*Ngauss-1))
			{
				printf("#ERROR:bad gauss number in line g %d Ngauss %d '%s'\n",g,*Ngauss,line);
				exit(1);
			}

			if (strncmp(word[1],"W",1)==0)
				Garray[g].Weight = atof(word[2]);
			else if (strncmp(word[1],"det",3)==0)
				Garray[g].det    = atof(word[2]);
			else if (strncmp(word[1],"M",1)==0)
			{
				Garray[g].M[0] = atof(word[2]);
				Garray[g].M[1] = atof(word[3]);
				Garray[g].M[2] = atof(word[4]);
			}
			else if ((strcmp(word[1],"iCovM")==0)&&(strcmp(word[2],"00")==0))
			{
				Garray[g].iSigM[0][0] = atof(word[3]);
				Garray[g].iSigM[0][1] = Garray[g].iSigM[1][0] = atof(word[5]);
				Garray[g].iSigM[0][2] = Garray[g].iSigM[2][0] = atof(word[7]);
			}
			else if ((strcmp(word[1],"iCovM")==0)&&(strcmp(word[2],"11")==0))
			{
				Garray[g].iSigM[1][1] = atof(word[3]);
				Garray[g].iSigM[1][2] = Garray[g].iSigM[2][1] = atof(word[5]);
				Garray[g].iSigM[2][2] = atof(word[7]);
				Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].SigM, Garray[g].iSigM, &det);
				Garray[g].Cons = 1.0/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));
			}
		} /* REMARK GAUSS */

	} /* while */

	fclose(fp);

	if (*Ngauss<=0)
	{
		printf("#ERROR:gaussfile \"%s\" does not contain any gaussians.\n",fname);
		exit(1);
	}

	return(chain);

} /* end of Read_Gaussian3D_File() */

// Reads a Gaussian file allocating memory
char readGAUSS3D(char *fname,int *Ngauss,GAUSS3D **p_Garray)
{
	bool debug = true;
	FILE *fp;
	int  g,g0,w,L;
	char line[512],buff[512],word[10][100],chain;
	int  Wsta[100],Wend[100],Nword,gauss_num=0;
	double det;
	bool safety = false; // true if memory was allocated!
	GAUSS3D *Garray;

	*Ngauss = 0;
	g0 = -1;
	chain = ' ';

	fp = fopen(fname,"r");
	if (fp==NULL)
	{
		printf("Msg(readGAUSS3D): Can't open Gaussian file \"%s\"\n",fname);
		exit(1);
	}
	if(debug)
		printf("Msg(readGAUSS3D): Gaussian file name %s\n",fname);

	/*
>> EXAMPLE OF THE FILE FOR GAUSSIAN MIXTURE MODEL <<
HEADER 3D Gaussian Mixture Model
REMARK COMMAND GaussMix -ip 1finA.pdb -ng 2
REMARK DATE Jun 25,2007 16:49:6
REMARK NGAUSS 2
REMARK logLike_final -26994.083495
HETATM    0  GAU GAU A   0     -19.800 223.788 124.704 0.571 0.571
REMARK GAUSS   0 W 0.5708323788
REMARK GAUSS   0 det 475387.1231151074
REMARK GAUSS   0 M -19.800051 223.787840 124.704350
REMARK GAUSS   0 iCovM 00 0.014196 01 0.001106 02 0.006553
REMARK GAUSS   0 iCovM 11 0.017889 12 0.001085 22 0.011367
TER
	 */
	while (feof(fp)==0)
	{
		line[0] = '\0';
		fgets(line,511,fp);

		L = strlen(line);
		if ((L>0)&&(line[L-1]=='\n'))
			line[L-1] = '\0';

		// Memory allocation
		if (strncmp(line,"REMARK NGAUSS",13) == 0)
		{
			buff[0] = '\0';
			sscanf(line,"%*s %*s %s",buff);
			gauss_num = atoi(buff);
			if(gauss_num > 0)
			{
				Garray = (GAUSS3D *) malloc(sizeof(GAUSS3D) * gauss_num);
				if(!Garray)
				{
					printf("Msg(readGAUSS3D): Sorry, memory allocation failed! Forcing exit!\n");
					exit(1);
				}
				safety = true;
			}
			else
			{
				printf("Msg(readGAUSS3D): Sorry, invalid REMARK NGAUSS! Forcing exit!\n");
				exit(3);
			}
		}


		if (strncmp(line,"HETATM",6) == 0)
			chain = line[21];

		if (strncmp(line,"REMARK GAUSS",12) == 0)
		{
			if(!safety)
			{
				printf("Msg(readGAUSS3D): Sorry, REMARK NGAUSS was not found! Forcing exit!\n");
				exit(2);
			}

			Get_Part_Of_Line(buff,line,12,512);
			Split_to_Words(buff,' ',&Nword,Wsta,Wend,100);
			for (w=0;w<Nword;++w)
				Get_Part_Of_Line(word[w],buff,Wsta[w],Wend[w]);
			g = atoi(word[0]);
			if (g0!=g)
				*Ngauss += 1;
			g0 = g;
			if (g!=(*Ngauss-1))
			{
				printf("#ERROR:bad gauss number in line g %d Ngauss %d '%s'\n",g,*Ngauss,line);
				exit(1);
			}

			if (strncmp(word[1],"W",1)==0)
				Garray[g].Weight = atof(word[2]);
			else if (strncmp(word[1],"det",3)==0)
				Garray[g].det    = atof(word[2]);
			else if (strncmp(word[1],"M",1)==0)
			{
				Garray[g].M[0] = atof(word[2]);
				Garray[g].M[1] = atof(word[3]);
				Garray[g].M[2] = atof(word[4]);
			}
			else if ((strcmp(word[1],"iCovM")==0)&&(strcmp(word[2],"00")==0))
			{
				Garray[g].iSigM[0][0] = atof(word[3]);
				Garray[g].iSigM[0][1] = Garray[g].iSigM[1][0] = atof(word[5]);
				Garray[g].iSigM[0][2] = Garray[g].iSigM[2][0] = atof(word[7]);
			}
			else if ((strcmp(word[1],"iCovM")==0)&&(strcmp(word[2],"11")==0))
			{
				Garray[g].iSigM[1][1] = atof(word[3]);
				Garray[g].iSigM[1][2] = Garray[g].iSigM[2][1] = atof(word[5]);
				Garray[g].iSigM[2][2] = atof(word[7]);
				Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].SigM, Garray[g].iSigM, &det);
				Garray[g].Cons = 1.0/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));
			}
		} /* REMARK GAUSS */

	} /* while */

	fclose(fp);

	if (*Ngauss<=0)
	{
		printf("Msg(readGAUSS3D): gaussfile \"%s\" does not contain any gaussians.\n",fname);
		exit(1);
	}

	if (*Ngauss != gauss_num)
	{
		printf("Msg(readGAUSS3D): Number of gaussians mismatch (%d vs %d) in gaussfile \"%s\"\n",*Ngauss, gauss_num ,fname);
		exit(1);
	}

	if(debug)
		printf("Msg(readGAUSS3D): Number of gaussians readed %d\n",gauss_num);

	*p_Garray = Garray;
	return(chain);

} /* end of readGAUSS3D() */



void Get_Part_Of_Line(char *part,char *line,int s,int e)
{
 int i,E,L;
 L = strlen(line)-1;
 if (line[L] == '\n') L -= 1;
 if (e>L) E = L; else E = e;
 for (i=s;i<=E;++i) part[i-s] = line[i];
 part[E-s+1] = '\0';

} /* end of Get_Part_of_Line() */


//	char *str;          /* Input String */
//	char splsym;        /* Symbol for split */
//	int *Nword;         /* Number of words  */
//	int Wsta[];         /* Start point of str for a wowd */
//	int Wend[];         /* End point of str for a wowd */
//	int Nwordmax;       /* Maxium number of Nword  */
void Split_to_Words(char *str,char splsym,int *Nword,int Wsta[],int Wend[],int Nwordmax)
{
 /* [Example]
  str = "abc/d/ef//ghi//"
  splsym = '/'
   -->
  Nword 4
  (Wsta,Wend) = {(0,2), (4,4), (6,7), (10,12)}
 */
 int i,L;

 L = strlen(str);
 *Nword = 0; i = 0;

 while ((i<L)&&(*Nword < Nwordmax))
 {
  if (str[i]!=splsym)
   { Wsta[*Nword] = i;
     while ((str[i]!=splsym)&&(i<=(L-1))) { ++i; }
     Wend[*Nword] = i-1;
     ++(*Nword); }
  ++i;
 }

} /* end of Split_to_Words() */

char *Get_Date_String()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec"};
 static char string[64];
 now_t = time(NULL);
 loc_t = localtime(&now_t);

 sprintf(string,"%s %d,%d %d:%d:%d",
  Mon[loc_t->tm_mon],loc_t->tm_mday,loc_t->tm_year+1900,
  loc_t->tm_hour,loc_t->tm_min,loc_t->tm_sec);
 return(string);

} /* end of Get_Date_String() */


// Visualize a 3x3 matrix
void viewmatrix3D(double matrix[3][3],char *name)
{
	// Covar
	printf("%s %lf %lf %lf\n",name,
			matrix[0][0], matrix[0][1], matrix[0][2]);
	printf("%s %lf %lf %lf\n",name,
			matrix[1][0], matrix[1][1], matrix[1][2]);
	printf("%s %lf %lf %lf\n",name,
			matrix[2][0], matrix[2][1], matrix[2][2]);
}

// Visualize a Gaussian array
void viewgauss(GAUSS3D *GaussArray, int ngauss)
{
	for(int g=0;g<ngauss;g++)
	{
		printf(">g %d M %lf %lf %lf\n",g, GaussArray[g].M[0], GaussArray[g].M[1], GaussArray[g].M[2]);
		printf(" det %lf cons %lf weight %lf\n",
				GaussArray[g].det,
				GaussArray[g].Cons,
				GaussArray[g].Weight);
		// Covar
		viewmatrix3D(GaussArray[g].SigM,"SigM");
		// 1/Covar
		viewmatrix3D(GaussArray[g].iSigM,"iSigM");
	}
}

// Writes a VMD file with vectors representing the 3D-Gaussians in "gauss"
// Warning: the eigenvectors should be already computed using "gaussian_sizes()"
void write_gauss_evecs(GAUSS3D *gauss,int num_gauss,char *file)
{
	bool debug = true;
	FILE *Fout;
	char *color[3];
	char color1[]="red";
	char color2[]="green";
	char color3[]="blue";
	float point[3],point2[3];
	float ratio=0.3;
	float radius=0.15; // cylinder radius
	color[0]=color1;
	color[1]=color2;
	color[2]=color3;

	// Open output VMD file (common)
	if ((Fout = fopen(file, "w")) == NULL)
	{
		printf("Msg(write_gauss_sizes): Intput mode file %s error!! \n",file);
		exit(1);
	}
	fprintf(Fout,"molecule new\n");
	fprintf(Fout,"display resetview\n");

	for(int j=0; j<3; j++)
	{
		fprintf(Fout,"draw color %s\n", color[j]);
		for(int i=0; i<num_gauss; i++)
		{

//			point[0] = gauss[i].M[0] + gauss[i].evec[j][0] * (1-ratio);
//			point[1] = gauss[i].M[1] + gauss[i].evec[j][1] * (1-ratio);
//			point[2] = gauss[i].M[2] + gauss[i].evec[j][2] * (1-ratio);
//
//			// Drawing cylinder
//			fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
//					gauss[i].M[0], gauss[i].M[1], gauss[i].M[2],
//					point[0], point[1], point[2], radius);

			point2[0] = gauss[i].M[0] + gauss[i].evec[j][0];
			point2[1] = gauss[i].M[1] + gauss[i].evec[j][1];
			point2[2] = gauss[i].M[2] + gauss[i].evec[j][2];

			// Drawing cone (now it's an arrow!)
			fprintf(Fout,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
//					point[0], point[1], point[2],
					gauss[i].M[0], gauss[i].M[1], gauss[i].M[2],
					point2[0], point2[1], point2[2], 2*radius);
		}
	}

	if(debug)
		printf("Msg(write_gauss_evecs): Written VMD file %s\n",file);

	fclose(Fout);
}

// Writes a VMD file with vectors representing the 3D-Gaussians bounding boxes
// Warning: the eigenvectors should be already computed using "gaussian_sizes()"
void write_gauss_sizes(GAUSS3D *gauss,int num_gauss,char *file)
{
	bool debug = true;
	FILE *Fout;
	char *color[3];
	char color1[]="red";
	char color2[]="green";
	char color3[]="blue";
	float point[3],point2[3],size[3];
	float ratio=0.1;
	float radius=0.15; // cylinder radius
	color[0]=color1;
	color[1]=color2;
	color[2]=color3;

	// Open output VMD file (common)
	if ((Fout = fopen(file, "w")) == NULL)
	{
		printf("Msg(write_gauss_sizes): Intput mode file %s error!! \n",file);
		exit(1);
	}
	fprintf(Fout,"molecule new\n");
	fprintf(Fout,"display resetview\n");

	for(int n=0; n<3; n++)
	{
		fprintf(Fout,"draw color %s\n", color[n]);
		for(int i=0; i<num_gauss; i++)
		{
			for(int j=0; j<3; j++)
			{
				if(j==n)
					size[j] = gauss[i].size[j];
				else
					size[j] = 0.0;
			}

			point[0] = gauss[i].M[0] + size[0] * (1-ratio);
			point[1] = gauss[i].M[1] + size[1] * (1-ratio);
			point[2] = gauss[i].M[2] + size[2] * (1-ratio);

			// Drawing cylinder
			fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
					gauss[i].M[0], gauss[i].M[1], gauss[i].M[2],
					point[0], point[1], point[2], radius);

			point2[0] = gauss[i].M[0] + size[0];
			point2[1] = gauss[i].M[1] + size[1];
			point2[2] = gauss[i].M[2] + size[2];

			// Drawing cone (now it's an arrow!)
			fprintf(Fout,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
					point[0], point[1], point[2],
					point2[0], point2[1], point2[2], 3*radius);

		}
	}

	if(debug)
		printf("Msg(write_gauss_sizes): Written VMD file %s\n",file);

	fclose(Fout);
}
