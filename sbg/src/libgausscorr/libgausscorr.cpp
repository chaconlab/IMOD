/*************************************************************************
 *                 libnmafit's LIBRARY: libnmafit.cpp                    *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Fitting/Molphing-related functions.                                 *
 *   (needs "libnma")                                                    *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include "gausscorr.h"

// Its code is in "MacroInfo.cpp" (libpdb)
int diagonalize_symmetric(double matrix[3][3], double eigen_vec[3][3], double eigenval[3]);

// Makes Gaussian array with standard-deviation "sd" from a Macromolecule "mol"
// (Allocating memory)
GAUSS3D *pdb2gauss(Macromolecule *mol,double sd)
{
	GAUSS3D *molgauss,*ptr;
	Tcoor pos;
	int num_atoms = mol->get_num_atoms();
	pdbIter *iter = new pdbIter(mol);

	// Allocate Gaussian array (one gaussian per atom)
	if( !(molgauss = (GAUSS3D *)malloc( sizeof(GAUSS3D) * num_atoms ) ) )
	{
		printf("Msg(pdb2gauss): Sorry, memory allocation failed\n");
		exit(1);
	}

	for( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )
	{
		(iter->get_atom())->getPosition(pos);
		ptr = &molgauss[iter->pos_atom];
		ptr->M[0] = pos[0];
		ptr->M[1] = pos[1];
		ptr->M[2] = pos[2];
		ptr->SigM[0][0] = sd;
		ptr->SigM[1][1] = sd;
		ptr->SigM[2][2] = sd;
		ptr->SigM[0][1] = 0.0;
		ptr->SigM[0][2] = 0.0;
		ptr->SigM[1][0] = 0.0;
		ptr->SigM[1][2] = 0.0;
		ptr->SigM[2][0] = 0.0;
		ptr->SigM[2][1] = 0.0;

		ptr->Weight = 1.0/(double)num_atoms;
		Cal_Inverse_Matrix3D_by_Cramer_Rule(ptr->iSigM, ptr->SigM, &(ptr->det));
		ptr->Cons = 1/(pow(2*M_PI,1.5)*sqrt(ptr->det));
	}

	delete iter;
	return molgauss;
}

// Updates Gaussian array positions "M" field from Macromolecule "mol" positions
// (Memory should be already allocated!)
void pdb2gauss_update(Macromolecule *mol, GAUSS3D *molgauss)
{
	GAUSS3D *ptr;
	Tcoor pos;
	int num_atoms = mol->get_num_atoms();
	pdbIter *iter = new pdbIter(mol);

	for( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )
	{
		(iter->get_atom())->getPosition(pos);
		ptr = &molgauss[iter->pos_atom];
		ptr->M[0] = pos[0];
		ptr->M[1] = pos[1];
		ptr->M[2] = pos[2];

	}
	delete iter;
}


// Compute the inverse of the "A" matrix into "InvA" matrix
// (Determinant of matrix A will be calculated)
void Cal_Inverse_Matrix3D_by_Cramer_Rule(double InvA[3][3],double A[3][3],double *Det)
{
 /*
 <Cramer's Rule>
 Inv[A] = 1/|A| transpose[Delta_ij]
 Delta_ij = (-1)**(i+j) * |Aij|
 Aij = matrix removing i-th row and j-th column.
 */
 int i,j;
 double det;

 det =  A[0][0]*A[1][1]*A[2][2]
      + A[0][2]*A[1][0]*A[2][1]
      + A[0][1]*A[1][2]*A[2][0]
      - A[0][2]*A[1][1]*A[2][0]
      - A[0][0]*A[1][2]*A[2][1]
      - A[0][1]*A[1][0]*A[2][2];

 /* printf("#det = %lf\n",det); */

 InvA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]);
 InvA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]);
 InvA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]);

 InvA[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
 InvA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]);
 InvA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]);

 InvA[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
 InvA[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]);
 InvA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]);

 for (i=0;i<3;++i)
	 for (j=0;j<3;++j)
		 InvA[i][j] /= det;

 *Det = det;

} /* end of  Cal_Inverse_Matrix3D_by_Cramer_Rule() */

// Compute the inverse of the "A" matrix into "InvA" matrix (symmetric)
// (Determinant of matrix A will be calculated)
void Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(double InvA[3][3],double A[3][3],double *Det)
{
 /*
 <Cramer's Rule>
 Inv[A] = 1/|A| transpose[Delta_ij]
 Delta_ij = (-1)**(i+j) * |Aij|
 Aij = matrix removing i-th row and j-th column.
 */
 int i,j;
 double det;


 det =  A[0][0]*A[1][1]*A[2][2] + A[0][2]*A[1][0]*A[2][1]
      + A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0]
      - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2];

 InvA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1])/det;
 InvA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1])/det;
 InvA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1])/det;

 InvA[1][0] = InvA[0][1];
 InvA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0])/det;
 InvA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0])/det;

 InvA[2][0] = InvA[0][2];
 InvA[2][1] = InvA[1][2];
 InvA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0])/det;

 *Det = det;
} /* end of  Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric() */


double Overlap_Integral_Bwn_Two_Gauss3Ds(GAUSS3D *gA,GAUSS3D *gB)
{
	int i;
	double K[3][3],invK[3][3],DetK,ov;
	double nm[3];
	double qform;

	/*
   Overlap = 1/[(2pi)**3/2 * sqrt(|P+Q|)] *
             exp[-1/2 * tr_(n-m) inv(P+Q) (n-m)]
            = 1/(2pi)**3/2 /sqrt(|K|)*exp[-1/2 * (n-m) invK (n-m)]
    where P = gA->SigM, Q = gB->SigM, m = gA->M, n = gB->M.
     K   = P + Q
  invK   = inv(P + Q)
	 */

	Add_Matrix3D_Symmetric(K,gA->SigM,gB->SigM);
	Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(invK,K,&DetK);

//	viewmatrix3D(K,"K");
//	viewmatrix3D(invK,"invK");

	for (i=0;i<3;++i)
		nm[i] =  gB->M[i] - gA->M[i];
	qform =  Quadratic_Form_3D(invK,nm);

//	printf("nm= %f %f %f  DetK= %f  qform= %f  one_over_2pi_32= %f\n",nm[0],nm[1],nm[2],DetK,qform,one_over_2pi_32);

	ov = exp(-0.5*qform) * one_over_2pi_32 / sqrt(DetK);
//	ov = exp(-0.5*qform) * (one_over_2pi_32 / sqrt(DetK) );
	return(ov);
} /* end of double Overlap_Integral_Bwn_Two_Gauss3Ds() */

//// Computes the Overlap Integral (Correlation) between two gaussian sets
//double Overlap_Integral_Bwn_GAUSS_MOLECULEs(GAUSS3D *Ga,GAUSS3D *Gb,int nga,int ngb,double factor)
//{
//	bool debug = true;
//	double sfact=2.0;
//	int i,j;
//	double ov, ov_all;
//	float dist[3];
//
//	if(factor == 0.0)
//		sfact = 9999999.0; // the initial factor is "exact"
//
//	ov_all = 0.0;
//	for (i=0; i<nga; i++)
//		for (j=0; j<ngb; j++)
//		{
//			// distance check
//			dist[0] = Ga[i].M[0] - Gb[j].M[0];
//			dist[1] = Ga[i].M[1] - Gb[j].M[1];
//			dist[2] = Ga[i].M[2] - Gb[j].M[2];
//			if( sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]) < sfact*Ga[i].SigM[0][0])
//			{
//				ov = Overlap_Integral_Bwn_Two_Gauss3Ds(&(Ga[i]),&(Gb[j]));
//				ov_all +=  ov;
//			}
//
////			if(debug)
////				printf("Ga[%d].weight= %8.6f  Gb[%d].weight= %8.6f  ov= %8.6f\n",i,Ga[i].Weight,j,Gb[j].Weight,ov);
//		}
//
//	if(factor == 0.0)
//		return(ov_all);
//	else
//		return(1.0-(ov_all/factor));
//} /* end of Overlap_Integral_Bwn_GAUSS_MOLECULEs() */

//// Computes the Overlap Integral (Correlation) between two gaussian sets
//double Overlap_Integral_Bwn_GAUSS_MOLECULEs(GAUSS3D *Ga,GAUSS3D *Gb,int nga,int ngb,double factor)
//{
//	bool debug = true;
//	int i,j;
//	double ov, ov_all;
//
//	ov_all = 0.0;
//	for (i=0; i<nga; i++)
//		for (j=0; j<ngb; j++)
//		{
//			ov = Overlap_Integral_Bwn_Two_Gauss3Ds(&(Ga[i]),&(Gb[j]));
////			ov_all += Ga[i].Weight * Gb[j].Weight * ov;
//			ov_all +=  ov;
//
////			if(debug)
////				printf("Ga[%d].weight= %8.6f  Gb[%d].weight= %8.6f  ov= %8.6f\n",i,Ga[i].Weight,j,Gb[j].Weight,ov);
//		}
//
//	if(factor == 0.0)
//		return(ov_all);
//	else
//		return(1.0-(ov_all/factor));
//} /* end of Overlap_Integral_Bwn_GAUSS_MOLECULEs() */

// Computes the Overlap Integral (Correlation) between two gaussian sets
// A-gaussian array define bounding boxes
double Overlap_Integral_Bwn_GAUSS_MOLECULEs(GAUSS3D *Ga,GAUSS3D *Gb,int nga,int ngb,double factor)
{
	bool debug = true;
	int i,j,num=0;
	double ov, ov_all;
	float bound_high[3],bound_low[3];

	ov_all = 0.0;
	for (i=0; i<nga; i++)
	{
		// A-i gaussian position (define bounding boxes size)
		// Warning: ".size[]" should be positive!
		bound_high[0] = Ga[i].M[0] + Ga[i].size[0];
		bound_high[1] = Ga[i].M[1] + Ga[i].size[1];
		bound_high[2] = Ga[i].M[2] + Ga[i].size[2];
		bound_low[0] = Ga[i].M[0] - Ga[i].size[0];
		bound_low[1] = Ga[i].M[1] - Ga[i].size[1];
		bound_low[2] = Ga[i].M[2] - Ga[i].size[2];

		for (j=0; j<ngb; j++)
		{
			if( Gb[j].M[0] < bound_high[0] && Gb[j].M[0] > bound_low[0] &&
				 Gb[j].M[1] < bound_high[1] && Gb[j].M[1] > bound_low[1] &&
					Gb[j].M[2] < bound_high[2] && Gb[j].M[2] > bound_low[2] )
					{

						ov = Overlap_Integral_Bwn_Two_Gauss3Ds(&(Ga[i]),&(Gb[j]));
						//			ov_all += Ga[i].Weight * Gb[j].Weight * ov;
						ov_all +=  ov;

						//			if(debug)
						//				printf("Ga[%d].weight= %8.6f  Gb[%d].weight= %8.6f  ov= %8.6f\n",i,Ga[i].Weight,j,Gb[j].Weight,ov);
						num++; // counts evaluated pairs
					}
		}
	}

	if(debug)
		printf("max pairs: %d  current pairs: %d  reduced_to: %f%%\n",nga*ngb,num,100.0*(float)num/(float)(nga*ngb));

	if(factor == 0.0)
		return(ov_all);
	else
		return(1.0-(ov_all/factor));
} /* end of Overlap_Integral_Bwn_GAUSS_MOLECULEs() */

// Adds two symmetric matrices  C = A + B
void Add_Matrix3D_Symmetric(double C[3][3],double A[3][3],double B[3][3])
{
 C[0][0] = A[0][0]+B[0][0];
 C[0][1] = A[0][1]+B[0][1];
 C[0][2] = A[0][2]+B[0][2];

 C[1][0] = C[0][1];
 C[1][1] = A[1][1]+B[1][1];
 C[1][2] = A[1][2]+B[1][2];

 C[2][0] = C[0][2];
 C[2][1] = C[1][2];
 C[2][2] = A[2][2]+B[2][2];
} /* end of Add_Matrix3D_Symmetric() */

// C = A + B
void Add_Matrix3D(double C[3][3],double A[3][3],double B[3][3])
{
 int i,j;
 for (i=0;i<3;++i)
  for (j=0;j<3;++j) C[i][j] = A[i][j] + B[i][j];
} /* end of Add_Matrix3D() */

// Compute: transpose[v] * A * v
double Quadratic_Form_3D(double A[3][3],double v[3])
{
 double qform;

 qform = A[0][0]*v[0]*v[0] + A[1][1]*v[1]*v[1] + A[2][2]*v[2]*v[2];
 qform += 2.0*(A[0][1]*v[0]*v[1] + A[0][2]*v[0]*v[2] + A[1][2]*v[1]*v[2]);
 return(qform);

} /* end of Quadratic_Form_3D() */



double Corr_Coeff_Bwn_Two_GAUSS3D_Arrays(int NgaussA,GAUSS3D *gAarray,int NgaussB,GAUSS3D *gBarray)
{
 double ovAB,ovAA,ovBB,CC;

 ovAB =  Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussA,gAarray,NgaussB,gBarray);
 ovAA =  Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussA,gAarray,NgaussA,gAarray);
 ovBB =  Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussB,gBarray,NgaussB,gBarray);

 if ((ovAA>0.0)&&(ovBB>0.0))
 {
   CC = ovAB/sqrt(ovAA*ovBB);
 }
 else CC = 0.0;
 /* printf("ovAB %e ovAA %e ovBB %e CC %lf %e\n",ovAB,ovAA,ovBB,CC); */
 return(CC);
} /* end of Corr_Coeff_Bwn_Two_GAUSS3D_Arrays() */



// Kawabata's version
double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(int NgaussA,GAUSS3D *gAarray,int NgaussB,GAUSS3D *gBarray)
{
 int i,j;
 double ov, ov_all;
 ov_all = 0.0;
 for (i=0;i<NgaussA;++i)
 {
  for (j=0;j<NgaussB;++j)
  {
    ov = Overlap_Integral_Bwn_Two_GAUSS3Ds(&(gAarray[i]),&(gBarray[j]));
    ov_all += gAarray[i].Weight * gBarray[j].Weight * ov;
  }
 }
 return(ov_all);
} /* end of Overlap_Integral_Bwn_GAUSS3D_Arrays() */

//// Mon's version (with bounding boxes)
//// B-gaussians have bounding boxes
//double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(int NgaussA,GAUSS3D *gAarray,int NgaussB,GAUSS3D *gBarray)
//{
// int i,j;
// double ov, ov_all, dist;
// float posA[3],posB[3],diff[3];
// ov_all = 0.0;
//
// for (i=0;i<NgaussA;i++)
// {
//	 // A-i gaussian position
//	 posA[0] = gAarray[i].M[0];
//	 posA[1] = gAarray[i].M[1];
//	 posA[2] = gAarray[i].M[2];
//	 for (j=0;j<NgaussB;j++)
//	 {
//		 diff[0] = posA[0]-posB[0];
//		 diff[1] = posA[1]-posB[1];
//		 diff[2] = posA[2]-posB[2];
//		 dist = sqrt( diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
//
//		 // Bounding box checking
//		 if( gBarray[j].M[0] )
//			 ov = Overlap_Integral_Bwn_Two_GAUSS3Ds(&(gAarray[i]),&(gBarray[j]));
//		 ov_all += gAarray[i].Weight * gBarray[j].Weight * ov;
//	 }
// }
//
// return(ov_all);
//} /* end of Overlap_Integral_Bwn_GAUSS3D_Arrays() */


double Overlap_Integral_Bwn_Two_GAUSS3Ds(GAUSS3D *gA,GAUSS3D *gB)
{
 int i;
 double K[3][3],invK[3][3],DetK,ov;
 double nm[3];
 double qform;
 /*
   Overlap = 1/[(2pi)**3/2 * sqrt(|P+Q|)] *
             exp[-1/2 * tr_(n-m) inv(P+Q) (n-m)]
            = 1/(2pi)**3/2 /sqrt(|K|)*exp[-1/2 * (n-m) invK (n-m)]
    where P = gA->SigM, Q = gB->SigM, m = gA->M, n = gB->M.
     K   = P + Q
  invK   = inv(P + Q)
 */
 Add_Matrix3D(K,gA->SigM,gB->SigM);
 Cal_Inverse_Matrix3D_by_Cramer_Rule(invK,K,&DetK);
 for (i=0;i<3;++i)
	 nm[i] =  gB->M[i] - gA->M[i];
 qform =  Quadratic_Form_3D(invK,nm);
 /* printf("nm %lf %lf %lf qform %lf\n",nm[0],nm[1],nm[2],qform);  */
 ov = exp(-0.5*qform) * one_over_2pi_32 / sqrt(DetK);
 return(ov);
} /* end of double Overlap_Integral_Bwn_Two_GAUSS3Ds() */


// Computing maximum size of gaussians to define their bounding boxes
// (Using covariance matrix eigenvectors/values)
// factor --> bounding box scaling (default = 1.0)
void gaussian_sizes(GAUSS3D *gauss,int num_gauss, float factor)
{
	bool debug=true;
	int i,m,n;
	double eigen_vec[3][3];
	double eigenval[3];

	for(i=0; i<num_gauss; i++)
	{
		// Diagonalize a 3x3 matrix & sort eigenval by size
		diagonalize_symmetric(gauss[i].SigM,gauss[i].evec,eigenval);
		for(n=0; n<3; n++) // x,y,z
			gauss[i].size[n] = 0.0; // max sizes array initialization
		for(m=0; m<3; m++) // eigenvector
		{
			if(debug)
				printf("gauss= %d  evec= %d --> %f %f %f\n",i,m,gauss[i].evec[m][0],gauss[i].evec[m][1],gauss[i].evec[m][2]);
			for(n=0; n<3; n++) // x,y,z
			{
				gauss[i].evec[m][n] *= sqrt(eigenval[m]);
				if(fabs(gauss[i].size[n]) < fabs(gauss[i].evec[m][n]))
					gauss[i].size[n] = fabs( gauss[i].evec[m][n] );
			}
			if(debug)
				printf("gauss= %d  evec= %d --> %f %f %f\n",i,m,gauss[i].evec[m][0],gauss[i].evec[m][1],gauss[i].evec[m][2]);
		}
		for(n=0; n<3; n++) // x,y,z
			gauss[i].size[n] *= factor; // scaling the bounding box
		if(debug)
			printf("gauss= %d  size --> %f %f %f\n",i,gauss[i].size[0],gauss[i].size[1],gauss[i].size[2]);
	}
}
