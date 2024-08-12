/*************************************************************************
 *                   libnmafit's HEADER: nmafit.h                        *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   libnmafit's main header.                                            *
 *   (defines, common definitions, data-types, data-structures, externs) *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

//#include <iostream>
#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
#include <libtools/include/timer.h> // Htimer -structures
//#include "libvolume/floatops.h" // ADP's floating point operations library (MUST BE ON TOP !!! ???)
//#include "libnma/libnma_misc.h" // Mon's NMA Internal Coordinates related library
//#include "libnma/libnma_io.h" // Mon's NMA Input-Output library
//#include "libnma/libnma_cg.h" // Mon's NMA Coarse-Graining library
//#include "libnma/libnma_deriv.h" // Mon's NMA Derivatives library
//#include "libnma/libnma_kinetic.h" // Mon's NMA Kinetic Energy Matrix library
//#include "libnma/libnma_hessian.h" // Mon's NMA Hessian Matrix library
//#include "libnma/libnma_diag.h" // Mon's LAPACK Matrix Diagonalization calls
//#include "libnma/libnma_move.h" // Mon's NMA IC-Motion library

#define one_over_2pi_32  0.0634936359342410

struct GAUSS3D
{
 int num;
 double M[3];         /* Mean Position */
 double SigM[3][3];   /* Covariance Matrix (Sigma) */
 double iSigM[3][3];  /* Inverce of Sigma Matrix   */
 double det;          /* Determinant of SigM       */
 double Cons;         /* 1/{2pi)**3/2 * sqrt(det)} */
 double Weight;       /* Weight for Gaussian Mixture */
 double evec[3][3]; // SigM's eigenvectors
 double size[3];    // Gaussian x,y,z size (from eigenvectors)
};

// LIBGAUSSCORR.CPP
void viewmatrix3D(double matrix[3][3],char *name);
void viewgauss(GAUSS3D *GaussArray, int ngauss);
GAUSS3D *pdb2gauss(Macromolecule *mol,double sd);
void pdb2gauss_update(Macromolecule *mol, GAUSS3D *molgauss);
void Cal_Inverse_Matrix3D_by_Cramer_Rule(double InvA[3][3],double A[3][3],double *Det);
void Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(double InvA[3][3],double A[3][3],double *Det);
double Quadratic_Form_3D(double A[3][3],double v[3]);
void Add_Matrix3D_Symmetric(double C[3][3],double A[3][3],double B[3][3]);
void Add_Matrix3D(double C[3][3],double A[3][3],double B[3][3]);
double Overlap_Integral_Bwn_GAUSS_MOLECULEs(GAUSS3D *Ga,GAUSS3D *Gb,int nga,int ngb,double factor=0.0);
// From GaussOI.c
double Corr_Coeff_Bwn_Two_GAUSS3D_Arrays(int NgaussA,GAUSS3D *gAarray,int NgaussB,GAUSS3D *gBarray);
double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(int NgaussA,GAUSS3D *gAarray,int NgaussB,GAUSS3D *gBarray);
double Overlap_Integral_Bwn_Two_GAUSS3Ds(GAUSS3D *gA,GAUSS3D *gB);

// LIBGAUSSIO.CPP
void Write_Gaussian3D_File(char *fname,int Ngauss,GAUSS3D *Garray,char  chain,char *comment);
char Read_Gaussian3D_File(char *fname,int *Ngauss,GAUSS3D *Garray);
char readGAUSS3D(char *fname,int *Ngauss,GAUSS3D **p_Garray);
double Corr_Coeff_Bwn_Two_GAUSS3D_Arrays();
double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays();
double Overlap_Integral_Bwn_Two_GAUSS3Ds();
void Get_Part_Of_Line(char *part,char *line,int s,int e);
void Split_to_Words(char *str,char splsym,int *Nword,int Wsta[],int Wend[],int Nwordmax);
char *Get_Date_String();

// Computing maximum size of gaussians to define their bounding boxes
// (Using covariance matrix eigenvectors/values)
void gaussian_sizes(GAUSS3D *gauss,int num_gauss,float factor=1.0);
// Writes a VMD file with vectors representing the 3D-Gaussians in "gauss"
// Warning: the eigenvectors should be already computed using "gaussian_sizes()"
void write_gauss_evecs(GAUSS3D *gauss,int num_gauss,char *file);
// Writes a VMD file with vectors representing the 3D-Gaussians in "gauss"
// Warning: the eigenvectors should be already computed using "gaussian_sizes()"
void write_gauss_sizes(GAUSS3D *gauss,int num_gauss,char *file);
