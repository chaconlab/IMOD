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

#include <iostream>
//#include "libpdb/ResIni.h"
#include "libvolume/include/floatops.h" // ADP's floating point operations library (MUST BE ON TOP !!! ???)
#include "libtools/include/timer.h" // Htimer -structures
#include "libnma/include/libnma_time.h" // Real-timer (Santi's)
//#include "libfrm/Correlation.h" // ADP_EM needed library
#include "libnma/include/nma.h" // Mon's NMA (Cartesian+Dihedral) library
#include "libnma/include/libnma_misc.h" // Mon's NMA Internal Coordinates related library
//#include "libnma/libnma_io.h" // Mon's NMA Input-Output library
//#include "libnma/libnma_cg.h" // Mon's NMA Coarse-Graining library
//#include "libnma/libnma_deriv.h" // Mon's NMA Derivatives library
//#include "libnma/libnma_kinetic.h" // Mon's NMA Kinetic Energy Matrix library
//#include "libnma/libnma_hessian.h" // Mon's NMA Hessian Matrix library
//#include "libnma/libnma_diag.h" // Mon's LAPACK Matrix Diagonalization calls
//#include "libnma/libnma_move.h" // Mon's NMA IC-Motion library

//#include "libtools/Mersenne.h" // Mersenne Twister Random Number generator

/*= DEFINITIONS ================================================================================*/


// #define FILE_NAME 101
// #define ZERO -99999.0
#ifndef ZERO
#define ZERO 0.0
#endif

/*= GLOBAL VARIABLES ===========================================================================*/
typedef struct mode
{
  double *mode;
  double energy;
  double prob;
  double weight;
  int n;
} mode;

typedef struct Mapvox // used by dffe_scoreX()
{
  float diff;
  int i,j,z; // voxel indices
  float pi,pj,pz; // voxel coordinates
} Mapvox;

typedef struct Mapvox2 // used by dffe_scoreX()
{
  int diff;
  int i,j,z; // voxel indices
  float pi,pj,pz; // voxel coordinates
} Mapvox2;

typedef struct Neighbor // used by dffe_scoreX()
{
  int offset; // offset
  float dist; // distance
} Neighbor;

typedef struct tree
{
  tree *left,*right,*parent;
  bool available;
  float sumr,suml;
  int index;
} tree;


/*= FUNCTIONS ==================================================================================*/

// ***********************************************************************************
// Docks First volume to the Second (Target, Reference) and outputs the transformation
// ***********************************************************************************
void docker(vlVolume *vol_lo, vlVolume *vol_hi, float *rotation, float *mov, float *traslation, int bw, int rho_max);
// Projects PDB to VOLUME & Filters it apropriately (Pablo's emt_pdb2vol code)
void pdb2map(Macromolecule *in_pdb, vlVolume **volume, float resolution, float grid_size, int filter_method, bool model_3BB2R=true);

// INVERSE PARABOLIC INTERPOLATION
// Returns the expected minimum displacement coordinate (abscissa, delta_x)
// Minimum = orig.(fx) - parab_interpol()
double parab_iterpol(double fx, double fa, double fb, double offset);
// Rotates and traslates a Macromolecule, checking whether its necessary! (default --> full rigid body rotation)
void rotrans(Macromolecule *mol, float x, float y, float z, float a, float b, float c, int n_chain=-1);

// Implementation of ecuation (3) in Kovacs, Cavasotto and Abagyan (2005)
// Creates an ordered list with the normal mode "relevances" and indexes
void relevance( Macromolecule *mol, int loop_i, int loop_f, tri *props, double *hess_matrix, int size, double angle, float **table );

// Weigted Random Sampling - WithOut Replacement (WRS-WOR)
// Based on the Wong & Easton (1980) 's method from Olken & Rotem (1995)
// (allocates memory to store the selected modes array, if *modes==NULL)
void select_mode_WOR(float *profile, int size, int num, mode **modes, float cutoff);
// It selects "num" modes according to a given probability profile,
// and it sets their weights too!
// "Without replacement" was done DELETING each selected item from list.
// (allocates memory to store the selected modes array, if *modes==NULL)
// profile --> (float *)
// eigval --> (double *) Eigenvalues array for SCV weighting...
void select_mode_prob(float *profile, double *eigval, int size, int num, mode **modes, float cutoff, bool randweight=true);
// It selects "num" modes according to a given probability profile.
// (allocates memory to store the selected modes array, *modes==NULL)
void select_mode_prob2(float *profile, int size, int num, mode **modes, float cutoff);

// It selects the "num" biggest "profile value" modes.
// (allocates memory to store the selected modes array, if *modes==NULL)
void select_mode_big(float *profile, int size, int num, mode **modes);

// It Merges "nev" modes, selected in "nm_props" (FAST)
// (warning: not allocates memory!)
void merge_modes(double *uu, double *hess_matrix, int size, mode *nm_props, int nev);
// It Merges "nev" modes, selected in "nm_props" (FAST)
// (warning: not allocates memory!)
// type = 0 --> phi,psi,...
// type = 1 --> phi,psi,chi,...
// type = 2 --> phi,chi,psi...
//void merge_modes(Macromolecule *mol, NM *uu, tri *props, double *hess_matrix, int size, mode *nm_props, int nev, int type = 1);

// Randomly weights the already present NM weight.
void rand_mode_weight(mode *modes, int num);

// It selects a peak according to a given probability profile.
int select_index(double *profile, int size, float cutoff);

//// Computes dihedral angle diferences between two macromolecules
//// The result will be stored inside "NM" struct
//void dihedral_diff(Macromolecule *mol1, Macromolecule *mol2, NM *uu);

// Sorts indexes from a profile based on its profile values (big-->small)
// (allocates memory to store the selected modes array, if *indexes==NULL)
void sort_profile(float *profile, int **indexes, int size);

// Overwrites the "mol_ref" atomic coordinates over "mol" ones.
void replace_pos(Macromolecule *mol, Macromolecule *mol_ref);

// Projects a Macromolecule into a volume object without unnecessary memory allocations.
// The convolution Kernel and Two appropiately paded volumes must be supplied (vol & dummy).
// "dim_vox" is only needed by "Gausian filter", otherwise it will be ignored
// nthreads --> Number of threads used in parallel. (=0 for non-parallel)
#ifdef USE_PTHREAD // Enables PThread parallel routines
void project_pdb(Macromolecule *mol,vlVolume **pvol,vlVolume **pdummy,float *kernel,int dim_vox,int method, float res, bool fast=false, int nthreads=0, convoluteK_data *threads_data=NULL);
#else
void project_pdb(Macromolecule *mol,vlVolume **pvol,vlVolume **pdummy,float *kernel,int dim_vox,int method, float res, bool fast=false);
#endif

// Accepts or rejects a candidate based on Simulated Annealing
bool anneal( double actual, double next, double temp );

// Projects a Map over a Macromolecule, and outputs a profile with the densities
// within a cube from each residue atom. (The profile has N-residues elements)
// radius = sphere radius in Amstrongs (defines neighborhood)
// size = kernel size (to save time)
// Profile memory is allocated if *p_profile=NULL.
void project_map2pdb(vlVolume *map, Macromolecule *mol,double **p_profile,float radius,int size);

// Smooths a profile, allocating memory for the new profile.
// method = 0 --> Windowed average
// method = 1 --> Gaussian Kernel
double *smooth_profile(double *profile, int size, int kernel, int method);

// Normalizes a profile, (overwriting it!)
// New values will range from "min" to "max" (lineal transformation)
void norm_profile(double *profile, int size, double min_new, double max_new);
void norm_profile(float *profile, int size, float min_new, float max_new);

// Computes normalized cross-correlation (from two un-normalized vectors)
double corr_vector( double *a, double *b, int size );

// B-factor profile computation function
// Cartesian eigenvectors needed.
// It does not allocate "profile" memory. "nm"-->[0,#modes-1]
// ("CA" atom only!)
void bfact_profile(Macromolecule *mol, double *evec, double *profile, int nm);

// Computes "sigma" from a vector (array)
float sig_vector(float *v, int size);

// Creates, whether necessary, an array with "num_atoms" elemets
// copied from an array with "num_res) elements (the same values will be
// set for the same residue atoms)
// "p_out" = NULL --> will lead to memory allocation
void array_res2atom(Macromolecule *mol, double *array, double **p_out);

// Normalizes making that the total voxel summation = "max"
void norm_vol(vlVolume *vol,float max,float sig=0);

// Mon made (6/10/2008)
// Difference map Energy
// map = model map (already allocated)
// target = target map (already allocated)
// offsets = array with the voxel offsets with any chance to be >0 (voxels inside target map)
// n_max = number of elements in "offsets" array (# of voxels that could be >0)
// ("diff" dffe_method)
float dffe_diff2(vlVolume *map, vlVolume *target, int *offsets, int n_max);

// Mon made (25/11/2008)
// UNDER CONSTRUCTION
// It checks both maps maximum dimensions, and makes a new map to where both maps fit.
vlVolume *fitsize(vlVolume *map, vlVolume *map2);

// Pose refinement of a Macromolecule into a Map using: INVERSE PARABOLIC INTERPOLATION
// Based on: Numerical Recipes, pag 402. 10.2 Parabolic Interpolation (1D x 6)
// Returns: correlation --> if pose improvement is significative
//          -1.0 --> if pose improvement is neglected (due to "refine_tol")
// Uses: project_pdb() and correlation_frame()
float pose_ref(
  Macromolecule **p_mol, // p_mol --> Scoring and moving Macromolecule (the same for both)
  vlVolume **p_volF, // Map successfully moved and updated (if return > 0)
  vlVolume **p_volR, // Reference map (untouched)
  vlVolume **p_vol_dummy, // Map buffer
  double corr0, // Initial correlation (fx - point, it saves some time)
  int *chain_order, // array with non-repeated chain indices for multi-chain pose refinement
  int n_chain, // number of chains (saves some time...)
  float *kernel,int dim_vox,int filter,float model_thr,vlDim pad,float resolution, // filtration related (1)
//  bool corr_nozero,float avgR,float sigR, // filtration related (2)
  bool corr_nozero,double avgR,double sigR, // filtration related (2)
  float shake_mov, float shake_rot, // step-size for translation and rotation
  float refine_tol, // significative-ness threshold
  Macromolecule *(*fptr)(Macromolecule *)=NULL, // (OPTIONAL)
	// Then, "mol" will be used for moveing and "mol2" & "fptr(mol)" for scoring.
	// fptr = Macromolecule *(*fptr)(Macromolecule *)  (this is mandatory in: CA-model)
  Macromolecule **p_mol2=NULL // (OPTIONAL) p_mol2 --> Initial Scoring Macromolecule
	// (Note that "mol2" will be a reference to an atoms sub-set in "mol")
#ifdef USE_PTHREAD // Enables PThread parallel routines
  , int nthreads=0 // Number of threads used in parallelization
#endif
  );

//// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
//// allocating profile memory if (*p_profile==NULL).
//// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
//// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
//// Comparison for Flexible Proteins and Predicted Protein Structures".
//// Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
//void gaussian_weight(Macromolecule *mol, Macromolecule *mol2, double **p_profile, double c = 5.0);
//
//// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
//// allocating profile memory if (*p_profile==NULL).
//// (17/7/2012) --> Taking into account sequence alignment (at atomic level) for different size molecules.
////                 *The Weights profile its relative to the first molecule (mol)
//// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
//// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
////                   Comparison for Flexible Proteins and Predicted Protein Structures".
////                   Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
//void gaussian_weight(Macromolecule *mol, Macromolecule *mol2, double **p_profile, bool *mask1, bool *mask2, double c = 5.0);
//
