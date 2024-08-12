/************************************************************************
 *                           NMAVIEW                                     *
 ************************************************************************
 * This program is part of iMOD: http://chaconlab.org/imod/index.html   *
 * (c) Jose Ramon Lopez-Blanco, Jose Ignacio Garzon and Pablo Chacon.   *
 * IQFR-CSIC's Structural Bioinformatics Group (2004-2011).             *
 ************************************************************************
 *                                                                       *
 *   Tool for appling the CCS mode displacements to pdb models           *
 *   also draws arrows in VMD to represent the displacements.            *
 *                                                                       *
 *************************************************************************
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <stdio.h> // needed by some compilers
#include "cmdl/CmdLine.h"
#include "libnma/include/libnma_io.h" // Mon's NMA Input-Output library
#include "libnma/include/libnma_affine.h" // Mon's NMA Affine-models library

char version[]="1.08"; // version code
float move_delta = 1.0;
float percent_value;
float factor = 1.0;
bool delta_switch=false;
bool novmd_switch = false;
bool numNM_switch=false;
bool switch_max=true;
bool plot_sphere = false; // Plot spheres
bool plot_valid = false; // Plot valid points
bool plot_paths = false; // Plot valid paths before choosing the best one
// bool plot_best = true; // Plot best arrows
bool plot_axis = false; // Plot best axis
bool plot_affine = false; // Plot affine model arrows
bool pdb_clusters = false; // Output clusters as the b-factors of pdbs
bool server = false; // =true, to enable "server mode" (maximize the default automated options)
bool NClusters_auto = false; // =true enable automatic selection of number of clusters
int NM_considered;
int operation = 1; // operation method
int arrow_level=0; // 0,1,2,3 --> atom, residue, segment, chain (respectively)
int nresidues=1; // Number of residues to average modes components (used only by "--residues [nresidues]" option)
float kthr = 0.2;
float kthr2 = 0.02;
float radius=0.15; // cylinder radius
float asize=25; // size of the axis cylinder

// Clustering related global variables
int requested = 10; // Number of requested clusters
float adjacency = 6.0; // Adjacency cutoff (in Angstroms)
int extend_size = 5; // Number of atoms threshold to consider cluster extension (clusters with less atoms will be extended). (Default=5)
float extend_cutoff = 10; // Distance (in Angstroms) to extend small clusters. (Default=10A)

// Arrow-illustration related global variables
int format = 2; // Output format for Arrows: 1=VMD, 2=JMol`s PMESH   (just arrows, the remaining visualization items will be only in VMD format)
int nsphere = 300; // Number of sphere points
double rsphere = 50; // Radius of the sphere (for visualization only)

// Sampling valid points variables
float integration_step = 500;
float delta_sampling = 2.0; // Valid points sampling increment (in Angstroms)
float far2path_thr = 12.0; // Distance threshold to determine when a trajectory goes too far from the cluster.
float close2path_thr = 7.0; // Distance to detect collision between path points and both current sampled points and any atom of the whole molecule.
float maxLengthFactor = M_PI;
float surf_dist = 5.0; // Distance from the "surface"
float close2vector_thr = 10.0; // Threshold to determine which points are close to a given vector.
float close2clusters_thr = 5.0; // threshold to determine which points are close to any different cluster atom.

// Arrow Drawing variables
int arrow_number = 1; // Number of arrows
float delta_path = 1.0; // Arrow and Path increment (in Angstroms)
float ribbon_width = 1.0; // Ribbon width in Angstroms
float arrow_width = 2.0; // Arrow width in Angstroms
float tip_width = 4.0; // Arrow's Tip width (in Angstroms)
float tip_length = 6.0; // Arrow's Tip length (in Angstroms)
float arrow_aspect = 0.7; // Arrow aspect ratio (non_arrow_tip / arrow_length) <-- (It will be redefined later if "tip_length" > 0)
float arrow_to_path = 0.8; // Arrow to path ratio
float min_aspect = 0.4; // Minimal aspect ratio when arrow length is smaller than tip length
float arrows_apart = 0.75; // Separation (in Angstroms) between arrows and paths
bool path_middle_lines = false; // Enable or disable path middle lines

// Input variables
char file_pdb[FILE_NAME],file_nma[FILE_NAME],file_out[FILE_NAME],file_fix[FILE_NAME],name[FILE_NAME],text[FILE_NAME]; // file names
char colortext[20]; // color name string

using namespace TCLAP;

// Input parser
void parseOptions(int argc, char** argv);

// Selecting valid points for a given cluster (those not occluded by any neighboring atom)
// First, the homogeneously-distributed spherical unit vectors are centered at current cluster CoM.
// Next, the positions of the atoms are projected into every unit vector. Only those atoms close to the vectors are considered.
// Then, the maximum projection found for each vector determines the potential valid point.
// Finally, any potential valid point non-occluded by any atom is separated some distance outwards and stored as valid sampling point.
//
//  cluster --> Current cluster
//  ncluster --> Number of cluster atoms
//  c --> Cluster index
//  com --> Cluster center of mass
//  p_samples --> Pointer to the valid points array (OUTPUT)
//  p_nsamples --> Pointer to the number of valid points (OUTPUT)
//  sp --> Spherical distribution of points
//  nsphere --> Number of spherical distribution of points
//  coords --> Full molecule coordinates
//  coords2 --> Full molecule coordinates to be centered at cluster center
//  atom_owner --> Atomic cluster ownership array
//  num_atoms --> Macromolecule number of atoms
//  surf_dist = 5.0; // Distance from the "surface"
//  close2vector_thr = 10.0; // Threshold to determine which points are close to a given vector.
//  close2clusters_thr = 5.0; // threshold to determine which points are close to any different cluster atom.
//  f_name --> Output filename for the centered spherical distribution of points (OPTIONAL)
void getValid(int *cluster, int ncluster, int c, double *com, float **p_samples, int *p_nsamples, float *sp, int nsphere, float *coords, float *coords2, int *atom_owner, int num_atoms, float surf_dist, float close2vector_thr, float close2clusters_thr, char *f_name=NULL);

// Compute the length-per-unit given one affine model and one samples set.
//  lengthPerUnit --> Array with the lengths per unit (it should have been already allocated)
//  samples --> Valid sample points set
//  nsamples --> Number of valid sample points
//  affine --> Affine model
//  integration_step --> Affine model integration step
void getLengthPerUnit(float *lengthPerUnit, float *samples, int nsamples, double *affine, float integration_step);

// Get the best valid point (the one that leads to the best path)
void getBest(double *affine, float *samples, int nsamples, int *cluster, int ncluster, int c, float *comf, float *lengthPerUnit, float *coords, float far2path_thr, float close2path_thr, float delta_sampling, int *atom_owner, int num_atoms, float integration_step, int *p_i_best, float *p_jmax_best, float *p_jmin_best, float *p_pathLength_best, char *paths_name=NULL);

// Recursive definition of determinate using expansion by minors.
double Determinant(double **a,int n);

// Draws arrow/s and path from affine model
//
// format --> 1: VMD, 2: JMol
// sample --> Best valid sample coordinates pointer
// affine --> Affine model matrix (4x4)
// jmin_best --> Initial path point relative to affine model)
// jmax_best --> Final path point (relative to affine model)
// pathLength_best --> Best path length
// lengthPerUnit --> Length per unit for best valid sample point
// comf --> Cluster's Center of Mass (weighted by vector field)
// delta_path --> Arrow sampling (in Angstroms)
// integration_step --> Affine model integration step
// f_name --> Best arrow & path file name (without extension, it will be aded according to "format")// arrow_number --> Number of arrows
// arrow_length --> Arrow length
// ribbon_width --> Ribbon width in Angstroms
// arrow_width --> Arrow width in Angstroms
// tip_width --> Arrow's Tip width (in Angstroms)
// tip_length --> Arrow's Tip length (in Angstroms)
// arrow_aspect --> // Arrow aspect ratio (non_arrow_tip / arrow_length) <-- (It will be redefined later if "tip_length" > 0)
// arrow_to_path --> Arrow to path ratio
// min_aspect --> Minimal aspect ratio when arrow length is smaller than tip length
// arrows_apart --> Separation (in Angstroms) between arrows and paths
// path_middle_lines --> Enable or disable path middle lines
void drawArrow(int format, float *sample, double *affine, float jmin_best, float jmax_best, float pathLength_best, float lengthPerUnit, float *comf,
		float delta_path, float integration_step,	char *f_name, int arrow_number, float arrow_length, float ribbon_width, float arrow_width, float tip_width,
		float tip_length, float arrow_aspect, float arrow_to_path, float min_aspect, float arrows_apart, bool path_middle_lines);

void getMoreOwners(Macromolecule *mol, int **p_res_owner, char **p_chain_owner);

void writeCluster(int *cluster, int natoms, int *res_owner, char *chain_owner, char *text);

// Source code of simple quick sort implementation using array ascending order in c programming language
// Taken from: http://www.cquestions.com/2008/01/c-program-for-quick-sort.html
void quicksort(int *x, int first, int last);


///////////////////////////////////////////////////
int main(int argc, char **argv)
{
	bool debug = false;
	float max_delta,amp,max_amp=0,correction,amp_thr=0;
	float pointx,pointy,pointz;
	float pointx2,pointy2,pointz2;
	int atoms,i,j,k,max_amp_atom,n_nm,offset;
	int num_atoms, num_atomsB, num_res;
	trs *nm;
	double *evec,*eval;
	FILE *Fout;
	Tcoor pos,pos2;
	file_fix[0]='\n'; // Signaling not introduced any file

	printf("imodview>\nimodview> Welcome to the Normal Modes and Springs visualization Tool  v%s\nimodview>\n",version);

	parseOptions(argc,argv);

	printf("imodview>\nimodview> PDB file %s\nimodview> Nm file %s\n",file_pdb,file_nma);
	printf("imodview> Move %f\n",move_delta);

	// Initialize aminoacids
	init_aminoacids();

	// Read pdb
	Macromolecule *pdb = new Macromolecule( file_pdb );
	pdb->readPDB( file_pdb );
	pdb->info(stdout);
	num_atoms = pdb->get_num_atoms();
	num_res = pdb->get_num_fragments();

	// Open output VMD file (common)
	if ((Fout = fopen(file_out, "w")) == NULL)
	{
		printf("imodview> Intput mode file %s error!! \n",file_out);
		exit(1);
	}

	if(format == 1) // if VMD
	{
		fprintf(Fout,"molecule new\n");
		fprintf(Fout,"display resetview\n");
		fprintf(Fout,"draw color %s\n", colortext);
	}

	int num_nms1=0,num_atoms1=0,num_comps1=0;
	int nipa,cont_springs=0;

	switch(operation)
	{
	case 1: // DRAW ARROWS
		// Reading "ptraj" NMs format
		trs *nms1;
		printf("Reading NMs according to ""ptraj"" format\n");
		nm = read_ptraj(file_nma,num_atoms,&num_nms1);
		n_nm = num_nms1;
		// Some input checking...
		if(NM_considered >= num_nms1 || NM_considered < 0)
		{
			printf("\nPlease enter a avalid NM index (from 1 to %d).\n\n",num_nms1);
			exit(1);
		}

		// Maximum amplitude check
		for(i=0;i<num_atoms;i++)
		{
			amp = sqrt( nm[i + NM_considered * num_atoms].x * nm[i + NM_considered * num_atoms].x
					+ nm[i + NM_considered * num_atoms].y * nm[i + NM_considered * num_atoms].y
					+ nm[i + NM_considered * num_atoms].z * nm[i + NM_considered * num_atoms].z ); // 3D amplitude
			if ( amp > max_amp ) // searches the maximum amplitude to normalize
			{
				max_amp = amp;
				max_amp_atom = i;
			}
		}
		printf("imodview> Maximum NM amplitude (3D): %lf\t for atom %d\n",max_amp,max_amp_atom);

		// Sets the minimum amplitude threshold to be showed
		amp_thr = max_amp * ( percent_value / 100 );

		// Normalizes & Corrects the delta factor to reach the desired maximum displacement
		if(switch_max)
			factor = (1.0/max_amp)*move_delta; // now, factor normalizes the maximum displacement to 1A.

		// Open output VMD file (common)
		printf("imodview> Mode selected: %d\n",NM_considered+1);
		printf("imodview> Showing displacements with 3D amplitude bigger than %.1f %%.\n",percent_value);
		printf("imodview> (i.e. %.5f).\n",amp_thr);
		printf("imodview> Save VMD/JMol arrows in: %s\n", file_out);

		// SERVER MODE: Automatic overriding and/or selection of parameters.
		if(server)
		{
			int too_many_res = 500; // "too many residues" threshold
			arrow_level = 1; // residue-wise arrows computation
			printf("imodview> Server mode ON: Overriding arrow_level = %d\n", arrow_level);
			if(num_res > too_many_res) // if too many residues...
			{
				nresidues = (int) ceil((float)num_res/too_many_res); // override the number of residues for arrow averaging...
				printf("imodview> Server mode ON: Overriding number of residues for arrow averaging = %d\n", nresidues);
			}
		}

		// Iterador para recorrer atomos
		pdbIter *iter_frag,*iter_seg,*iter_ch,*iter;
		Fragment *frag;
		Segment *seg;
		Chain *ch;
		double com[3],nmav[3]; // Center of mass, and NM AVerage.
		float dummy[3],dummy2[3],dummy3[3];
		offset = 0;
		switch(arrow_level)
		{
		case 0: // ATOM
			iter = new pdbIter( pdb );
			for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screen atoms
			{
				i = iter->pos_atom; // it's shorter... (I'm a lazy guy)

				amp = sqrt(  nm[i + NM_considered * num_atoms].x * nm[i + NM_considered * num_atoms].x
						+ nm[i + NM_considered * num_atoms].y * nm[i + NM_considered * num_atoms].y
						+ nm[i + NM_considered * num_atoms].z * nm[i + NM_considered * num_atoms].z ); // 3D amplitude
				if ( amp > amp_thr ) // If the atom moves significantly (user defined)
				{
					// get initial position
					(iter->get_atom())->getPosition(pos);

					// get final point
					pointx2=pos[0]+nm[i + NM_considered * num_atoms].x*factor;
					pointy2=pos[1]+nm[i + NM_considered * num_atoms].y*factor;
					pointz2=pos[2]+nm[i + NM_considered * num_atoms].z*factor;

					if(format == 1) // VMD
					{
						// Drawing cylinder
						pointx=pos[0]+nm[i + NM_considered * num_atoms].x*factor*0.7;
						pointy=pos[1]+nm[i + NM_considered * num_atoms].y*factor*0.7;
						pointz=pos[2]+nm[i + NM_considered * num_atoms].z*factor*0.7;
						fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
								pos[0],pos[1],pos[2],pointx,pointy,pointz,radius);
						// Drawing cone (now it's an arrow!)
						fprintf(Fout,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
								pointx,pointy,pointz,pointx2,pointy2,pointz2,radius*3);
					}
					else // JMol
					{
						fprintf(Fout,"draw ar%d arrow FIXED {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} color %s width %.4f\n",iter->pos_atom+1,
								pos[0],pos[1],pos[2],pointx2,pointy2,pointz2,colortext,radius*3);
					}
				}
			}
			break;
		case 1: // FRAGMENT
			iter_frag = new pdbIter( pdb );
			i=0;
			j=0;
			com[0] = com[1] = com[2] = 0.0;
			nmav[0] = nmav[1] = nmav[2] = 0.0;
			for(iter_frag->pos_fragment=0; !iter_frag->gend_fragment(); iter_frag->next_fragment())
			{
				frag = iter_frag->get_fragment();
				iter = new pdbIter( frag );

				// Iters each fragment
				for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screen atoms
				{
					offset = i + NM_considered * num_atoms;
					// get mode (atom vector)
					nmav[0] += nm[offset].x;
					nmav[1] += nm[offset].y;
					nmav[2] += nm[offset].z;
					// get position
					(iter->get_atom())->getPosition(pos);
					com[0] += pos[0];
					com[1] += pos[1];
					com[2] += pos[2];

					i++;
					j++;
				}

				if( (iter_frag->pos_fragment + 1) % nresidues == 0 )
				{
					// computing averages...
					nmav[0] /= j;
					nmav[1] /= j;
					nmav[2] /= j;
					com[0] /= j;
					com[1] /= j;
					com[2] /= j;

					// Drawing arrow
					amp = sqrt( nmav[0] * nmav[0] + nmav[1] * nmav[1] + nmav[2] * nmav[2] ); // 3D amplitude
					if ( amp > amp_thr ) // If the fragment moves significatively (user defined)
					{
						// get final point
						pointx2=com[0]+nmav[0]*factor;
						pointy2=com[1]+nmav[1]*factor;
						pointz2=com[2]+nmav[2]*factor;

						if(format == 1) // VMD
						{
							// Drawing cylinder
							pointx=com[0]+nmav[0]*factor*0.7;
							pointy=com[1]+nmav[1]*factor*0.7;
							pointz=com[2]+nmav[2]*factor*0.7;
							fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
									com[0],com[1],com[2],pointx,pointy,pointz,radius);

							// Drawing cone (now it's an arrow!)
							fprintf(Fout,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
									pointx,pointy,pointz,pointx2,pointy2,pointz2,radius*3);
						}
						else // JMol
						{
							fprintf(Fout,"draw ar%d arrow FIXED {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} color %s width %.4f\n",iter_frag->pos_fragment+1,
									com[0],com[1],com[2],pointx2,pointy2,pointz2,colortext,radius*3);
						}
					}

					com[0] = com[1] = com[2] = 0.0;
					nmav[0] = nmav[1] = nmav[2] = 0.0;
					j=0; // counts fragment atoms
				}

				delete iter;
			}
			delete iter_frag;
			break;
		case 2: // SEGMENT
			iter_seg = new pdbIter( pdb );
			i=0;
			for(iter_seg->pos_segment=0; !iter_seg->gend_segment(); iter_seg->next_segment())
			{
				seg = iter_seg->get_segment();
				iter = new pdbIter( seg );
				com[0] = com[1] = com[2] = 0.0;
				nmav[0] = nmav[1] = nmav[2] = 0.0;

				// Iters each fragment
				j=0; // counts fragment atoms
				for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screen atoms
				{
					offset = i + NM_considered * num_atoms;
					// get mode (atom vector)
					nmav[0] += nm[offset].x;
					nmav[1] += nm[offset].y;
					nmav[2] += nm[offset].z;
					// get position
					(iter->get_atom())->getPosition(pos);
					com[0] += pos[0];
					com[1] += pos[1];
					com[2] += pos[2];

					i++;
					j++;
				}
				// computing averages...
				nmav[0] /= j;
				nmav[1] /= j;
				nmav[2] /= j;
				com[0] /= j;
				com[1] /= j;
				com[2] /= j;

				// Drawing arrow
				amp = sqrt( nmav[0] * nmav[0] + nmav[1] * nmav[1] + nmav[2] * nmav[2] ); // 3D amplitude
				if ( amp > amp_thr ) // If the fragment moves significatively (user defined)
				{
					// get final point
					pointx2=com[0]+nmav[0]*factor;
					pointy2=com[1]+nmav[1]*factor;
					pointz2=com[2]+nmav[2]*factor;

					if(format == 1) // VMD
					{
						// Drawing cylinder
						pointx=com[0]+nmav[0]*factor*0.7;
						pointy=com[1]+nmav[1]*factor*0.7;
						pointz=com[2]+nmav[2]*factor*0.7;
						fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
								com[0],com[1],com[2],pointx,pointy,pointz,radius);

						// Drawing cone (now it's an arrow!)
						fprintf(Fout,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
								pointx,pointy,pointz,pointx2,pointy2,pointz2,radius*3);
					}
					else // JMol
					{
						fprintf(Fout,"draw ar%d arrow FIXED {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} color %s width %.4f\n",iter_seg->pos_segment+1,
								com[0],com[1],com[2],pointx2,pointy2,pointz2,colortext,radius*3);
					}
				}

				delete iter;
			}
			delete iter_seg;
			break;
		case 3: // CHAIN
			iter_ch = new pdbIter( pdb );
			i=0;
			for(iter_ch->pos_chain=0; !iter_ch->gend_chain(); iter_ch->next_chain())
			{
				ch = iter_ch->get_chain();
				iter = new pdbIter( ch );
				com[0] = com[1] = com[2] = 0.0;
				nmav[0] = nmav[1] = nmav[2] = 0.0;

				// Iters each fragment
				j=0; // counts fragment atoms
				for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screen atoms
				{
					offset = i + NM_considered * num_atoms;
					// get mode (atom vector)
					nmav[0] += nm[offset].x;
					nmav[1] += nm[offset].y;
					nmav[2] += nm[offset].z;
					// get position
					(iter->get_atom())->getPosition(pos);
					com[0] += pos[0];
					com[1] += pos[1];
					com[2] += pos[2];

					i++;
					j++;
				}
				// computing averages...
				nmav[0] /= j;
				nmav[1] /= j;
				nmav[2] /= j;
				com[0] /= j;
				com[1] /= j;
				com[2] /= j;

				// Drawing arrow
				amp = sqrt( nmav[0] * nmav[0] + nmav[1] * nmav[1] + nmav[2] * nmav[2] ); // 3D amplitude
				if ( amp > amp_thr ) // If the fragment moves significatively (user defined)
				{
					// get final point
					pointx2=com[0]+nmav[0]*factor;
					pointy2=com[1]+nmav[1]*factor;
					pointz2=com[2]+nmav[2]*factor;

					if(format == 1) // VMD
					{
						// Drawing cylinder
						pointx=com[0]+nmav[0]*factor*0.7;
						pointy=com[1]+nmav[1]*factor*0.7;
						pointz=com[2]+nmav[2]*factor*0.7;
						fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
								com[0],com[1],com[2],pointx,pointy,pointz,radius);

						// Drawing cone (now it's an arrow!)
						fprintf(Fout,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
								pointx,pointy,pointz,pointx2,pointy2,pointz2,radius*3);
					}
					else // JMol
					{
						fprintf(Fout,"draw ar%d arrow FIXED {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} color %s width %.4f\n",iter_ch->pos_chain+1,
								com[0],com[1],com[2],pointx2,pointy2,pointz2,colortext,radius*3);
					}
				}
				delete iter;
			}
			delete iter_ch;
			break;
		}
		break;

		case 2: // DRAW SPRINGS
			twid *ipas;

			printf("imodview> Reading force constants file: %s",file_nma);
			// Read Force constants file (Kfile) allocating memory
			// Warning: if "coord" not provided, distances ".d" will not be updated!
			read_Kfile(&ipas,&nipa,file_nma);
			printf(" (%d nipas)\n",nipa);

			iter = new pdbIter(pdb);
			for(k=0; k<nipa; k++)
			{
				if(ipas[k].C > kthr && ipas[k].C < kthr2)
				{
					iter->pos_atom=ipas[k].k;
					iter->get_atom()->getPosition(pos); // get i-atom position
					iter->pos_atom=ipas[k].l;
					iter->get_atom()->getPosition(pos2); // get j-atom position
					fprintf(Fout,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %f resolution 20 filled 1\n",
							pos[0],pos[1],pos[2],pos2[0],pos2[1],pos2[2],ipas[k].C*radius);
					cont_springs++; // counts depicted springs
				}
			}

			printf("imodview> Saved %d springs in file: %s\nimodview> END\n",cont_springs,file_out);
			break;

		case 3: // AFFINE MODEL based clustering and arrow visualization
		{
			printf("imodview> AFFINE MODEL based clustering and arrow visualization.\n"
					"imodview>\n"
					"imodview> Implemented and optimized from:\n"
					"imodview> \"Automated Illustration of Molecular Flexibility\"\n"
					"imodview> Aaron Bryden, George Phillips Jr., Michael Gleicher\n"
					"imodview> IEEE Transactions on Visualization and Computer Graphics, 18:1 132-145 (2012)\n"
					"imodview>\n");

			// Reads an entire ptraj file allocating memory
			int check_natoms;
			read_ptraj( file_nma, &evec, &eval, &check_natoms, &num_nms1, &num_comps1);

			if(check_natoms != num_atoms)
			{
				fprintf(stderr,"Different number of atoms between input structure (%d) and eigenvectors file (%d)\nForcing exit!\n",num_atoms,check_natoms);
				exit(1);
			}

			float *coords; // Non-modified atomic coordinates
			float *coords2; // Cluster centered atomic coordinates
			pdb->coordMatrix(&coords); // Get 1D array with atomic coordinates
			pdb->coordMatrix(&coords2); // Get 1D array with atomic coordinates

			printf("imodview> %d modes readed, with %d components, using mode %d\n",num_nms1,num_comps1,NM_considered+1);

			// Some input checking...
			if(NM_considered >= num_nms1 || NM_considered < 0)
			{
				printf("\nPlease enter a avalid NM index (from 1 to %d).\n\n",num_nms1);
				exit(1);
			}
			double *evecx = evec + num_comps1*NM_considered;

			int *rigid_bodies = NULL; // Rigid bodies ownership... (as a function of residue index)
			int nrigid = 0; // Number of rigid bodies detected
			if(file_fix[0] != '\n')
			{
				rigid_bodies = (int *) allocate(sizeof(int) * num_res, "rigid_bodies memory allocation failed!");
				// Reads the "fix" file and outputs a "number of residues sized" array of "ints" with the indices
				// of those "rigid bodies" as a function of the residues. It searches for chunks of continuously fixed
				// sequences of residues and assigns correlative indices to the corresponding "rigid bodies".
				nrigid = read_fixRes( file_fix, pdb, rigid_bodies);
				printf("imodview> Readed %d different rigid-bodies from %s\n",nrigid,file_fix);
			}

			int ***rclusters = NULL;
			int **rnatoms = NULL;
			float **rerrors = NULL;
			printf("imodview> HIERARCHICAL CLUSTERING... requested clusters: %d\n",requested);
			if(requested > num_atoms)
			{
				fprintf(stderr,"Warning: More clusters requested (%d) than atoms available (%d)! Setting it to %d!\n",requested,num_atoms,num_atoms);
				requested = num_atoms-1;
			}
			doHiearchicalExpMatClustering(pdb, coords, requested, evecx, adjacency, &rclusters, &rnatoms, &rerrors, rigid_bodies, nrigid, extend_size, extend_cutoff);

			// COMPUTING CLUSTERS SPEED (To determine the relative arrow length)
			float speed;
			float speed_max;  // maximum cluster speed
			float **rspeeds;
			float *rspeeds_max; // maximum speeds per requested number of clusters
			rspeeds = (float **) allocate(sizeof(float *) * requested, "");
			rspeeds_max = (float *) allocate(sizeof(float) * requested, ""); // Allocating  memory for maximum-speeds array
			for(int i=0; i<requested; i++) // screen requested clusters
			{
				rspeeds[i] = (float *) allocate(sizeof(float) * (i+1), "");  // Allocating  memory for speeds
				speed_max = 0.0;
				for(int j=0; j<=i; j++) // screen clusters
				{
					speed = 0.0; // speed summation for current cluster
					for(int m=0; m<rnatoms[i][j]; m++)
						speed += magnitude(evecx + 3*rclusters[i][j][m],3); // speeds are proportional to eigenvector magnitude
					speed /= rnatoms[i][j]; // average speed per atom
					rspeeds[i][j] = speed;
					if(speed > speed_max)
						speed_max = speed; // store the maximum cluster speed
				}
				rspeeds_max[i] = speed_max; // store the maximum cluster speed per requested clusters
			}

			//	  // Array to store the cluster index of a given atom
			//	  int *atom_owner = (int *) allocate(sizeof(int) * num_atoms, "");
			//	  for(int i=0; i<num_atoms; i++)
			//		  atom_owner[i] = -1; // -1 means unassigned!

			// Setting which cluster any atom belongs to ("atom_owner")
			//	  int r; // requested number of clusters
			//	  int c; // cluster index
			//	  for(r = 0; r < requested; r++) // Screen requested number of clusters
			//		  for(c = 0; c <= r; c++) // Screen its clusters
			//			  for(int i = 0; i<rnatoms[r][c]; i++) // screen cluster's atoms
			//				  atom_owner[ rclusters[r][c][i] ] = c; // Set atom owner (cluster index an atom belongs to)

			// Get Residue-index and Chain-id arrays...
			int *res_owner = NULL;
			char *chain_owner = NULL;
			getMoreOwners(pdb, &res_owner, &chain_owner);
			//			// Show chain and res owner
			//			for(int i=0; i<num_atoms; i++) // screen cluster atoms
			//				fprintf(stderr,"res= %4d  chain= %c\n",res_owner[i], chain_owner[i]);

			// SAVING AND DUMPING CLUSTERING OUTPUT
			float toterr;
			double *atomcluster;
			if(pdb_clusters)
				atomcluster = (double *) allocate(sizeof(double) * num_atoms, ""); // Allocating  memory for atom cluster ownership array

			double *error_elbow, *error_elbow_D1, *error_elbow_D2, *error_elbow_S;;

			error_elbow = (double *) allocate(sizeof(double) * requested, ""); // Allocating  memory for atom cluster ownership array
			error_elbow_D1 = (double *) allocate(sizeof(double) * requested, ""); // Allocating  memory for atom cluster ownership array
			error_elbow_D2 = (double *) allocate(sizeof(double) * requested, ""); // Allocating  memory for atom cluster ownership array
			error_elbow_S = (double *) allocate(sizeof(double) * requested, ""); // Allocating  memory for atom cluster ownership array


			for(int i=0; i<requested; i++)
			{
				// Sorting the atom indices of every cluster
				for(int n=0; n<=i; n++)
					quicksort(rclusters[i][n], 0, rnatoms[i][n]-1);

				// // Show all clusters and their atoms (for a requested number of clusters)
				// showClusters(rclusters[i], i+1, rnatoms[i], rerrors[i]);


				toterr = 0.0;
				for(int n=0; n<=i; n++) {
					toterr += rerrors[i][n] * rnatoms[i][n];
					//  Show speeds and error for each cluster
					// printf("%2d %2d %6d %10.4e (%5.3f) %10.4e\n",i+1,n+1,rnatoms[i][n],rspeeds[i][n],rspeeds[i][n]/rspeeds_max[i],rerrors[i][n]);
				}
				// Show total error for each cluster
				error_elbow[i]=toterr;
			}

			error_elbow_D1[0] = 0;
			for(int i=1; i<requested-1; i++) {
				error_elbow_D1[i] = error_elbow[i]-error_elbow[i+1];
			}
			error_elbow_D2[0] = 0; error_elbow_D2[1] = 0; error_elbow_S[0] = 0;  error_elbow_S[1] = 0;
			double elbow_max = -1.0; int best_elbow;
			for(int i=2; i<requested-1; i++) {
				error_elbow_D2[i] = error_elbow_D1[i]-error_elbow_D2[i+1];
				error_elbow_S[i]=(error_elbow_D2[i]-error_elbow_D1[i+1])/(i/3.0);
				if (error_elbow_S[i] > elbow_max ) { elbow_max =  error_elbow_S[i]; best_elbow =i;}
			}

			FILE *fel;
			if( (fel = fopen("clst_elbow.txt", "w")) == NULL)
			{
				printf("imodview>  clst_elbow.txt writing error!! \n");
				exit(1);
			}
			printf("imodview>  Best cluster size %d (%7.3f) alternatives: ",best_elbow+1, error_elbow_S[best_elbow]);
			double elbow_fact=error_elbow[requested-1]/elbow_max*0.9;
			for(int i=0; i<requested-1; i++) {
				if ((i!=best_elbow) && (error_elbow_S[i] > elbow_max * 0.2))
					printf("%d(%7.3f) -",i+1,error_elbow_S[i]);
				if (error_elbow_S[i] > elbow_max * 0.95)
					fprintf(fel,"%3d %7.3f %7.3f %7.3f %7.3f %7.3f\n", i+1,error_elbow[i], error_elbow_D1[i], error_elbow_D2[i], error_elbow_S[i]*elbow_fact, ceil(error_elbow[0]) );
				else
					fprintf(fel,"%3d %7.3f %7.3f %7.3f %7.3f %7.3f\n", i+1,error_elbow[i], error_elbow_D1[i], error_elbow_D2[i], error_elbow_S[i]*elbow_fact, 0.0 );
			}
			fclose(fel);


			// arc distance
			double P1x,P2x,P1y,P2y, aL, bL, cL, MD_arc, temp_arc;
			int best_arc[3];
			for(int j=0; j<3; j++) {
				P1x=j;
				P1y=error_elbow[j];
				P2x=requested-1;
				P2y=error_elbow[requested-1];
				aL = P2y - P1y;
				bL = -(P2x - P1x);
				cL = P1y * P2x - P1x * P2y;
				// printf("f(x)=(%f*x+%f)/%f\n",  aL, cL, bL);
				MD_arc = 0.0;
				for(int i=j; i<requested; i++) {
					temp_arc= 	fabs((aL * i + bL * error_elbow[i] + cL)) /  (sqrt(aL * aL + bL * bL));
					// printf("%d %f\n", i+1, temp_arc);
					if (temp_arc > MD_arc ) { MD_arc = temp_arc; best_arc[j]=i+1;}
				}
			}

			printf("----->  ");
			for(int j=0; j<3; j++) {
				printf("(%d)->%d ", j+1, best_arc[j]);
			}
			printf("\n");


			int N_start, N_end;

			if (NClusters_auto) {
				N_start = best_elbow;
				N_end = N_start +1;
			} else {
				N_start = 0;
				N_end = requested;
			}
			// ************************
			// SAVING PDBs & TXT
			// ************************

			for(int i=N_start; i<N_end; i++)
			{

				if(pdb_clusters)
				{
					// Write pdb files with cluster information in b-factors
					sprintf(text,"%scluster_n%d_r%d.pdb",name,NM_considered+1,i+1);
					//					for(int n=i; n>=0; n--) // n --> screen clusters
					for(int n=0; n<=i; n++) // n --> screen clusters
						for(int m=0; m<rnatoms[i][n]; m++) // m --> atom index
							atomcluster[rclusters[i][n][m]] = 99 * (double) n / (double)(i+1); // set atom cluster ownership
					pdb->exchange_Pdbfact(atomcluster);
					pdb->writePDB(text); // write pdb with cluster information in B-factors
				}

				// Write plain-text file with cluster residues selection (for Jmol)
				for(int n=0; n<=i; n++) // screen clusters
				{
					sprintf(text,"%scluster_n%d_r%d_c%d.txt",name,NM_considered+1,i+1,n+1);
					writeCluster(rclusters[i][n],rnatoms[i][n],res_owner,chain_owner,text);
				}
			}
			free(res_owner);
			free(chain_owner);
			if(pdb_clusters)
				free(atomcluster);




			// ************************
			// ARROW-related things....
			// ************************

			// Array to store the cluster index of a given atom
			int *atom_owner = (int *) allocate(sizeof(int) * num_atoms, "");
			for(int i=0; i<num_atoms; i++)
				atom_owner[i] = -1; // -1 means unassigned!

			// REMOVE THIS LATER, AND USE "COORD" SIMPLE ARRAYS INSTEAD!!!
			pdbIter *it = new pdbIter(pdb); // molecule iterator

			// Generating a given number of spiral points uniformly distributed over the surface of a sphere.
			float *sp = NULL; // Sphere points (float triplets), NULL enables memory allocation
			spherePoints(1.0, nsphere, &sp);

			// MAIN ARROW-VISUALIZATION LOOP
			//
			int ncluster; // Number of atoms in current cluster
			int nsamples; // Number of valid samples for current cluster
			int *cluster; // Current cluster
			float *samples = NULL; // valid sample points for visualization
			double *affine = NULL;
			int r; // requested number of clusters
			int c; // cluster index


			//			double **affine_c;
			//			// Pablo max number of clusters 50
			//			affine_c = (double **) allocate(sizeof(double *)*50,"");
			//			for(int ii=0; ii<50; ii++)
			//				affine_c[ii] = (double *) allocate(sizeof(double)*16,"");
			//			double error = 0.0;
			//			double dis_error =0.0;
			//			double temp_e;
			//			double error_C=0;
			//			double error_CF=0;
			//			double Si = 0;
			//			double Si_C = 0;
			//			   float a; //dissimilarity to assigned cluster
			//			    float b; //dissimilarity to nearest cluster not assigned
			//			float s; //silhouette
			//
			//			double bi;
			//
			//			for(r = 0; r < requested; r++) // Screen requested number of clusters
			//			{
			//				Si_C = 0;
			//				error_C=0;
			//				error_CF=0;
			//				// PRECOMPUTE ALL affine matrices
			//				for(int c = 0; c <= r; c++) // Screen clusters
			//				{
			//					// Select cluster
			//					cluster = rclusters[r][c];
			//					ncluster = rnatoms[r][c];
			//					// COMPUTE AFFINE MODEL
			//					// This function is intended to compute the Affine transformation: v = M·x + t
			//					// (it uses homogeneous coordinates, i.e. a 4 elements vector to define a 3D position)
			//					getExpMat(coords, num_atoms, evec+num_comps1*NM_considered, cluster, ncluster, &affine);
			//					for(int ind=0; ind<16; ind++)
			//						affine_c[c][ind]=affine[ind];
			//				}
			//				// Screen clusters
			//				for(int c = 0; c <= r; c++)
			//				{
			//					cluster = rclusters[r][c];
			//					ncluster = rnatoms[r][c];
			//
			//					// be the average distance between i {\displaystyle i} i and all other data points in the same cluster,
			////					double ai =0.0;
			////					for(int i = 0; i < ncluster; i++)
			////					{
			////						it->pos_atom = cluster[i]; // Get atom
			////						(it->get_atom())->getPosition(pos); // get atom position
			////						ai += getError(pos, evecx+it->pos_atom*3, affine_c[c]); // Get atom-wise error (without any normalization)
			////					} // end cluster number
			////					ai /= (ncluster-1.0);
			//
			//					// bi is the smallest average distance of i to all points in any other cluster, of which i is not a member.
			//					bi =100000000000.0;
			//					for(int c2 = 0; c2 <= r; c2++) // Screen other clusters
			//					{
			//						if (c2!=c)	 {
			//							dis_error =0.0;
			//							for(int i = 0; i < ncluster; i++)
			//							{
			//								it->pos_atom = cluster[i];          // Get atom
			//								(it->get_atom())->getPosition(pos); // get atom position
			//								dis_error += getError(pos, evecx+it->pos_atom*3, affine_c[c2]);
			//							}
			//							dis_error /= (ncluster-1.0);;
			//							if (dis_error < bi ) bi = dis_error;
			//						//	{printf ("%5d %5d  %lf bi %lf  ai %lf\n",  c, c2, dis_error, bi, ai  ); }
			//						}
			//					}
			//					// COMPUTE AFFINE MODEL ERROR & DISSIMILATY
			//					error = 0.0; Si = 0;
			//					for(int i = 0; i < ncluster; i++)
			//					{
			//						it->pos_atom = cluster[i]; // Get atom
			//						(it->get_atom())->getPosition(pos); // get atom position
			//
			//						// This calculates the Error (formula XX in Bryden's paper).
			//						// The Error is the distance between the eigenvector component and the affine model evaluation.
			//						temp_e = getError(pos, evecx+it->pos_atom*3, affine_c[c]); // Get atom-wise error (without any normalization)
			//						error += temp_e;
			//						// dissimilariy
			//						if (ai < bi)
			//							temp_e = 1.0 - ai/bi;
			//						else if (ai > bi)
			//							temp_e = bi/ai -1.0;
			//						else temp_e=0.0;
			//						Si += temp_e;
			//						 if ((i<10)&&(c!=0)) {printf ("%5d %5d %2d %lf %lf %lf  %lf\n",  i, ncluster, c,  ai, bi, temp_e, Si ); getchar();}
			//					} // end cluster number
			//					Si_C += Si;
			//
			////				if (r==c)	printf("%2d %2d %5d  %lf (%lf)  %lf (%lf)   %lf -> %lf  MIn %lf %lf   -> %lf %d %d\n",
			////											r+1,c+1, ncluster, error,factor*error,error/ncluster,factor*error/ncluster,
			////												dis_error, ai, bi, Si/(ncluster-1.0), Si_C/(num_atoms-1.0), ncluster,  num_atoms);
			//				// getchar();
			//					error_C += error;
			//					error_CF += factor*error;
			//
			//
			//				}  // screen clusters
			//				printf("%2d %lf %lf %lf %lf   -> %lf\n",r+1, error_C,error_CF,error_C/num_atoms,factor*error_CF/num_atoms, Si_C/(num_atoms-1.0));
			//
			//
			//			}



			for(int r=N_start; r<N_end; r++)
				//			for(r = 0; r < requested; r++) // Screen requested number of clusters
			{
				// Setting which cluster any atom belongs to ("atom_owner")
				for(c = 0; c <= r; c++) // Screen clusters
					for(int i = 0; i<rnatoms[r][c]; i++) // screen cluster's atoms
						atom_owner[ rclusters[r][c][i] ] = c; // Set atom owner (cluster index an atom belongs to)


				for(int c = 0; c <= r; c++) // Screen clusters
				{
					//					printf("imodview> Processing cluster %d (for %d requested clusters)\n",c+1,r+1);

					// Select cluster
					cluster = rclusters[r][c];
					ncluster = rnatoms[r][c];

					// COMPUTE AFFINE MODEL
					// This function is intended to compute the Affine transformation: v = M·x + t
					// (it uses homogeneous coordinates, i.e. a 4 elements vector to define a 3D position)
					getExpMat(coords, num_atoms, evec+num_comps1*NM_considered, cluster, ncluster, &affine);


					// Shows a standard format matrix (standard output)
					if(debug)
						show_matrix_standard(affine, 4, "Affine model matrix for selected cluster");

					// Compute the rotational and non-rotational "kinetic energies"

					// NOTE: Think about determining the center of the motion using the affinity matrix and considering as center point
					//       that having the smallest displacement. Write one equation per component (x,y,z) differentiate and equal to zero...
					//       // Check this function:
					//       multMatrix_4x4_x_3f(affine, pos, av); // Assumes bottom row of matrix is [0,0,0,1]

					//			  // Load "affine" into "A"
					//			  double **A,**Ad;
					//			  A = (double **) allocate(sizeof(double *)*3,"");
					//			  Ad = (double **) allocate(sizeof(double *)*3,"");
					//			  for(int i=0; i<3; i++)
					//			  {
					//				  A[i] = (double *) allocate(sizeof(double)*3,"");
					//				  Ad[i] = (double *) allocate(sizeof(double)*3,"");
					//				  for(int j=0; j<3; j++)
					//					  if(i==j)
					//						  A[i][j] = affine[j*4+i]-1;
					//					  else
					//						  A[i][j] = affine[j*4+i];
					//			  }
					//
					//			  float coma[3];
					//			  for(int k=0; k<3; k++)
					//			  {
					//				  for(int i=0; i<3; i++) // rows
					//					  for(int j=0; j<3; j++) // cols
					//					  {
					//						  if(j != k)
					//							  Ad[i][j] = affine[j*4+i];
					//						  else
					//							  Ad[i][j] = affine[3*4+i];
					//					  }
					//				  coma[k] = Determinant(Ad,3)/Determinant(A,3);
					//			  }
					//			  printf("coma= %f %f %f\n",coma[0],coma[1],coma[2]);

					// COMPUTE THE RIGID-BODY TRANSFORMATION MATRIX AND THE "OTHER" MOTION MATRIX
					double affineRot[16],affineOther[16];
					for(int i=0; i<16; i++)
					{
						affineRot[i] = 0.0;
						affineOther[i] = 0.0;
					}
					//	  expMatRot(1,2) = (expMat(1,2)-expMat(2,1))/2.0;
					//	  expMatRot(0,2) = (expMat(0,2)-expMat(2,0))/2.0;
					//	  expMatRot(0,1) = (expMat(0,1)-expMat(1,0))/2.0;
					affineRot[2*4+1] = (affine[2*4+1] - affine[1*4+2])/2.0;
					affineRot[2*4+0] = (affine[2*4+0] - affine[0*4+2])/2.0;
					affineRot[1*4+0] = (affine[1*4+0] - affine[0*4+1])/2.0;
					//	  expMatRot(2,1) = -(expMat(1,2)-expMat(2,1))/2.0;
					//	  expMatRot(2,0) = -(expMat(0,2)-expMat(2,0))/2.0;
					//	  expMatRot(1,0) = -(expMat(0,1)-expMat(1,0))/2.0;
					affineRot[1*4+2] = -affineRot[2*4+1];
					affineRot[0*4+2] = -affineRot[2*4+0];
					affineRot[0*4+1] = -affineRot[1*4+0];
					//	  expMatRot(0,3) = expMat(0,3);
					//	  expMatRot(1,3) = expMat(1,3);
					//	  expMatRot(2,3) = expMat(2,3);
					affineRot[3*4+0] = affine[3*4+0];
					affineRot[3*4+1] = affine[3*4+1];
					affineRot[3*4+2] = affine[3*4+2];

					// The difference should be the non-rigid transformation...
					for(int i=0; i<16; i++)
						affineOther[i] = affine[i] - affineRot[i];

					//draw the axis of rotation direction
					//		MathVec3D<float> axisVec;
					//		axisVec[0] = -(expMat(1,2));
					//		axisVec[1] = expMat(0,2);
					//		axisVec[2] = -expMat(0,1);

					// GET AXIS OF ROTATION FROM AFFINE MODEL
					float axisVec[3],normAxis[3];
					//	  affine = affineRot; // MEGAÑAPA!!!
					//	  affine = affineOther; // MEGAÑAPA!!!

					//			  axisVec[0] = -affine[2*4+1]; // (1,2)
					//			  axisVec[1] = affine[2*4+0]; // (0,2)
					//			  axisVec[2] = -affine[1*4+0]; // (0,1)
					axisVec[0] = -affineRot[2*4+1]; // (1,2)
					axisVec[1] = affineRot[2*4+0]; // (0,2)
					axisVec[2] = -affineRot[1*4+0]; // (0,1)

					normalize(axisVec,normAxis);
					if(debug)
						fprintf(stdout,"axisVec= %e %e %e   normAxis= %f %f %f\n",axisVec[0],axisVec[1],axisVec[2],normAxis[0],normAxis[1],normAxis[2]);

					// DRAW AFFINE MODEL ROTATION AXIS
					if(plot_axis)
					{
						FILE *f_axis;
						//				  sprintf(text,"axis_r%02d_c%02d.vmd",r+1,c+1);
						sprintf(text,"%saxis_n%d_r%d_c%d.vmd",name,NM_considered+1,r+1,c+1);
						if( (f_axis = fopen(text, "w")) == NULL)
						{
							printf("imodview> Affine-model-axis-VMD-file writing error!! \n");
							exit(1);
						}
						fprintf(f_axis,"molecule new\ndraw color red\ndraw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius 0.200 resolution 20 filled 1\ndisplay resetview"
								,-normAxis[0]*asize+com[0],-normAxis[1]*asize+com[1],-normAxis[2]*asize+com[2],normAxis[0]*asize+com[0],normAxis[1]*asize+com[1],normAxis[2]*asize+com[2]);
						//				  fprintf(f_axis,"molecule new\ndraw color red\ndraw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius 0.200 resolution 20 filled 1\ndisplay resetview"
						//						  ,-normAxis[0]*asize+comf[0],-normAxis[1]*asize+comf[1],-normAxis[2]*asize+comf[2],normAxis[0]*asize+comf[0],normAxis[1]*asize+comf[1],normAxis[2]*asize+comf[2]);
						//				  fprintf(f_axis,"molecule new\ndraw color red\ndraw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius 0.200 resolution 20 filled 1\ndisplay resetview"
						//						  ,-normAxis[0]*asize+coma[0],-normAxis[1]*asize+coma[1],-normAxis[2]*asize+coma[2],normAxis[0]*asize+coma[0],normAxis[1]*asize+coma[1],normAxis[2]*asize+coma[2]);
						fclose(f_axis);
					}

					// COMPUTING THE AFFINE MODEL NORMALIZATION FACTOR
					int i;
					double accum = 0.0;
					double *eveca;
					max_amp = 0.0;
					int affine_norm_method = 1; // Set affine model normalization method

					if(affine_norm_method == 0) // Do not normalize at all!
					{
						if(debug)
							printf("imodview> NO Normalization!\n");
						factor = 1.0;
					}
					else if(affine_norm_method == 1) // Normalize to set the mode component of maximum amplitude = 1A.
					{
						for(int n = 0; n < ncluster; n++)
						{
							i = cluster[n]; // index for atoms belonging to current cluster
							amp = magnitude(evecx+i*3);
							if ( amp > max_amp ) // searches the maximum amplitude to normalize
							{
								max_amp = amp;
								max_amp_atom = i;
							}
						}
						factor = (1.0/max_amp)*move_delta; // now, factor normalizes the maximum displacement to 1A.
						if(debug)
							printf("imodview> Cluster's maximum amplitude %lf\t for atom %d. factor= %f  delta= %.2f A\n",max_amp,max_amp_atom,factor,move_delta);
					}
					else if(affine_norm_method == 2) // Normalize to the cluster's mode norm
					{
						for(int i = 0; i < ncluster; i++)
						{
							eveca = (evecx+i*3);
							accum += eveca[0]*eveca[0] + eveca[1]*eveca[1] + eveca[2]*eveca[2];
						}
						max_amp = sqrt(accum); // N-dimensional vector norm
						factor = (1.0/max_amp)*move_delta; // now, factor normalizes to make the norm = 1.
						if(debug)
							printf("imodview> Norm of cluster's mode components: %lf. factor= %f  delta= %.2f A\n",max_amp,factor,move_delta);
					}


					//			  speed = 0.0; // speed summation for current cluster
					//			  float speedRot = 0.0;
					//			  float speedOther = 0.0;

					// DRAW ARROWS FROM AFFINE MODEL
					if(plot_affine)
					{
						float av[3]; // vector computed from affine model matrix
						float p1[3],p2[3];
						FILE *f_cluster;
						if(format == 1) // VMD
							sprintf(text,"%saffine_n%d_r%d_c%d.vmd",name,NM_considered+1,r+1,c+1);
						else // Jmol
							sprintf(text,"%saffine_n%d_r%d_c%d.txt",name,NM_considered+1,r+1,c+1);
						if( (f_cluster = fopen(text, "w")) == NULL)
						{
							printf("imodview> Affine-model-Cluster-VMD-file writing error!! \n");
							exit(1);
						}
						if(format == 1) // VMD
							fprintf(f_cluster,"molecule new\ndraw color %s\n",colortext);

						for(int i = 0; i < ncluster; i++)
						{
							it->pos_atom = cluster[i]; // Get atom
							(it->get_atom())->getPosition(pos); // get atom position

							// Apply affine transformation to given coordinates (pos)
							multMatrix_4x4_x_3f(affine, pos, av); // Assumes bottom row of matrix is [0,0,0,1]

							// get final point
							p2[0] = pos[0] + av[0]*factor;
							p2[1] = pos[1] + av[1]*factor;
							p2[2] = pos[2] + av[2]*factor;

							// Draw computed (and scaled) velocities...

							//					  speed += pow(magnitude(av,3),2); // speeds are proportional to eigenvector magnitude
							//					  // Apply affineRot transformation to given coordinates (pos)
							//					  multMatrix_4x4_x_3f(affineRot, pos, av); // Assumes bottom row of matrix is [0,0,0,1]
							//					  speedRot += pow(magnitude(av,3),2);
							//					  // Apply affineOther transformation to given coordinates (pos)
							//					  multMatrix_4x4_x_3f(affineOther, pos, av); // Assumes bottom row of matrix is [0,0,0,1]
							//					  speedOther += pow(magnitude(av,3),2);

							if(format == 1) // VMD
							{
								// Drawing cylinder
								p1[0] = pos[0] + av[0]*factor*0.7;
								p1[1] = pos[1] + av[1]*factor*0.7;
								p1[2] = pos[2] + av[2]*factor*0.7;
								fprintf(f_cluster,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
										pos[0],pos[1],pos[2],p1[0],p1[1],p1[2],radius);
								// Drawing cone (it's an arrow!)
								fprintf(f_cluster,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
										p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],radius*3);
							}
							else // JMol
							{
								fprintf(f_cluster,"draw ar%d arrow FIXED {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} color %s width %.4f\n",i+1,
										pos[0],pos[1],pos[2],p2[0],p2[1],p2[2],colortext,radius*3);
							}
						}
						if(format == 1) // VMD
							fprintf(f_cluster,"display resetview\n");
						fclose(f_cluster);

						//				  printf("r= %d  c= %d  speed= %f  speedRot= %f (%f)  speedOther= %f (%f)\n",r,c,speed,speedRot,speedRot/speed,speedOther,speedOther/speed);
					}

					// Computing cluster center of mass
					float pos[3];
					com[0] = com[1] = com[2] = 0.0;
					for(int i = 0; i<ncluster; i++)
					{
						it->pos_atom = cluster[i];
						(it->get_atom())->getPosition(pos);
						com[0] += pos[0];
						com[1] += pos[1];
						com[2] += pos[2];
					}
					com[0] /= ncluster;
					com[1] /= ncluster;
					com[2] /= ncluster;

					// SELECTING VALID SAMPLING POINTS (those not occluded by any neighboring atom)
					//  First, the homogeneously-distributed spherical unit vectors are centered at current cluster CoM.
					//  Next, the positions of the atoms are projected into every unit vector. Only those atoms close to the vectors are considered.
					//  Then, the maximum projection found for each vector determines the potential valid point.
					//  Finally, any potential valid point non-occluded by any atom is separated some distance outwards and stored as valid sampling point.
					if(plot_sphere)
					{
						//				  sprintf(text,"sphere_r%02d_c%02d.vmd",r+1,c+1);
						sprintf(text,"%ssphere_n%d_r%d_c%d.vmd",name,NM_considered+1,r+1,c+1);
						getValid(cluster, ncluster, c, com, &samples, &nsamples, sp, nsphere, coords, coords2, atom_owner, num_atoms, surf_dist, close2vector_thr, close2clusters_thr, text);
					}
					else
						getValid(cluster, ncluster, c, com, &samples, &nsamples, sp, nsphere, coords, coords2, atom_owner, num_atoms, surf_dist, close2vector_thr, close2clusters_thr);

					// THIS PLOTS THE VALID POINTS (DO NOT DELETE)
					if(plot_valid)
					{
						FILE *f_points;
						sprintf(text,"%svalid_n%d_r%d_c%d.vmd",name,NM_considered+1,r+1,c+1);
						if( (f_points = fopen(text, "w")) == NULL)
						{
							printf("imodview> Valid points file writing error!! \n");
							exit(1);
						}
						fprintf(f_points,"molecule new\ndraw color %s\n",colortext);
						for(int i=0;i<nsamples;i++)
							fprintf(f_points,"draw sphere \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",samples[i*3],samples[i*3+1],samples[i*3+2],0.6);
						fprintf(f_points,"display resetview\n");
						fclose(f_points);
					}

					// COMPUTING THE LENGTH-PER-UNITs FOR EVERY VALID SAMPLING POINT
					float *lengthPerUnit; // Array with the length-per-unit of every valid point
					lengthPerUnit = (float *) allocate(sizeof(float) * nsphere, ""); // allocating the maximum memory possible...
					// Compute the length-per-unit given one affine model and one samples set.
					getLengthPerUnit(lengthPerUnit, samples, nsamples, affine, integration_step);

					// COMPUTING VELOCITY-WEIGHTED CENTER OF MASS FOR CURRENT CLUSTER
					// This improves the scoring used in best-valid-point determination
					float tot_weight = 0.0;
					float weight;
					float comf[3]; // Cluster (weighted) center of mass (float)
					comf[0] = 0.0;
					comf[1] = 0.0;
					comf[2] = 0.0;
					for(int i = 0; i < ncluster; i++)
					{
						// Get weight
						weight = pow(magnitude(evec+num_comps1*NM_considered+cluster[i]*3),2);
						tot_weight += weight;
						// Compute CoM
						it->pos_atom = cluster[i];
						(it->get_atom())->getPosition(pos);
						comf[0] += weight * pos[0];
						comf[1] += weight * pos[1];
						comf[2] += weight * pos[2];
					}
					comf[0] /= tot_weight;
					comf[1] /= tot_weight;
					comf[2] /= tot_weight;
					if(debug)
						printf("weighted_com= %f %f %f    standard_com= %f %f %f\n",comf[0],comf[1],comf[2],com[0],com[1],com[2]);

					// LOOK FOR THE BEST VALID POINT
					float jmax_best,jmin_best; // Max/Min path amplitudes for the best trajectory
					float pathLength_best; // Best path length
					int i_best; // Best sample index
					if(plot_paths)
					{
						//				  sprintf(text,"paths_r%02d_c%02d.vmd",r+1,c+1);
						sprintf(text,"%spaths_n%d_r%d_c%d.vmd",name,NM_considered+1,r+1,c+1);
						// Get the best valid point (i.e. the one that leads to the best path)
						getBest(affine, samples, nsamples, cluster, ncluster, c, comf, lengthPerUnit, coords, far2path_thr, close2path_thr, delta_sampling, atom_owner, num_atoms, integration_step, &i_best, &jmax_best, &jmin_best, &pathLength_best, text);
					}
					else
						// Get the best valid point (i.e. the one that leads to the best path)
						getBest(affine, samples, nsamples, cluster, ncluster, c, comf, lengthPerUnit, coords, far2path_thr, close2path_thr, delta_sampling, atom_owner, num_atoms, integration_step, &i_best, &jmax_best, &jmin_best, &pathLength_best);

					float *sample = samples + i_best*3; // best sample pointer

					//
					// PLOT THE BEST PATH AND THE CORRESPONDING ARROW/S
					//

					// Set arrow length for current cluster
					float arrow_length = arrow_to_path * pathLength_best * rspeeds[r][c] / rspeeds_max[r]; // Arrow length in Angstroms

					// Set aspect ratio (for constant arrow-tip size)
					if(arrow_length > tip_length)
						arrow_aspect = 1-(tip_length/arrow_length);
					else // if tip is bigger than the arrow itself...
						arrow_aspect = min_aspect;

					// Draws arrow/s and path from affine model
					//			  sprintf(text,"best_r%02d_c%02d",r+1,c+1);
					sprintf(text,"%sarrow_n%d_r%d_c%d",name,NM_considered+1,r+1,c+1); // file-extension will be set by "drawArrow()"
					//			printf("imodview> Drawing arrow and trajectory for:  r= %d  c= %d\n",r+1,c+1);
					//			  printf("OUTSIDE: pathLength_best= %f  arrow_length= %f  delta_path= %f\n",pathLength_best,arrow_length,delta_path);

					drawArrow(format, samples + i_best*3, affine, jmin_best, jmax_best, pathLength_best, lengthPerUnit[i_best], comf,
							delta_path, integration_step,	text, arrow_number, arrow_length, ribbon_width, arrow_width, tip_width,
							tip_length, arrow_aspect, arrow_to_path, min_aspect, arrows_apart, path_middle_lines);

					free(lengthPerUnit);
				}


			}
			break;
		}
		default:
			printf("imodview> Please introduce a valid operation mode!\n");
			exit(0);
			break;
	}

	fclose(Fout);
	return 0;
}
// Main ENDS...

void parseOptions(int argc, char** argv)
{
	std::string temp;
	CmdLine cmd("imodview","Cartesian Normal Modes and Springs Visualization tool", version );

	try
	{
		//
		// Define required arguments no labeled
		//
		UnlabeledValueArg<std::string> pdb("pdb","PDB input file. Warning: It must match exactly that used in NMA.","default","pdb");
		cmd.add( pdb );
		UnlabeledValueArg<std::string> nms("nms","Input eigenvectors (ptraj) or force constants (Kfile) file. Warning: Only Cartesian modes allowed.","default","ptraj/Kfile");
		cmd.add( nms );
		UnlabeledValueArg<std::string> out("out","Ouput VMD file.","default","vmd");
		cmd.add( out );

		SwitchArg Server("","server", "Enable \"server mode\" to maximize the automatic selection of parameters. (default=disabled)", false);
		cmd.add( Server );

		// Arrow drawing related parameters
		SwitchArg PlotSphere("","save_sphere", "Plot the points homogeneously distributed over a sphere centered on each cluster. (default=disabled)", false);
		cmd.add( PlotSphere );
		SwitchArg PlotValid("","save_valid", "Plot the valid points, i.e. those not colliding with any atom. (default=disabled)", false);
		cmd.add( PlotValid );
		SwitchArg PlotPaths("","save_paths", "Plot the full paths obtained from every valid point. (default=disabled)", false);
		cmd.add( PlotPaths );
		SwitchArg PlotAxis("","save_axis", "Plot the cluster axis or rotation assuming rigid body motion. (default=disabled)", false);
		cmd.add( PlotAxis );
		SwitchArg PlotAffine("","save_affine", "Plot the arrows obtained from the calculated affine model (in arbitrary units). (default=disabled)", false);
		cmd.add( PlotAffine );
		ValueArg<int> Format("","format", "Arrows output format: 1=VMD, 2=JMol(PMESH). Only for arrows. (Default=JMol)",false,2,"int");
		cmd.add( Format );

		// Clustering related parameters
		SwitchArg PDBClusters("","save_pdb_clusters", "Save pdbs with cluster information stored in the B-factors. (default=disabled)", false);
		cmd.add( PDBClusters );
		ValueArg<float> ExtendCutoff("","extend_cutoff", "Distance (in Angstroms) to extend small clusters. (Default=10A)",false,10,"float");
		cmd.add( ExtendCutoff);
		ValueArg<int> ExtendSize("","extend_size", "Number of atoms threshold to consider cluster extension. Only clusters with less or equal number atoms will be extended. (Default=5)",false,5,"int");
		cmd.add( ExtendSize );
		ValueArg<float> Adjacency("","adj_cutoff", "Adjacency cutoff (in Angstroms). (Default=6A)",false,6,"float");
		cmd.add( Adjacency);
		ValueArg<std::string> FixFile("f","fixFile", "Input ASCII file defining the ICs that were fixed during NMA (default=disabled). "
				"If modes were computed removing arbitrary ICs, the user must introduce here the file generated by iMode's --save_fixfile option.",false,"fixstring","string");
		cmd.add( FixFile );
		ValueArg<int> NClusters("","nclusters", "Number of requested clusters. (Default=15)",false,15,"int");
		cmd.add( NClusters );

		//
		// Define labeled arguments
		//
		ValueArg<std::string> Name("o","name", "Output files basename (default=disabled).",false,"","string");
		cmd.add( Name );
		ValueArg<float> KThr2("","kthr2", "Only those springs with force constants below this threshold will be shown (default=disabled).",false,999999.0,"float");
		cmd.add( KThr2 );
		ValueArg<float> KThr("","kthr", "Only those springs with force constants above this threshold will be shown (default=disabled).",false,0.0,"float");
		cmd.add( KThr );
		ValueArg<float> percent("","pthr", "Minimum percentual amplitude (from maximum) to show arrows (default=0, all arrows).",false,0.0,"float");
		cmd.add( percent );
		ValueArg<int> NResidues("","nresidues", "Sets the number of residues for averaging (only used by \"--level 1\" option (default=1).",false,1,"int");
		cmd.add( NResidues );
		ValueArg<int> Level("","level", "Sets the averaging level to compute arrows, 0=atoms, 1=residues, 2=segments, 3=chains (default=0).",false,0,"int");
		cmd.add( Level );
		ValueArg<float> Radius("","thick", "Arrow/spring thickness factor (default=0.05).",false,0.05,"float");
		cmd.add( Radius );
		ValueArg<float> max("","max", "Maximum arrow length [A] (default=10).",false,10.0,"float");
		cmd.add( max );
		ValueArg<int> Operation("","op", "Sets the operation method, 1=Arrows, 2=Springs, 3=Bryden's (default=1).",false,1,"int");
		cmd.add( Operation );
		ValueArg<std::string> color("","color", "Set color. All VMD colors available (default=black).",false,"black","string");
		cmd.add( color );
		ValueArg<int> numNM("n","nev", "Normal Mode index (1,2,...,size) (default=1).",false,1,"int");
		cmd.add( numNM );

		// Parse the command line2.
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		strcpy(file_pdb,((temp=pdb.getValue()).c_str())); // Gets pdb name
		strcpy(file_nma,((temp=nms.getValue()).c_str())); // Gets input file name
		strcpy(file_out,((temp=out.getValue()).c_str())); // Gets output file name
		// strcpy(name,((temp=Name.getValue()).c_str())); // Gets Basename
		if(Name.isSet())
			sprintf(name,"%s_",(temp=Name.getValue()).c_str()); // Gets Basename
		else
			name[0] = '\0';

		strcpy(colortext,((temp=color.getValue()).c_str())); // Gets the color name.
		move_delta = max.getValue(); // default displacement
		percent_value = percent.getValue(); // default minimum percentual threshold to visualice.
		operation = Operation.getValue();
		nresidues = NResidues.getValue();
		arrow_level = Level.getValue();
		kthr = KThr.getValue();
		kthr2 = KThr2.getValue();
		radius = Radius.getValue();
		NM_considered = --numNM.getValue();
		if (NClusters .isSet())
			requested = NClusters.getValue(); // Number of requested clusters
		else {
			NClusters_auto = true;
			requested = 15;
		}


		if(FixFile.isSet())
			strcpy(file_fix,((temp=FixFile.getValue()).c_str())); // Gets Fix-file
		adjacency = Adjacency.getValue(); // Adjacency cutoff (in Angstroms)
		extend_size = ExtendSize.getValue(); // Number of atoms threshold to consider cluster extension (clusters with less atoms will be extended). (Default=5)
		extend_cutoff = ExtendCutoff.getValue(); // Distance (in Angstroms) to extend small clusters. (Default=10A)
		pdb_clusters = PDBClusters.getValue(); // Output clusters as the b-factors of pdbs

		format = Format.getValue(); // Arrows output format: 1=VMD, 2=JMol(PMESH). Only for arrows. (Default=JMol)
		plot_sphere = PlotSphere.getValue(); // Plot spheres
		plot_valid = PlotValid.getValue(); // Plot valid points
		plot_paths = PlotPaths.getValue(); // Plot valid paths before choosing the best one
		plot_axis = PlotAxis.getValue(); // Plot best axis
		plot_affine = PlotAffine.getValue(); // Plot affine model arrows
		server = Server.getValue(); // Server mode
	}
	catch ( ArgException& e )
	{
		std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl;
	}
}


// Selecting valid points for a given cluster (those not occluded by any neighboring atom)
// First, the homogeneously-distributed spherical unit vectors are centered at current cluster CoM.
// Next, the positions of the atoms are projected into every unit vector. Only those atoms close to the vectors are considered.
// Then, the maximum projection found for each vector determines the potential valid point.
// Finally, any potential valid point non-occluded by any atom is separated some distance outwards and stored as valid sampling point.
//
//  cluster --> Current cluster
//  ncluster --> Number of cluster atoms
//  c --> Cluster index
//  com --> Cluster center of mass
//  p_samples --> Pointer to the valid points array (OUTPUT)
//  p_nsamples --> Pointer to the number of valid points (OUTPUT)
//  sp --> Spherical distribution of points
//  nsphere --> Number of spherical distribution of points
//  coords --> Full molecule coordinates
//  coords2 --> Full molecule coordinates to be centered at cluster center
//  atom_owner --> Atomic cluster ownership array
//  num_atoms --> Macromolecule number of atoms
//  surf_dist = 5.0; // Distance from the "surface"
//  close2vector_thr = 10.0; // Threshold to determine which points are close to a given vector.
//  close2clusters_thr = 5.0; // threshold to determine which points are close to any different cluster atom.
//  f_name --> Output filename for the centered spherical distribution of points (OPTIONAL)
void getValid(int *cluster, int ncluster, int c, double *com, float **p_samples, int *p_nsamples, float *sp, int nsphere, float *coords, float *coords2, int *atom_owner, int num_atoms, float surf_dist, float close2vector_thr, float close2clusters_thr, char *f_name)
{
	float dummy[3];

	float *samples;
	if(*p_samples == NULL)
		*p_samples = (float *) allocate( sizeof(float)*nsphere*3, ""); // Allocating the maximum memory possible
	samples = *p_samples; // already allocated

	// Pre-centering atomic coordinates in the corresponding cluster (to speed up radial-vectors projection)
	for(int i=0; i<ncluster; i++)
	{
		coords2[3*cluster[i] ] -= com[0];
		coords2[3*cluster[i] + 1] -= com[1];
		coords2[3*cluster[i] + 2] -= com[2];
	}

	// THIS PLOTS THE SPHERE at cluster center (DO NOT DELETE)
	if(f_name != NULL) // if output file requested
	{
		FILE *f_sphere;
		if( (f_sphere = fopen(f_name, "w")) == NULL)
		{
			printf("imodview> Affine-model-axis-VMD-file writing error!! \n");
			exit(1);
		}
		fprintf(f_sphere,"molecule new\ndraw color %s\n",colortext);
		for(int i=0;i<nsphere;i++)
			fprintf(f_sphere,"draw sphere \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",sp[i*3]*rsphere+com[0],sp[i*3+1]*rsphere+com[1],sp[i*3+2]*rsphere+com[2],0.5);
		fprintf(f_sphere,"display resetview\n");
		fclose(f_sphere);
	}

	// Project cluster's atoms positions over the spherical radial vectors to obtain valid sample points for visualization
	float close2clusters_thr2 = pow(close2clusters_thr,2); // Squared threshold to determine which points are close to any different cluster atom.
	float proj_max; // maximum projection over the unitary radial vectors
	float proj; // current projection
	float dist; // current minimum distance
	bool good;
	int nsamples = 0; // valid samples number
	for(int i=0; i<nsphere; i++) // screen unitary radial vectors
	{
		proj_max = 0.0; // reset maximum projection
		for(int j=0; j<ncluster; j++) // screen cluster atoms
		{
			// Get the minimum distance between atomic coordinates and the unitary radial vectors:
			// Distance between point P and a line passing by (0,0,0) with direction B:    d = | P x B | / |B|
			// (i.e. projecting every cluster atom over every radial vector)
			crossProduct(coords2 + 3*cluster[j], sp + 3*i, dummy);
			dist = magnitude(dummy);

			if(dist <= close2vector_thr) // if atoms are near current vector
			{
				// Get the projection over current vector
				proj = dotProduct(sp + i*3, coords2 + 3*cluster[j]); // projection between unitary radial vector and centered atomic positions

				// and search the maximum
				if(proj > 0 && proj > proj_max)
					proj_max = proj;
			}
		}

		good = true;
		if(proj_max > 0.0)
		{
			proj_max += surf_dist; // to separate valid sample points from cluster surface

			// Get the maximum projection point over current radial vector
			dummy[0] = sp[i*3]     * proj_max + com[0];
			dummy[1] = sp[i*3 + 1] * proj_max + com[1];
			dummy[2] = sp[i*3 + 2] * proj_max + com[2];

			// Check whether it is too close to any atom from a different cluster
			for(int k=0; k<num_atoms && good; k++) // screen all atoms (somewhat inefficient...)
				if(atom_owner[k] != c) // if atom does not belong to current cluster
				{
					// Check distance
					if(sqrDist(dummy,coords+3*k) < close2clusters_thr2) // if too close from other clusters...
						good = false;
				}

			if(good) // Good point found (if current projection point is not close to other clusters...)
			{
				// then add projection over the radial vector to the valid sample points array
				samples[nsamples*3] 	  = dummy[0];
				samples[nsamples*3 + 1] = dummy[1];
				samples[nsamples*3 + 2] = dummy[2];
				nsamples++;
			}
		}
	}

	// Moving back the atomic coordinates into the original position (for next cluster calculations)
	for(int i=0; i<ncluster; i++)
	{
		coords2[3*cluster[i] ] += com[0];
		coords2[3*cluster[i] + 1] += com[1];
		coords2[3*cluster[i] + 2] += com[2];
	}

	*p_nsamples = nsamples; // outputs number of valid samples

	if(nsamples == 0)
		fprintf(stderr,"Warning: No valid sample points found! Check thresholds!\n");

}

// Compute the length-per-unit given one affine model and one samples set.
//  lengthPerUnit --> Array with the lengths per unit (it should have been already allocated)
//  samples --> Valid sample points set
//  nsamples --> Number of valid sample points
//  affine --> Affine model
//  integration_step --> Affine model integration step
void getLengthPerUnit(float *lengthPerUnit, float *samples, int nsamples, double *affine, float integration_step)
{
	bool debug = false;
	float *sample;
	float *buffer;
	double affine_new[16];
	float *pos1,*pos2;
	pos1 = (float *) allocate(sizeof(float) * 3, "");
	pos2 = (float *) allocate(sizeof(float) * 3, "");

	//	  // THIS PLOTS THE ARC-LENGHT ESTIMATION POINTS (DO NOT DELETE)
	//	  FILE *f_arc;
	//	  if( (f_arc = fopen("arc.vmd", "w")) == NULL)
	//	  {
	//		  printf("imodview> Arc points file writing error!! \n");
	//		  exit(1);
	//	  }
	//	  fprintf(f_arc,"molecule new\ndraw color %s\n",colortext);

	// "ArcLengthPerUnit" COMPUTATION
	for(int i=0; i<nsamples; i++) // Screening all valid sample points
	{
		lengthPerUnit[i] = 0.0;
		sample = samples + 3*i; // current sample point
		for(int j=0; j<3; j++)
			pos1[j] = sample[j];
		for(float j=0.2; j<=1; j+=0.2) // 5 interpolation points are enough...
		{
			// Integrate affine model into affine_new
			getExp(affine,affine_new,j*integration_step,10);
			//	show_matrix_standard(affine_new, 4, "Affine_new");

			// Apply affine transformation to sample point
			multMatrix_4x4_x_3f(affine_new, sample, pos2); // Assumes bottom row of matrix is [0,0,0,1]

			lengthPerUnit[i] += Dist(pos1,pos2);

			//			  // THIS PLOTS THE ARC-LENGHT ESTIMATION POINTS (DO NOT DELETE)
			//			  if(dump_draw)
			//				  fprintf(f_arc,"draw sphere \"%8.5f %8.5f %8.5f\" radius %.4f resolution 10\n",pos2[0],pos2[1],pos2[2],0.2);

			// Swapping positions
			buffer = pos1;
			pos1 = pos2;
			pos2 = buffer;
		}
		if(debug)
			printf("Sample point %4d  LengthPerUnit= %f\n",i,lengthPerUnit[i]);
	}
	//	  // THIS PLOTS THE ARC-LENGHT ESTIMATION POINTS (DO NOT DELETE)
	//	  fprintf(f_arc,"display resetview\n");
	//	  fclose(f_arc);
}

// Get the best valid point (the one that leads to the best path)
void getBest(double *affine, float *samples, int nsamples, int *cluster, int ncluster, int c, float *comf, float *lengthPerUnit, float *coords, float far2path_thr, float close2path_thr, float delta_sampling, int *atom_owner, int num_atoms, float integration_step, int *p_i_best, float *p_jmax_best, float *p_jmin_best, float *p_pathLength_best, char *paths_name)
{
	bool debug = false;
	bool collision = false; // if true, it stops path integration
	float *pathpoints; // Current path points
	float pathLength; // Current path length
	float pathCenter[3]; // Current path center
	int maxpathpoints; // Maximum number of points per path
	int npathpoints; // Current number of points in path trajectory
	float jmax,jmin; // Maximum and minimum path amplitude (in affine matrix units and without multiplying by integration_step)
	float jmax_best,jmin_best; // Max/Min path amplitudes for the best trajectory
	float score,score_best; // Scores...
	float pathLength_best; // Best path length
	int i_best = -1; // Best sample index
	float dist_min;
	float *buffer;
	double affine_new[16];
	float *sample; // sample point coordinates pointer
	float *pos1,*pos2;

	float far2path_thr2 = pow(far2path_thr,2);
	float close2path_thr2 = pow(close2path_thr,2);
	int pathgap = (int) ceil(close2path_thr/delta_sampling); // Number of path points excluded from collision checking

	pos1 = (float *) allocate(sizeof(float) * 3, "");
	pos2 = (float *) allocate(sizeof(float) * 3, "");

	// DETERMINE THE MAXIMUM PATH LENGTH (from maximum cluster inter-atomic distance)
	float dist;
	float maxLength = -1.0; // Maximum path length integration amplitude (in Angstroms)
	for(int i=0; i<ncluster; i++) // screen cluster atoms
		for(int j=i+1; j<ncluster; j++) // screen cluster atoms
		{
			dist = sqrDist(coords+cluster[i]*3, coords+cluster[j]*3);
			if(dist > maxLength)
				maxLength = dist;
		}
	maxLength = sqrt(maxLength) * maxLengthFactor; // from squared distance to distance...

	// Maximum number of points per path trajectory
	maxpathpoints = 2 * (int) ((float) ceil(1.5*maxLength) / (float) delta_sampling); // "2" due to + and - path senses...
	// Allocating memory for current path points (to check later for self path collision)
	pathpoints = (float *) allocate( sizeof(float) * maxpathpoints * 3, "");
	// Initialization
	for(int i=0; i < maxpathpoints * 3; i++)
		pathpoints[i] = 0.0;

	if(debug)
		printf("Maximum path length: %f  Maximum path points: %d  c= %d\n",maxLength,maxpathpoints,c);

	// THIS PLOTS THE LONGEST PATHS (DO NOT DELETE)
	FILE *f_paths;
	if(paths_name != NULL)
	{
		if( (f_paths = fopen(paths_name, "w")) == NULL)
		{
			printf("imodview> Arc points file writing error!! \n");
			exit(1);
		}
		fprintf(f_paths,"molecule new\ndraw color %s\n",colortext);
	}

	// TRACING VALID-POINTS PATHS
	score_best = 0.0; // reset best score
	for(int i=0; i<nsamples; i++) // Screening all valid sample points
	{
		sample = samples + i*3; // current sample point
		npathpoints = 0; // counts the number of path points
		pathLength = 0.0; // Measure current path length
		pathCenter[0] = 0.0; // Reset path center
		pathCenter[1] = 0.0;
		pathCenter[2] = 0.0;

		// "+" motion sense...
		collision = false;
		pos1[0] = sample[0]; // Set "pos1" to current motion start point
		pos1[1] = sample[1];
		pos1[2] = sample[2];
		for(float j=delta_sampling/lengthPerUnit[i]; j <= (maxLength+0.1)/lengthPerUnit[i] && !collision; j += delta_sampling/lengthPerUnit[i]) // Integrating at 1A length steps
		{
			// Integrate affine model into affine_new
			getExp(affine,affine_new,j*integration_step,10);

			// Apply affine transformation to sample point
			multMatrix_4x4_x_3f(affine_new, sample, pos2); // Assumes bottom row of matrix is [0,0,0,1]

			if(paths_name != NULL)
				fprintf(f_paths,"draw sphere \"%8.5f %8.5f %8.5f\" radius %.4f resolution 10\n",pos2[0],pos2[1],pos2[2],0.2);

			// Check distance between current point and all molecule atoms
			dist_min = 9999999.9; // reset minimal distance between trajectory point and cluster atoms
			for(int k=0; k<num_atoms && !collision; k++)
			{
				dist = sqrDist(pos2,coords+k*3);
				if(dist < close2path_thr2) // path occlusion test...
				{
					collision = true;
					// printf("Path point %d collided with atoms: sqrt(dist)= %f  pos2= %f %f %f  coords[k*3]= %f %f %f\n",npathpoints,sqrt(dist),pos2[0],pos2[1],pos2[2],coords[k*3],coords[k*3+1],coords[k*3+2]);
				}
				if(atom_owner[k] == c && dist < dist_min)
					dist_min = dist; // Update the minimal distance with current cluster atoms
			}

			if(!collision && dist_min > far2path_thr2)// if current trajectory point is too far from the nearest atom of current cluster
			{
				collision = true; // it does not collide, but should not be considered
				// printf("Path point %d too far from any cluster atom: sqrt(dist_mon)= %f  pos2= %f %f %f  min(%d)= %f %f %f\n",npathpoints,sqrt(dist_min),pos2[0],pos2[1],pos2[2],mink,min[0],min[1],min[2]);
			}

			// Check distance between current point and the initial ones
			for(int k=0; k < npathpoints-pathgap && !collision; k++)
			{
				dist = sqrDist(pos2,pathpoints+3*k);
				if(dist < close2path_thr2)
				{
					collision = true;
					// printf("Path point %d collided with any other path point: sqrt(dist)= %f\n",npathpoints,sqrt(dist));
				}
			}

			if(!collision) // Adding a valid path point
			{
				if(npathpoints < maxpathpoints) // This should avoid any possible memory overflow...
				{
					// Current path length
					pathLength += Dist(pos1,pos2);

					jmax = j; // this way it stores the last valid path amplitude
					pathpoints[3*npathpoints] = pos2[0];
					pathpoints[3*npathpoints + 1] = pos2[1];
					pathpoints[3*npathpoints + 2] = pos2[2];
					npathpoints++;
				}
				else
					collision = true; // stop adding points
			}

			// Swapping pos1 by pos2
			buffer = pos1;
			pos1 = pos2;
			pos2 = buffer;
		}

		// "-" motion sense...
		collision = false;
		pos1[0] = sample[0]; // Set "pos1" to current motion start point
		pos1[1] = sample[1];
		pos1[2] = sample[2];
		for(float j= -delta_sampling/lengthPerUnit[i]; j >= -(maxLength+0.1)/lengthPerUnit[i] && !collision; j -= delta_sampling/lengthPerUnit[i]) // Integrating at 1A length steps
		{
			// Integrate affine model into affine_new
			getExp(affine,affine_new,j*integration_step,10);

			// Apply affine transformation to sample point
			multMatrix_4x4_x_3f(affine_new, sample, pos2); // Assumes bottom row of matrix is [0,0,0,1]

			if(paths_name != NULL)
				fprintf(f_paths,"draw sphere \"%8.5f %8.5f %8.5f\" radius %.4f resolution 10\n",pos2[0],pos2[1],pos2[2],0.2);

			// Check distance between current point and all molecule atoms
			dist_min = 9999999.9; // reset minimal distance between trajectory point and cluster atoms
			for(int k=0; k<num_atoms && !collision; k++)
			{
				dist = sqrDist(pos2,coords+k*3);
				if(dist < close2path_thr2) // occlusion test...
					collision = true;
				if(atom_owner[k] == c && dist < dist_min)
					dist_min = dist; // Update the minimal distance with cluster atoms
			}

			if(!collision && dist_min > far2path_thr2)// if current trajectory point is too far from the nearest atom of current cluster
				collision = true; // it does not collide, but should not be considered

			// Check distance between current point and the initial one
			for(int k=pathgap; k < npathpoints-pathgap && !collision; k++)
			{
				dist = sqrDist(pos2,pathpoints+3*k);
				if(dist < close2path_thr2)
					collision = true;
			}

			if(!collision) // Adding a valid path point
			{
				if(npathpoints < maxpathpoints) // This should avoid any possible memory overflow...
				{
					// Current path length
					pathLength += Dist(pos1,pos2);

					jmin = j; // this way it stores the last valid path amplitude
					pathpoints[3*npathpoints] = pos2[0];
					pathpoints[3*npathpoints + 1] = pos2[1];
					pathpoints[3*npathpoints + 2] = pos2[2];
					npathpoints++;
				}
				else
					collision = true; // stop adding points
			}

			// Swapping positions
			buffer = pos1;
			pos1 = pos2;
			pos2 = buffer;
		}

		if(npathpoints > 0)
		{
			// Computing the maximum inter-trajectory-point distance
			score = 0.0;
			for(int m=0; m < npathpoints; m++)
				for(int n=m+1; n < npathpoints; n++)
				{
					dist = sqrDist(pathpoints+3*m,pathpoints+3*n);
					if(score < dist)
						score = dist;
				}
		}
		else
		{
			pathLength = 0.0;
			dist = 0.0;
			score = 0.0;
			jmin = 0.0;
			jmax = 0.0;
		}

		if(debug)
			printf("Sample %4d  dist= %f  length= %f perUnit= %f  jmax= %f  jmin= %f  npathpoints= %d  --> score= %f\n", i, dist, pathLength, lengthPerUnit[i], jmax, jmin, npathpoints, score);

		if(score > score_best)
		{
			i_best = i;
			jmax_best = jmax;
			jmin_best = jmin;
			pathLength_best = pathLength;
			score_best = score;
		}
	}

	if(paths_name != NULL)
	{
		fprintf(f_paths,"display resetview\n");
		fclose(f_paths);
	}

	// Some checking...
	if(i_best < 0) // No best points found!
	{
		printf("ERROR: No valid points exists!  nsamples= %d  i_best= %d  score_best= %f\n", nsamples, i_best, score_best);
		exit(2);
	}

	// Output
	*p_jmax_best = jmax_best;
	*p_jmin_best = jmin_best;
	*p_pathLength_best = pathLength_best;
	*p_i_best = i_best;

	// Free memory
	free(pathpoints);
	free(pos1);
	free(pos2);

}


// Recursive definition of determinate using expansion by minors.
double Determinant(double **a,int n)
{
	int i,j,j1,j2;
	double det = 0;
	double **m = NULL;

	if (n < 1)
	{ /* Error */
		exit(2);
	}
	else if (n == 1)
	{ /* Shouldn’t get used */
		det = a[0][0];
	}
	else if (n == 2)
	{
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	}
	else
	{
		det = 0;
		for (j1=0;j1<n;j1++)
		{
			m = (double **) malloc((n-1)*sizeof(double *));
			for (i=0;i<n-1;i++)
				m[i] = (double *) malloc((n-1)*sizeof(double));
			for (i=1;i<n;i++)
			{
				j2 = 0;
				for (j=0;j<n;j++)
				{
					if (j == j1)
						continue;
					m[i-1][j2] = a[i][j];
					j2++;
				}
			}
			det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
			for (i=0;i<n-1;i++)
				free(m[i]);
			free(m);
		}
	}
	return(det);
}

// Draws arrow/s and path from affine model
//
// format --> 1: VMD, 2: JMol
// sample --> Best valid sample coordinates pointer
// affine --> Affine model matrix (4x4)
// jmin_best --> Initial path point relative to affine model)
// jmax_best --> Final path point (relative to affine model)
// pathLength_best --> Best path length
// lengthPerUnit --> Length per unit for best valid sample point
// comf --> Cluster's Center of Mass (weighted by vector field)
// delta_path --> Arrow sampling (in Angstroms)
// integration_step --> Affine model integration step
// f_name --> Best arrow & path file name (without extension, it will be aded according to "format")
// arrow_number --> Number of arrows
// arrow_length --> Arrow length
// ribbon_width --> Ribbon width in Angstroms
// arrow_width --> Arrow width in Angstroms
// tip_width --> Arrow's Tip width (in Angstroms)
// tip_length --> Arrow's Tip length (in Angstroms)
// arrow_aspect --> // Arrow aspect ratio (non_arrow_tip / arrow_length) <-- (It will be redefined later if "tip_length" > 0)
// arrow_to_path --> Arrow to path ratio
// min_aspect --> Minimal aspect ratio when arrow length is smaller than tip length
// arrows_apart --> Separation (in Angstroms) between arrows and paths
// path_middle_lines --> Enable or disable path middle lines
void drawArrow(int format, float *sample, double *affine, float jmin_best, float jmax_best, float pathLength_best, float lengthPerUnit, float *comf,
		float delta_path, float integration_step, char *f_name, int arrow_number, float arrow_length, float ribbon_width, float arrow_width, float tip_width,
		float tip_length, float arrow_aspect, float arrow_to_path, float min_aspect, float arrows_apart, bool path_middle_lines)
{
	bool debug = false;
	bool triangles = true; // true --> Draws arrow with triangles, false --> with quads (Triangles look much better than quads in Jmol... why?)
	float tip_increment = delta_path * (tip_width)/(arrow_length*(1.0-arrow_aspect)); // Arrow's tip "vertical" increment between steps
	int ntip = 0; // tip step index
	int arrow_index = 0; // Current arrow index
	float *arrow_start; // Array with the length at which every arrow starts
	float *arrow_tip; // Array with the length at which the arrow-tip starts
	float *arrow_end; // Array with the length at which every arrow ends
	float tip1[3],tip2[3];
	float dummy[3],dummy2[3],dummy3[3];
	float p1[3],p2[3];
	bool tip = true;
	bool overlap = true; // Enable arrow/path overlapping
	bool arrow_starting; // Signaling arrow start
	double affine_new[16];
	float pathLength; // Current path length
	float *buffer;
	float *norm1,*norm2; // trajectory and "radius" normal vectors
	float *a1,*a2; // arrow moved from cluster center points
	float *pos1,*pos2;
	arrow_start = (float *) allocate(sizeof(float)*arrow_number, "");
	arrow_tip = (float *) allocate(sizeof(float)*arrow_number, "");
	arrow_end = (float *) allocate(sizeof(float)*arrow_number, "");
	norm1 = (float *) allocate(sizeof(float)*3, "");
	norm2 = (float *) allocate(sizeof(float)*3, "");
	a1 = (float *) allocate(sizeof(float)*3, "");
	a2 = (float *) allocate(sizeof(float)*3, "");
	pos1 = (float *) allocate(sizeof(float) * 3, "");
	pos2 = (float *) allocate(sizeof(float) * 3, "");

	// Compute arrow limits for "n" arrows of length "l" spaced by "a" and laying in a trajectory of length "L"
	//     a     l      a     l      a      l      a
	// |-------=====>-------=====>-------=====>-------|                  a = (L-n*l)/(n+1)
	//
	float a = (pathLength_best-arrow_number*arrow_length) / (arrow_number+1);
	for(int i=0; i < arrow_number; i++)
	{
		arrow_start[i] = a*(i+1) + arrow_length*i;
		arrow_tip[i] = a*(i+1) + arrow_length*(i + arrow_aspect);
		arrow_end[i] = (a+arrow_length)*(i+1);
	}

	//	PMESH: Format of the pmesh files required by Jmol:
	//
	//	The format of a pmesh file is relatively simple (example file):
	//
	//	100
	//	3.0000 3.0000 1.0000
	//	2.3333 3.0000 1.0000
	//	...(98 more like this)
	//	81
	//	5
	//	0
	//	10
	//	11
	//	1
	//	0
	//	...(80 more sets like this)
	//
	//	The first line defines the number of grid points defining the surface (integer, n)
	//
	//	The next n lines define the Cartesian coordinates of each of the grid points (n lines of x, y, z floating point data points)
	//
	//	The next line specifies the number of polygons, m, to be drawn (81 in this case).
	//
	//	The next m sets of numbers, one number per line, define the polygons.
	//	In each set, the first number, p, specifies the number of points in each set.
	//	Currently this number must be either 4 (for triangles) or 5 (for quadrilaterals).
	//	The next p numbers specify indexes into the list of data points (starting with 0).
	//	The first and last of these numbers must be identical in order to "close" the polygon.


	// WARNING: check this later!!! There is needed more memory than it should!!!
	// It depends on the sampling difference between path (coarse) sampling and arrow (fine) sampling.
	int realloc_block = 100; // memory reallocation block size

	// For arrows...
	int npoints = 0; // number of grid points defining the surface (integer, n)
	int npolys = 0; // number of polygons, m, to be drawn
	int nsteps_max = (int) ceil(arrow_length/delta_path); // estimated number of arrow steps...
	int allocated_points = nsteps_max * 2; // two points per step
	int allocated_polys;
	float *points = (float *) allocate(sizeof(float) * allocated_points * 3,""); // three coordinates per point
	int *polys;
	if(triangles)
	{
		allocated_polys = nsteps_max * 2; // two triangles per arrow step
		polys = (int *) allocate(sizeof(int) * allocated_polys * 3,""); // three points per triangle
	}
	else
	{
		allocated_polys = nsteps_max; // one quad per arrow step
		polys = (int *) allocate(sizeof(int) * allocated_polys * 4,"");  // four points per triangle
	}

	//	printf("pathLength_best= %f  arrow_length= %f  delta_path= %f  npoints_max= %f\n",pathLength_best,arrow_length,delta_path,npoints_max);

	// For paths...
	int npoints2 = 0; // number of grid points defining the surface (integer, n)
	int npolys2 = 0; // number of polygons, m, to be drawn
	int nsteps2_max = (int) ceil(pathLength_best/delta_path); // estimated number of path steps...
	int allocated_points2 = nsteps2_max * 2; // two points per path step
	int allocated_polys2 = nsteps2_max; // one square per path step
	float *points2 = (float *) allocate(sizeof(float) * allocated_points2 * 3,""); // two points per path step...
	int *polys2 = (int *) allocate(sizeof(int) * allocated_polys2 * 4,""); // four points each square...

	if(debug)
		//		printf("pathLength_best= %f  arrow_length= %f  delta_path= %f  npoints_max= %d  npoints2_max= %d\n",pathLength_best,arrow_length,delta_path,npoints_max,npoints2_max);
		printf("pathLength_best= %f  arrow_length= %f  delta_path= %f  nsteps_max= %d  nsteps2_max= %d\n",pathLength_best,arrow_length,delta_path,nsteps_max,nsteps2_max);

	FILE *f_best;
	FILE *f_path;
	char name[100];

	switch(format)
	{
	case 1:
		sprintf(name,"%s.vmd",f_name);
		if( (f_best = fopen(name, "w")) == NULL)
		{
			printf("imodview> Best arrow file writing error!! \n");
			exit(1);
		}
		fprintf(f_best,"molecule new\ndraw color %s\n",colortext);
		break;
	case 2:
		// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
		sprintf(name,"%s.pmesh",f_name);
		if( (f_best = fopen(name, "w")) == NULL)
		{
			printf("imodview> Best arrow file writing error!! \n");
			exit(1);
		}
		sprintf(name,"%s_path.pmesh",f_name);
		if( (f_path = fopen(name, "w")) == NULL)
		{
			printf("imodview> Best path file writing error!! \n");
			exit(1);
		}
		break;
	case 3:
		sprintf(name,"%s.ngl",f_name);
		if( (f_best = fopen(name, "w")) == NULL)
		{
			printf("imodview> Best arrow file writing error!! \n");
			exit(1);
		}
		sprintf(name,"%s_path.ngl",f_name);
		if( (f_path = fopen(name, "w")) == NULL)
		{
			printf("imodview> Best path file writing error!! \n");
			exit(1);
		}
		break;
	}

	//	  fprintf(f_paths,"draw color %s\ndraw material Transparent\n","iceblue");
	//	  fprintf(f_paths,"draw color %s\n",colortext);
	pathLength = 0.0;
	float j;
	int i_max = (jmax_best - jmin_best) / (delta_path/lengthPerUnit); // Maximum number of steps
	//	for(float j=jmin_best; j <= jmax_best; j += delta_path/lengthPerUnit) // Integrating at 1A length steps
	//  for(float j=jmin_best; j <= 0; j += delta_path/lengthPerUnit[i_best]) // Integrating at 1A length steps
	for(int i=0; i <= i_max; i++)
	{
		// First, it makes sure that enough memory is available
		if(npoints >= allocated_points - 4) // if it is running out of memory
		{
			points = (float *) reallocate(points, sizeof(float) * (allocated_points + realloc_block) * 3, ""); // three coordinates per point
			allocated_points += realloc_block;
		}
		if(npoints2 >= allocated_points2 - 4) // if it is running out of memory
		{
			points2 = (float *) reallocate(points2, sizeof(float) * (allocated_points2 + realloc_block) * 3, ""); // three coordinates per point
			allocated_points2 += realloc_block;
		}
		if(npolys >= allocated_polys - 4)
		{
			if(triangles)
				polys = (int *) reallocate(polys, sizeof(int) * (allocated_polys + realloc_block) * 3, ""); // three points per triangle
			else
				polys = (int *) reallocate(polys, sizeof(int) * (allocated_polys + realloc_block) * 4, "");  // four points per triangle
			allocated_polys += realloc_block;
		}
		if(npolys2 >= allocated_polys2 - 4)
		{
			polys2 = (int *) reallocate(polys2, sizeof(int) * (allocated_polys2 + realloc_block) * 4, "");  // four points per triangle
			allocated_polys2 += realloc_block;
		}

		j = jmin_best + i*delta_path/lengthPerUnit; // current integration point

		// Integrate affine model into affine_new
		getExp(affine,affine_new,j*integration_step,10);

		// Apply affine transformation to sample point
		multMatrix_4x4_x_3f(affine_new, sample, pos2); // Assumes bottom row of matrix is [0,0,0,1]

		if(j == jmin_best)
		{
			// Integrate affine model into affine_new
			getExp(affine,affine_new,(jmin_best-delta_path/lengthPerUnit)*integration_step,10);

			// Apply affine transformation to sample point
			multMatrix_4x4_x_3f(affine_new, sample, pos1); // Assumes bottom row of matrix is [0,0,0,1]
		}

		subtract(comf,pos2,dummy); // dummy --> vector from trajectory point to the cluster center
		subtract(pos2,pos1,dummy2); // dummy2 --> tangent vector to trajectory
		crossProduct(dummy,dummy2,dummy3); // dummy3 --> vector orthogonal to both trajectory and "radius"
		normalize(dummy3,norm2); // dummy --> normalized vector orthogonal to both trajectory and "radius"

		// Arrows should be slightly apart from trajectory...
		normalize(dummy,dummy2);
		multByScalar(dummy2,-arrows_apart);

		a2[0] = pos2[0] + dummy2[0];
		a2[1] = pos2[1] + dummy2[1];
		a2[2] = pos2[2] + dummy2[2];

		//		if(j == jmin_best)
		if(i == 0)
		{
			p1[0] = pos2[0]+norm2[0]*ribbon_width;
			p1[1] = pos2[1]+norm2[1]*ribbon_width;
			p1[2] = pos2[2]+norm2[2]*ribbon_width;
			p2[0] = pos2[0]-norm2[0]*ribbon_width;
			p2[1] = pos2[1]-norm2[1]*ribbon_width;
			p2[2] = pos2[2]-norm2[2]*ribbon_width;
			points2[3*npoints2]     = p1[0];
			points2[3*npoints2 + 1] = p1[1];
			points2[3*npoints2 + 2] = p1[2];
			npoints2++;
			points2[3*npoints2]     = p2[0];
			points2[3*npoints2 + 1] = p2[1];
			points2[3*npoints2 + 2] = p2[2];
			npoints2++;
		}

		//		if(j>jmin_best)
		if(i>0)
		{
			pathLength += Dist(pos1,pos2); // update current path length

			// PLOT TRAJECTORY
			if(overlap || (pathLength < arrow_start[arrow_index] || pathLength > arrow_end[arrow_index]))
			{
				switch(format)
				{
				case 1:
					// line {x1 y1 z1} {x2 y2 z2} [width w] [style $<$solid|dashed$>$]
					fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
							pos1[0]+norm1[0]*ribbon_width,pos1[1]+norm1[1]*ribbon_width,pos1[2]+norm1[2]*ribbon_width,
							pos2[0]+norm2[0]*ribbon_width,pos2[1]+norm2[1]*ribbon_width,pos2[2]+norm2[2]*ribbon_width, 2);
					fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
							pos1[0]-norm1[0]*ribbon_width,pos1[1]-norm1[1]*ribbon_width,pos1[2]-norm1[2]*ribbon_width,
							pos2[0]-norm2[0]*ribbon_width,pos2[1]-norm2[1]*ribbon_width,pos2[2]-norm2[2]*ribbon_width, 2);
					if(path_middle_lines)
						fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
								pos2[0]+norm2[0]*ribbon_width,pos2[1]+norm2[1]*ribbon_width,pos2[2]+norm2[2]*ribbon_width,
								pos2[0]-norm2[0]*ribbon_width,pos2[1]-norm2[1]*ribbon_width,pos2[2]-norm2[2]*ribbon_width, 2);
					break;

				case 2:
				case 3:
					// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
					p1[0] = pos2[0]+norm2[0]*ribbon_width;
					p1[1] = pos2[1]+norm2[1]*ribbon_width;
					p1[2] = pos2[2]+norm2[2]*ribbon_width;
					p2[0] = pos2[0]-norm2[0]*ribbon_width;
					p2[1] = pos2[1]-norm2[1]*ribbon_width;
					p2[2] = pos2[2]-norm2[2]*ribbon_width;
					points2[3*npoints2]     = p1[0];
					points2[3*npoints2 + 1] = p1[1];
					points2[3*npoints2 + 2] = p1[2];
					npoints2++;
					points2[3*npoints2]     = p2[0];
					points2[3*npoints2 + 1] = p2[1];
					points2[3*npoints2 + 2] = p2[2];
					npoints2++;
					polys2[4*npolys2]     = npoints2 - 4;
					polys2[4*npolys2 + 1] = npoints2 - 3;
					polys2[4*npolys2 + 2] = npoints2 - 1;
					polys2[4*npolys2 + 3] = npoints2 - 2;
					npolys2++;
					break;
				}
			}
			dummy2[0] = 0.0;
			dummy2[1] = 0.0;
			dummy2[2] = 0.0;

			// PLOT ARROW/S
			if(arrow_index < arrow_number && pathLength > arrow_start[arrow_index])
				if(pathLength > arrow_tip[arrow_index])
				{
					if(tip_width-(ntip+1)*tip_increment < (delta_path/10) ) // Draw tip
					{
						switch(format)
						{
						case 1:
							// Draw tip triangle
							fprintf(f_best,"draw triangle {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f}\n",
									a1[0]+norm1[0]*(tip_width-ntip*tip_increment), a1[1]+norm1[1]*(tip_width-ntip*tip_increment), a1[2]+norm1[2]*(tip_width-ntip*tip_increment),
									a1[0]-norm1[0]*(tip_width-ntip*tip_increment), a1[1]-norm1[1]*(tip_width-ntip*tip_increment), a1[2]-norm1[2]*(tip_width-ntip*tip_increment),
									a2[0], a2[1], a2[2]);
							// Draw tip triangle border lines
							fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
									a1[0]+norm1[0]*(tip_width-ntip*tip_increment), a1[1]+norm1[1]*(tip_width-ntip*tip_increment), a1[2]+norm1[2]*(tip_width-ntip*tip_increment),
									a2[0], a2[1], a2[2], 2);
							fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
									a1[0]-norm1[0]*(tip_width-ntip*tip_increment), a1[1]-norm1[1]*(tip_width-ntip*tip_increment), a1[2]-norm1[2]*(tip_width-ntip*tip_increment),
									a2[0], a2[1], a2[2], 2);
							break;
						case 2:
						case 3:
							// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
							points[3*npoints]     = a2[0];
							points[3*npoints + 1] = a2[1];
							points[3*npoints + 2] = a2[2];
							npoints++;
							if(triangles)
							{
								polys[3*npolys]     = npoints - 3;
								polys[3*npolys + 1] = npoints - 2;
								polys[3*npolys + 2] = npoints - 1;
							}
							else
							{
								polys[4*npolys]     = npoints - 3;
								polys[4*npolys + 1] = npoints - 2;
								polys[4*npolys + 2] = npoints - 1;
							}
							npolys++;
							break;
						}
						arrow_index++; // next arrow
						arrow_starting = true; // signaling new arrow
						ntip=0; // reset tip step index
					}
					else // Draw tip body
					{
						if(ntip == 0) // then draw tip's base lines
						{
							p1[0] = a1[0]+norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0];
							p1[1] = a1[1]+norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1];
							p1[2] = a1[2]+norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2];
							p2[0] = a1[0]-norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0];
							p2[1] = a1[1]-norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1];
							p2[2] = a1[2]-norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2];
							switch(format)
							{
							case 1:
								// draw first base line
								fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
										a1[0]+norm1[0]*(arrow_width-ntip*tip_increment)+dummy2[0], a1[1]+norm1[1]*(arrow_width-ntip*tip_increment)+dummy2[1], a1[2]+norm1[2]*(arrow_width-ntip*tip_increment)+dummy2[2],
										a1[0]+norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]+norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]+norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2], 2);
								// draw second base line
								fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
										a1[0]-norm1[0]*(arrow_width-ntip*tip_increment)+dummy2[0], a1[1]-norm1[1]*(arrow_width-ntip*tip_increment)+dummy2[1], a1[2]-norm1[2]*(arrow_width-ntip*tip_increment)+dummy2[2],
										a1[0]-norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]-norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]-norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2], 2);
								break;
							case 2:
							case 3:
								// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
								points[3*npoints]     = p1[0];
								points[3*npoints + 1] = p1[1];
								points[3*npoints + 2] = p1[2];
								npoints++;
								points[3*npoints]     = p2[0];
								points[3*npoints + 1] = p2[1];
								points[3*npoints + 2] = p2[2];
								npoints++;
								break;
							}
						}
						p1[0] = a2[0]+norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0];
						p1[1] = a2[1]+norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1];
						p1[2] = a2[2]+norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2];
						p2[0] = a2[0]-norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0];
						p2[1] = a2[1]-norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1];
						p2[2] = a2[2]-norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2];
						switch(format)
						{
						case 1:
							// draw first triangle
							fprintf(f_best,"draw triangle {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f}\n",
									a1[0]+norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]+norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]+norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2],
									a1[0]-norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]-norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]-norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2],
									a2[0]+norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0], a2[1]+norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1], a2[2]+norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2]);
							// draw second triangle
							fprintf(f_best,"draw triangle {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f}\n",
									a2[0]+norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0], a2[1]+norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1], a2[2]+norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2],
									a2[0]-norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0], a2[1]-norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1], a2[2]-norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2],
									a1[0]-norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]-norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]-norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2]);
							// draw first border line
							fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
									a1[0]+norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]+norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]+norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2],
									a2[0]+norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0], a2[1]+norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1], a2[2]+norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2], 2);
							// draw second border line
							fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
									a2[0]-norm2[0]*(tip_width-(ntip+1)*tip_increment)+dummy2[0], a2[1]-norm2[1]*(tip_width-(ntip+1)*tip_increment)+dummy2[1], a2[2]-norm2[2]*(tip_width-(ntip+1)*tip_increment)+dummy2[2],
									a1[0]-norm1[0]*(tip_width-ntip*tip_increment)+dummy2[0], a1[1]-norm1[1]*(tip_width-ntip*tip_increment)+dummy2[1], a1[2]-norm1[2]*(tip_width-ntip*tip_increment)+dummy2[2], 2);
							break;
						case 2:
						case 3:
							// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
							if(triangles)
							{
								points[3*npoints]     = p1[0];
								points[3*npoints + 1] = p1[1];
								points[3*npoints + 2] = p1[2];
								polys[3*npolys]     = npoints - 2;
								polys[3*npolys + 1] = npoints - 1;
								polys[3*npolys + 2] = npoints;
								npoints++;
								npolys++;
								points[3*npoints]     = p2[0];
								points[3*npoints + 1] = p2[1];
								points[3*npoints + 2] = p2[2];
								polys[3*npolys]     = npoints - 2;
								polys[3*npolys + 1] = npoints - 1;
								polys[3*npolys + 2] = npoints;
								npoints++;
								npolys++;
							}
							else
							{
								points[3*npoints]     = p1[0];
								points[3*npoints + 1] = p1[1];
								points[3*npoints + 2] = p1[2];
								npoints++;
								points[3*npoints]     = p2[0];
								points[3*npoints + 1] = p2[1];
								points[3*npoints + 2] = p2[2];
								npoints++;
								polys[4*npolys]     = npoints - 4;
								polys[4*npolys + 1] = npoints - 3;
								polys[4*npolys + 2] = npoints - 1;
								polys[4*npolys + 3] = npoints - 2;
								npolys++;
							}
							break;
						}
						ntip++;
					}
				}
				else // Draw arrow body
				{
					if(arrow_starting)
					{
						p1[0] = a1[0]+norm1[0]*arrow_width+dummy2[0];
						p1[1] = a1[1]+norm1[1]*arrow_width+dummy2[1];
						p1[2] = a1[2]+norm1[2]*arrow_width+dummy2[2];
						p2[0] = a1[0]-norm1[0]*arrow_width+dummy2[0];
						p2[1] = a1[1]-norm1[1]*arrow_width+dummy2[1];
						p2[2] = a1[2]-norm1[2]*arrow_width+dummy2[2];
						switch(format)
						{
						case 1:
							// draw base line
							fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
									a1[0]+norm1[0]*arrow_width+dummy2[0], a1[1]+norm1[1]*arrow_width+dummy2[1], a1[2]+norm1[2]*arrow_width+dummy2[2],
									a1[0]-norm1[0]*arrow_width+dummy2[0], a1[1]-norm1[1]*arrow_width+dummy2[1], a1[2]-norm1[2]*arrow_width+dummy2[2], 2);
							break;
						case 2:
						case 3:
							// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
							points[3*npoints]     = p1[0];
							points[3*npoints + 1] = p1[1];
							points[3*npoints + 2] = p1[2];
							npoints++;
							points[3*npoints]     = p2[0];
							points[3*npoints + 1] = p2[1];
							points[3*npoints + 2] = p2[2];
							npoints++;
							break;
						}
						arrow_starting = false; // signaling new arrow
					}

					p1[0] = a2[0]+norm2[0]*arrow_width+dummy2[0];
					p1[1] = a2[1]+norm2[1]*arrow_width+dummy2[1];
					p1[2] = a2[2]+norm2[2]*arrow_width+dummy2[2];
					p2[0] = a2[0]-norm2[0]*arrow_width+dummy2[0];
					p2[1] = a2[1]-norm2[1]*arrow_width+dummy2[1];
					p2[2] = a2[2]-norm2[2]*arrow_width+dummy2[2];

					switch(format)
					{
					case 1:
						// draw first triangle
						fprintf(f_best,"draw triangle {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f}\n",
								a1[0]+norm1[0]*arrow_width+dummy2[0], a1[1]+norm1[1]*arrow_width+dummy2[1], a1[2]+norm1[2]*arrow_width+dummy2[2],
								a1[0]-norm1[0]*arrow_width+dummy2[0], a1[1]-norm1[1]*arrow_width+dummy2[1], a1[2]-norm1[2]*arrow_width+dummy2[2],
								a2[0]+norm2[0]*arrow_width+dummy2[0], a2[1]+norm2[1]*arrow_width+dummy2[1], a2[2]+norm2[2]*arrow_width+dummy2[2]);
						// draw second triangle
						fprintf(f_best,"draw triangle {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f}\n",
								a2[0]+norm2[0]*arrow_width+dummy2[0], a2[1]+norm2[1]*arrow_width+dummy2[1], a2[2]+norm2[2]*arrow_width+dummy2[2],
								a2[0]-norm2[0]*arrow_width+dummy2[0], a2[1]-norm2[1]*arrow_width+dummy2[1], a2[2]-norm2[2]*arrow_width+dummy2[2],
								a1[0]-norm1[0]*arrow_width+dummy2[0], a1[1]-norm1[1]*arrow_width+dummy2[1], a1[2]-norm1[2]*arrow_width+dummy2[2]);
						// draw first border line
						fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
								a1[0]+norm1[0]*arrow_width+dummy2[0], a1[1]+norm1[1]*arrow_width+dummy2[1], a1[2]+norm1[2]*arrow_width+dummy2[2],
								a2[0]+norm2[0]*arrow_width+dummy2[0], a2[1]+norm2[1]*arrow_width+dummy2[1], a2[2]+norm2[2]*arrow_width+dummy2[2], 2);
						// draw second border line
						fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
								a2[0]-norm2[0]*arrow_width+dummy2[0], a2[1]-norm2[1]*arrow_width+dummy2[1], a2[2]-norm2[2]*arrow_width+dummy2[2],
								a1[0]-norm1[0]*arrow_width+dummy2[0], a1[1]-norm1[1]*arrow_width+dummy2[1], a1[2]-norm1[2]*arrow_width+dummy2[2], 2);
						break;
					case 2:
					case 3:
						// fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
						if(triangles)
						{
							points[3*npoints]     = p1[0];
							points[3*npoints + 1] = p1[1];
							points[3*npoints + 2] = p1[2];
							polys[3*npolys]     = npoints - 2;
							polys[3*npolys + 1] = npoints - 1;
							polys[3*npolys + 2] = npoints;
							npoints++;
							npolys++;
							points[3*npoints]     = p2[0];
							points[3*npoints + 1] = p2[1];
							points[3*npoints + 2] = p2[2];
							polys[3*npolys]     = npoints - 2;
							polys[3*npolys + 1] = npoints - 1;
							polys[3*npolys + 2] = npoints;
							npoints++;
							npolys++;
						}
						else
						{
							points[3*npoints]     = p1[0];
							points[3*npoints + 1] = p1[1];
							points[3*npoints + 2] = p1[2];
							npoints++;
							points[3*npoints]     = p2[0];
							points[3*npoints + 1] = p2[1];
							points[3*npoints + 2] = p2[2];
							npoints++;
							polys[4*npolys]     = npoints - 4;
							polys[4*npolys + 1] = npoints - 3;
							polys[4*npolys + 2] = npoints - 1;
							polys[4*npolys + 3] = npoints - 2;
							npolys++;
						}
						break;
					}
				}
		}
		else
		{
			if(path_middle_lines)
				switch(format)
				{
				case 1:
					fprintf(f_best,"draw line {%8.5f %8.5f %8.5f} {%8.5f %8.5f %8.5f} width %d style solid\n",
							pos2[0]+norm2[0]*ribbon_width,pos2[1]+norm2[1]*ribbon_width,pos2[2]+norm2[2]*ribbon_width,
							pos2[0]-norm2[0]*ribbon_width,pos2[1]-norm2[1]*ribbon_width,pos2[2]-norm2[2]*ribbon_width, 2);
					break;
				case 2:
				case 3:
					fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
					break;
				}
		}

		// Swapping positions
		buffer = pos1;
		pos1 = pos2;
		pos2 = buffer;
		// Swapping positions
		buffer = a1;
		a1 = a2;
		a2 = buffer;
		// Swapping normals
		buffer = norm1;
		norm1 = norm2;
		norm2 = buffer;
	}
	//  printf("i_best= %d  pathLength= %f  pathLength_best= %f \n",i_best,pathLength,pathLength_best);
	if(debug)
		printf("pathLength= %f  pathLength_best= %f \n",pathLength,pathLength_best);

	switch(format)
	{
	case 1:
		fprintf(f_best,"display resetview\n");
		fclose(f_best);
		break;
	case 2:
		//		fprintf(stderr,"JMOL NOT IMPLEMENTED YET\n");
		if(npoints <= allocated_points)
		{



			fprintf(f_best,"%d\n",npoints);
			for(int i=0; i<npoints; i++)
				fprintf(f_best,"%.4f %.4f %.4f\n",points[3*i],points[3*i+1],points[3*i+2]);
			fprintf(f_best,"%d\n",npolys);
			if(triangles)
				for(int i=0; i<npolys; i++)
					fprintf(f_best,"%d\n%d\n%d\n%d\n%d\n",4,polys[3*i],polys[3*i+1],polys[3*i+2],polys[3*i]);
			else
			{
				for(int i=0; i<npolys-1; i++)
					fprintf(f_best,"%d\n%d\n%d\n%d\n%d\n%d\n",5,polys[4*i],polys[4*i+1],polys[4*i+2],polys[4*i+3],polys[4*i]);
				fprintf(f_best,"%d\n%d\n%d\n%d\n%d\n",4,polys[4*(npolys-1)],polys[4*(npolys-1)+1],polys[4*(npolys-1)+2],polys[4*(npolys-1)]);
				// Last is a triangle
			}
			// fprintf(stdout,"OK:  npoints_max= %d  npoints= %d  npolys= %d!\n",npoints_max, npoints, npolys);
		}
		else
		{
			fprintf(stdout,"Memory chungo happened:  allocated_points= %d  npoints= %d  npolys= %d!\n",allocated_points, npoints, npolys);
			exit(2);
		}
		fclose(f_best);

		if(npoints2 <= allocated_points2)
		{
			fprintf(f_path,"%d\n",npoints2);
			for(int i=0; i<npoints2; i++)
				fprintf(f_path,"%.4f %.4f %.4f\n",points2[3*i],points2[3*i+1],points2[3*i+2]);
			fprintf(f_path,"%d\n",npolys2);
			for(int i=0; i<npolys2; i++)
				fprintf(f_path,"%d\n%d\n%d\n%d\n%d\n%d\n",5,polys2[4*i],polys2[4*i+1],polys2[4*i+2],polys2[4*i+3],polys2[4*i]);
			// fprintf(stdout,"OK:  npoints2_max= %d  npoints2= %d  npolys2= %d!\n",npoints2_max, npoints2, npolys2);
		}
		else
		{
			fprintf(stdout,"Memory chungo happened:  allocated_points2= %d  npoints2= %d  npolys2= %d!\n",allocated_points, npoints2, npolys2);
			exit(2);
		}
		fclose(f_path);

		break;


	case 3:

		if(npoints <= allocated_points)
		{

			float col[3];

			col[0]=0.5; col[1]=0; col[2]=0; // C1 =0


			fprintf(f_best,"\n\nshape.addMesh(");

			fprintf(f_best,"[ %.4f , %.4f , %.4f ",points[polys[0]],points[polys[0]+1],points[polys[0]+2]);
			fprintf(f_best,", %.4f , %.4f , %.4f ",points[polys[0+1]*3],points[polys[0+1]*3+1],points[polys[0+1]*3+2]);
			fprintf(f_best,", %.4f , %.4f , %.4f ",points[polys[0+2]*3],points[polys[0+2]*3+1],points[polys[0+2]*3+2]);

			for(int i=1; i<npolys; i++) {
				fprintf(f_best,", %.4f , %.4f , %.4f ",points[polys[3*i+0]*3+0],points[polys[3*i+0]*3+1],points[polys[3*i+0]*3+2]);
				fprintf(f_best,", %.4f , %.4f , %.4f ",points[polys[3*i+1]*3+0],points[polys[3*i+1]*3+1],points[polys[3*i+1]*3+2]);
				fprintf(f_best,", %.4f , %.4f , %.4f ",points[polys[3*i+2]*3+0],points[polys[3*i+2]*3+1],points[polys[3*i+2]*3+2]);
			}
			fprintf(f_best,"] , \n");

			fprintf(f_best,"[ %.4f , %.4f , %.4f ",col[0],col[1],col[2]);
			fprintf(f_best,", %.4f , %.4f , %.4f ",col[0],col[1],col[2]);
			fprintf(f_best,", %.4f , %.4f , %.4f ",col[0],col[1],col[2]);

			for(int i=1; i<npolys; i++) {
				fprintf(f_best,", %.4f , %.4f , %.4f ",col[0],col[1],col[2]);
				fprintf(f_best,", %.4f , %.4f , %.4f ",col[0],col[1],col[2]);
				fprintf(f_best,", %.4f , %.4f , %.4f ",col[0],col[1],col[2]);
			}
			fprintf(f_best,"] , \n\n");





		}

		fclose(f_best);


		fclose(f_path);

		break;













	}

	free(arrow_start);
	free(arrow_tip);
	free(arrow_end);
	free(norm1);
	free(norm2);
	free(a1);
	free(a2);
	free(pos1);
	free(pos2);
	free(points);
	free(polys);
}


void getMoreOwners(Macromolecule *mol, int **p_res_owner, char **p_chain_owner)
{
	bool debug = false;
	pdbIter *it_chain = new pdbIter(mol);
	pdbIter *it_frag;
	pdbIter *it_at;
	Atom *at;
	Fragment *frag;
	Chain *chain;
	int natoms = 0;
	int num_atoms = mol->get_num_atoms();
	int *res_owner; // Array telling which residue an atom belongs to
	char *chain_owner; // Array telling which chain an atom belongs to

	// Memory allocation
	if(*p_res_owner == NULL)
		res_owner = (int *) allocate( sizeof(int) * num_atoms, "" );
	if(*p_chain_owner == NULL)
		chain_owner = (char *) allocate( sizeof(char) * num_atoms, "" );
	*p_res_owner = res_owner; // output
	*p_chain_owner = chain_owner; // output

	for( it_chain->pos_chain = 0; !it_chain->gend_chain(); it_chain->next_chain() ) // screen chains
	{
		chain = it_chain->get_chain();
		it_frag = new pdbIter(chain);
		for( it_frag->pos_fragment = 0; !it_frag->gend_fragment(); it_frag->next_fragment() ) // screen fragments
		{
			frag = it_frag->get_fragment();
			it_at = new pdbIter(frag);
			for( it_at->pos_atom = 0; !it_at->gend_atom(); it_at->next_atom() ) // screen fragment's atoms
			{
				at = it_at->get_atom(); // get current atom
				if(debug)
					printf("atom %d belongs to fragment %d and chain %s\n", natoms, frag->getIdNumber(), chain->getName());
				res_owner[natoms] = frag->getIdNumber();
				chain_owner[natoms] = (chain->getName())[0];
				natoms++; // counts atoms
			}
		}
	}
}

void writeCluster(int *cluster, int natoms, int *res_owner, char *chain_owner, char *text)
{
	int ires;
	int ires_old = res_owner[ cluster[0] ];
	int ires0 = res_owner[ cluster[0] ];
	char ichain;
	char ichain_old = chain_owner[ cluster[0] ];

	FILE *f_file;
	if( (f_file = fopen(text, "w")) == NULL)
	{
		printf("File writing error!! \n");
		exit(1);
	}
	fprintf(f_file,"select ");

	for(int i=0; i<natoms; i++) // screen cluster atoms
	{
		ires = res_owner[ cluster[i] ];
		ichain = chain_owner[ cluster[i] ];

		if( !(ires == ires_old || ires == ires_old+1) || ichain != ichain_old ) // If not different residue or consecutive, or different chain
		{   // , then different segment or single residue.
			if(ires0 == ires_old)
				fprintf(f_file,"%d:%c,", ires0, ichain_old);
			else
				fprintf(f_file,"%d-%d:%c,", ires0, ires_old, ichain_old);
			ires0 = ires;
		}
		if(i == natoms-1) // last cluster atom
		{
			if(ires0 == ires) // if initial residue of current segment is equal to current residue, then single residue segment
				fprintf(f_file,"%d:%c", ires0, ichain_old);
			else // Otherwise current segment ends with current residue ("ires")
				fprintf(f_file,"%d-%d:%c", ires0, ires, ichain_old);
			ires0 = ires;
		}

		ires_old = ires;
		ichain_old = ichain;
	}

	fclose(f_file);
}

// Source code of simple quick sort implementation using array ascending order in c programming language
// Taken from: http://www.cquestions.com/2008/01/c-program-for-quick-sort.html
void quicksort(int *x, int first, int last)
{
	int pivot,j,temp,i;

	if(first<last)
	{
		pivot=first;
		i=first;
		j=last;

		while(i<j)
		{
			while( x[i] <= x[pivot] && i < last )
				i++;
			while( x[j] > x[pivot] )
				j--;
			if(i<j)
			{
				temp=x[i];
				x[i]=x[j];
				x[j]=temp;
			}
		}

		temp=x[pivot];
		x[pivot]=x[j];
		x[j]=temp;
		quicksort(x,first,j-1);
		quicksort(x,j+1,last);
	}
}
