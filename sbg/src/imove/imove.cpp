/*************************************************************************
 *                             iMOVE                                     *
 *************************************************************************
 * This program is part of iMOD: http://chaconlab.org/imod/index.html    *
 * (c) Jose Ramon Lopez-Blanco, Jose Ignacio Garzon and Pablo Chacon.    *
 * IQFR-CSIC's Structural Bioinformatics Group (2004-2011).              *
 *************************************************************************
 *                                                                       *
 *   It applies either ICS / CCS modes to a pdb model,                   *
 *   and outputs the deformed structures according to:                   *
 *   K-matrix, V/W-arrays, Simple rotation scheme, and linearly.         *
 *   (It takes into account Multiple-Chains and different CG-models)     *
 *                                                                       *
 *************************************************************************
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <stdio.h> // needed in 48Gb-ubuntuserver
#include "cmdl/CmdLine.h"
#include "libnma/include/libnma_misc.h" // Mon's NMA related library
#include "libnma/include/libnma_io.h" // Mon's NMA Input-Output library
#include "libnma/include/libnma_cg.h" // Mon's NMA Coarse-Graining library
#include "libnma/include/libnma_move.h" // Mon's NMA IC-Motion library
#include "libnma/include/libnma_deriv.h" // Mon's NMA Derivatives library
#include "libnma/include/libnma_time.h" // Real-timer (Santi's)

/*==============================================================================================*/
char version[]="1.14"; // version code
char file_pdb[FILE_NAME];
char file_ptraj[FILE_NAME];
char file_out[FILE_NAME];
char fix_file[FILE_NAME];
char proj_file[FILE_NAME];
char dummy_string[FILE_NAME];
char k_file[FILE_NAME];

float max_factor;
int type=0; // Chi/not-Chi dihedral angle
int model=2; // Coarse-Graining model
int mov_type; // Motion type
int typei = -1; // Internal coordinates type (output type of normal modes)
int modeli = -1; // Internal Coordinates model (output model of normal modes)
int ndiv,nev;
int nthreads = 1; // Number of threads used by move_vwMFA_par()
bool nomodel = false; // true = CG-model building and formating disabled (initial model)
bool saveformat_switch = false; // = true --> save formated input PDB
bool savemodel_switch = false; // = true --> save 3BB2R-model input PDB
bool fixFrag_switch = false;
bool fixIC_switch = false;
bool read_fixDH_switch = false;
bool verb = false;
bool nomass_switch = false;
bool defamp_switch = true; // by default, the RMSD thermal amplitude will be set.
bool lineal_switch = false; // = true, then lineal motion; otherwise "sinusoidal" bounce.
bool proj_switch = false; // Enables the generation if a PC-filtered MD-trajectory
bool harmonic_energy = false; // Enables harmonic energy calculation from a "Kfile"
int fixmodel = 0; // =0 --> no ICs fixation
int initial = 0; // first frame index (from 0 to N-1)
int final = 99999999; // last frame
double Ta; // Temperature
double amp_factor = 5; // Number of standard deviations to define maximum motion
double factor = 1; // Lineal factor to control the motion amplitude easily

// Boltzmann constant = 1.3806503 � 10-23 m2 kg s-2 K-1
// SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM   1.9872065e-3L
// 	This is the gas constant R in (kcal/mol)/K.
// SimTK_BOLTZMANN_CONSTANT_KCAL_ANGSTROM   SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM
// 	Boltzmann's constant in Kcal-Angstrom units of (kcal/mol)/kelvin; same as R.
double Kb = 0.00198717; // Kboltz (in Angstroms)
#define PI 3.14159265358979323846
#define MAXLINE 200 // Maximum number of characters per readed line

// *****************************************************
using namespace TCLAP;
void parseOptions(int argc, char** argv);

// Computes the harmonic energy of "mol" by using the provided force constants matrix "decint"
double molenergy(Macromolecule *mol,float *coord,twid *decint,int nipa);

// *****************************************************
// *                        MAIN                       *
// *****************************************************
int main(int argc, char **argv)
{
	bool debug = false;
	int num_res,num_atoms,index,iframe,iframe_old,i;
	Htimer ht_timer; // timer
	timerReal t1; // Santi's timer

	printf("imove>\nimove> Welcome to the Internal Coordinates Motion Tool v%s\nimove>\n",version);

	// Parsing Input
	parseOptions(argc,argv);

	printf("imove> Input PDB-file %s\n", file_pdb);
	printf("imove> Multi-PDB output file %s\n", file_out);
	printf("imove> Number of frames to be generated %d\n", ndiv);
	printf("imove> Normal Mode choosen %d\n", nev);
	nev--; // NM indices range from 0 to nmodes-1

	// CA  Conditions
	Condition *calpha = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	calpha->add( " CA " );
	Conditions *calpha2 = new Conditions();
	calpha2->add( calpha );

	// NCAC  Conditions ( N-, CA-, C- selection)
	Condition *ncac = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	ncac->add( " N  " );
	ncac->add( " CA " );
	ncac->add( " C  " );
	Conditions *ncac2 = new Conditions();
	ncac2->add( ncac );

	// Initialize aminoacids and nucleotids
	init_aminoacids();

	// Reading input pdb
	printf( "imove> Reading PDB file %s\n",file_pdb );
	Macromolecule * molr = new Macromolecule( "pdb" );
	molr->readPDB(file_pdb );
	molr->info(stdout);

	// Needed to match different CG-models
	// Warning: A "finer"-grained model is always expected as input!!!
	Macromolecule *molini;
	if(modeli != -1) // if an exchange in the CG-model is demanded
		molini = new Macromolecule(molr);

	if(!nomodel)
	{
		// Formating PDB first
		printf( "imove> Formatting residues names\n" );
		molr->format_residues(false,model);
		if(saveformat_switch)
		{
			sprintf(dummy_string,"imove_format.pdb");
			molr->writePDB( dummy_string );
		}
	}

	Macromolecule *molx,*molNCACx;
	pdbIter *iter;
	pdbIter *iterA,*iter2;
	float imass;

	// Needed further to convert normal modes into the appropriate atomic model
	molini = new Macromolecule(molr); // Initial (readed and formatted) macromolecule copy

	// Setting Coarse-Graining model
	switch(model)
	{
	case 0: // CA-only + (NH and CO)-terminal model
		printf( "imove> Coarse-Graining model: CA-model\n");

		// N,CA,C selection
		molNCACx = molr->select( ncac2 ); // Creates a new molecule copied & selected from molr

		// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
		// Warning, atoms are not copied, they're just pointers to the original atoms.
		// setmass = true --> adds masses to Occupancy and Bfactor, otherwise not left unchanged.
		// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
		//                   (unit mass is divided between NH or CO and their CAs)
		// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
		//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
		if(!nomodel) // Makes model (if it's needed)
			molx = cg_CA( molNCACx, true, nomass_switch ); // masses will be computed (warning mol is created by reference to the molNCACx Macromol.)
		else
			molx = cg_CA( molNCACx, false, nomass_switch ); // input pdb masses will be used

		// Saving "NCAC" model
		if(savemodel_switch)
		{
			sprintf(dummy_string,"imove_ncacx.pdb");
			molNCACx->writePDB( dummy_string );
		}
		break;

	case 4: // CA-only
		printf( "imove> Coarse-Graining model: CA-only\n");
		// if CA selection
		molx = molr->select( calpha2 );
		break;

	case 3: // N,CA,C-model
		printf( "imove> Coarse-Graining model: N,CA,C-model\n");

		// N,CA,C selection
		molx = molr->select( ncac2 );
		if(!nomodel) // Sets masses (if it's needed)
			mass_NCAC( molx, nomass_switch );
		break;

	case 1: // 3BB2R model
		printf( "imove> Coarse-Graining model: 3BB2R\n");
		molx = molr;
		if(!nomodel) // Makes 3BB2R model (if it's needed)
		{
			// CREATES a 3BB2R reduced model
			//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
			//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
			//     There are a few exceptions: for Ala, Gly and Pro,
			//     and for the 1st and last residues.
			printf("imove> Creating 3BB2R reduced model:\n");
			cg_3BBR2(molx);
		}
		break;

	case 2: // Full-Atom
		printf( "imove> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
		molx = molr;
		if(!nomodel) // Sets masses (if it's needed)
			mass_FA(molx, nomass_switch); // sets masses
		break;
	}

	// Saving current model
	if(savemodel_switch)
	{
		sprintf(dummy_string,"imove_model.pdb");
		molx->writePDB( dummy_string, false ); // renumbers the PDB
	}

	float *coord=NULL;
	double energy=0.0; // will store computed energy
	twid *decint; // force constants matrix
	int nipa; // number of interacting pairs of atoms
	if(harmonic_energy)
	{
		printf ("imove> Getting coordinates single row (for molx)\n");
		molx->coordMatrix( &coord );

		// Reading contacts form a force constants file (Klist)
		printf("imove> Reading force constants file for harmonic energy computations: %s",k_file);

		// Read Force constants file (Kfile) allocating memory
		// Warning: if "coord" not provided, distances ".d" will not be updated!
		read_Kfile(&decint, &nipa, k_file, coord);
		printf(" (%d nipas readed) ",nipa);

		printf ("imove> Free coordinates single row\n");
//		free(coord);
	}

	// Reading ptraj
	double *eveci,*evali;
	int n_atoms,n_vects,ncompi;
	read_ptraj(file_ptraj,&eveci,&evali,&n_atoms,&n_vects,&ncompi);
	// Showing some "ptraj" info...
	printf("imove> Ptraj info: %s  vectors=%d  components=%d\n",file_ptraj,n_vects,ncompi);

	// Each move should have an amplitude related to its RMSD at temperature T;
	// the Go's paper about variance in DAS is:
	// Go, N. Biophys. Chem., 53 (1990) 105-112. "A theorem on amplitudes of
	// thermal atomic fluctuations in large molecules ... calculated by NMA".
	//
	// Given we move in DAS directly (without Jacobian), we need an expression for the
	// DAS normal mode coordinate "ti". Lets see ec.(25):
	// 		MSD = <ti^2> = Kb*T/Li   (where Li is the i-th eigenvalue)
	// 		RMSD(ti) = sqrt( Kb*T/L )
	//
	// Assuming Boltzmann distribution and thermal equilibrium we can determine the
	// population above (or below) any given energy at any temperature "T", lets do it:
	// Boltzmann distribution function is:
	// 		f(E) = exp( -E/KbT )/KbT   (normalized to a probability of 1)
	// Calculating its integral between E=0 and E=E' we get the fraction
	// of population below energy E' (Fb):
	// 		F(E) = -exp( -E/KbT )  ---->  Fb(E') = 1 - exp( -E'/KbT )
	// Thus, the fraction above E' (Fa) is:
	// 		Fa(E') = exp( -E'/KbT )
	// If E' = s*Kb*T , we have:
	// 		Fa(E'= s*Kb*T) = exp( -s )
	// In this case, we can compute the following table:
	//  s(amp_factor)     Fa(E')
	//  1                0.36788
	//  2                0.13534
	//  3                0.04979
	//  4                0.01832
	//  5                0.00673795 **
	//  6                0.00248752
	// 10                0.0000454
	// **By default, we choose s=5 to define the maximum feasible end to end motion.
	// (about 0.67% of the molecules would have this energy level)
	// Note that given the equipartition theorem, each system vibrational DoF
	// will have  <Evib> = KbT , kinetic-->KbT/2 and potential-->KbT/2

// THE MATH (21/9/2010):
//	The potential energy associated to a normal mode "i" is:
//
//	  Vi = 0.5 * eval_i * Qi^2                                 (1)
//
//	where, "eval_i" is the i-th eigenvalue and "Qi" is the normal
//	coordinate (just a scalar).	Note eval_i has force constant units.
//
//	The classical Equipartition theorem says that each system
//	vibrational degree of freedom will have an average energy about:
//
//	  <E> = Kb * T                                             (2)
//
//	  where, Kb is the Boltzmann constant and T the absolute temperature.
//	This way,  using eq.1 and 2, the maximum average amplitude (Qi_max)
//	will be:
//
//	  Qi_max = (2*Kb*T/eval_i)^0.5                            (3)
//
//	  Given the solutions of the generalized eigenvalue problem yield
//	harmonic solutions, we generate the movie by applying the simple
//	harmonic motion to the normal mode coordinate of interest:
//
//	  Qi(k) = Qi_max * sin( (k*PI)/(n_fram-1) )               (4)
//
//	  where, Qi(k) is the normal mode coordinate amplitude for a given
//	frame, k is the frame index (ranging from -(n_fram-1)/2 to
//	+(n_fram-1)/2), PI is 3.14..., and n_fram is the total number of
//	frames in the movie.
//
//	  With the normal mode amplitudes, we can compute the IC amplitudes
//	to be applied to each IC (y_a):
//	  y_a(k) =  Qi(k) * q_i                                   (5)
//
//	  where, q_i is the eigenvector (in IC space). Note the vector y_a(k)
//	represents the IC coordiante amplitudes to be applied in each movie
//	frame.
//
//	  Obtaining the cartesian coordinates is straightforward. One just
//	have to rotate/translate each IC component of the vector y_i in the
//	initial structure.

	max_factor = factor * sqrt( amp_factor * Kb * Ta / evali[nev] ) ; // not needed to divide by 2 because: increment = max/(ndiv-1)
	printf("imove> Maximum internal coordinate amplitude: %f\n", max_factor);
//	fprintf(stderr,"max_factor(thermal)= %f\n",max_factor);

	// Computing increment
	double increment = max_factor/(ndiv-1);
	double increment_old = 0.0;
	double step;

	// Deleting old output file
	FILE *file = fopen(file_out,"w");
	fclose(file);

	// Some variables declarations
	int size;
	Atom *atom;
	Segment *seg;
	Residue *res;
	pdbIter *iter_atoms;
	pdbIter *iter_frags;
	pdbIter *iter_seg;
	int i_atom;
	iter_seg = new pdbIter( molr, true, true, true, true ); // Iterator to screen segments

	// Initializing some internal coords. related residue properties
	// (# atoms, # units, # internal coordinates, etc...)
	tri *propsx;
	float matrix4[4][4];
	bool *addrotx=NULL; // Should be added ROTATIONs? (with fixing)
	int *unat; // needed

	switch(model)
	{
	case 0:
		properCA(molNCACx,&propsx,&unat);
		break;
	case 3:
		properCA(molx,&propsx,&unat);
		break;
	case 1:
	case 2:
		properMFA(molx,&propsx,&unat,type,model);
		break;
	}

	num_res = iter_seg->num_fragment();
	num_atoms = iter_seg->num_atom();
	printf( "imove> Selected model number of residues: %d\n", num_res );
//	printf( "imove> Selected model number of (pseudo)atoms: %d\n", num_atoms );

	// Counting number of dihedral and internal coordinates (number of degrees of freedom)
	if(mov_type == 3) // if cartesian motion
	{
		size = molx->num_atoms() * 3;
	}
	else
	{
		// Theoretic number of DoFs= 3*T + 3*R -6 + Dihedral
		// (Before fixing, i.e. the fix-file format DoFs...)
		int old_size,seg_atoms=0;
		size = 0;
		old_size = size; // temp (Dihedral ICs)
		for( iter_seg->pos_fragment = 0; !iter_seg->gend_fragment(); iter_seg->next_fragment() ) // screen residues
			size += propsx[iter_seg->pos_fragment].nan; // Each residue may have different number of dihedral-angles!
		printf( "imove> Number of Dihedral angles: %d\n", size );
		old_size = size;
		for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // screen segments
		{
			seg_atoms = (iter_seg->get_segment())->num_atoms();
			if(seg_atoms > 1)
				size += 6;
			else
				size += 3;
		}
		printf( "imove> Rotational/Translational ICs (Non-Eckart): %d\n", size-old_size );
		size -= 6; // Eckart conditions
		printf( "imove> Predicted number of ICs (Eckart): %d\n", size );
	}

	// Mobile Internal Coordinates selection variables
	bool *fixedx=NULL;
	int old_size;

	// Needed to match different CG-models
	// Warning: A "finer"-grained model is allways expected as input pdb!!!
	Macromolecule *moli,*molNCAC;
	floating *evec;
	if(modeli != -1)
	{
		// Setting IC-Coarse-Graining model
		switch(modeli)
		{
		case 0: // CA-only + (NH and CO)-terminal model
			printf( "imove> Coarse-Graining model: CA-IC\n");

			// N,CA,C selection
			molNCAC = molini->select( ncac2 );

			// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
			// Warning, atoms are not copied, they're just pointers to the original atoms.
			// setmass = true --> adds masses to Occupancy and Bfactor, otherwise not left unchanged.
			// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
			//                   (unit mass is divided between NH or CO and their CAs)
			// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
			//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
			moli = cg_CA( molNCAC, true, nomass_switch ); // masses will be computed
			break;

		case 3: // N,CA,C-model
			printf( "imove> Coarse-Graining model: N,CA,C-model\n");
			// N,CA,C selection
			moli = molini->select( ncac2 );
			mass_NCAC( moli, nomass_switch );
			break;
		case 1: // 3BB2R model
			printf( "imove> Coarse-Graining model: 3BB2R\n");
			moli = new Macromolecule( molini );
			printf("imove> Creating 3BB2R reduced model:\n");
			cg_3BBR2(moli);
			break;
		case 2: // Full-Atom
			printf( "imove> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
			moli = molini;
			mass_FA(moli); // sets masses
			break;
		}
		num_atoms = moli->get_num_atoms();
		printf( "imove> Selected model number of (pseudo)atoms: %d\n", num_atoms );

		tri *propsi;
		switch(modeli)
		{
		case 0:
		case 3:
			properCA(moli,&propsi,&unat);
			break;
		case 1:
		case 2:
			properMFA(moli,&propsi,&unat,typei,modeli);
			break;
		}

		// Saving "i" model
		if(savemodel_switch)
		{
			sprintf(dummy_string,"imove_i.pdb");
			moli->writePDB( dummy_string );
		}

		// Computing Total number of degrees of freedom
		int sizei = 0;
		pdbIter *iteri = new pdbIter( moli ); // iter to screen fragments (residues)
		for( iteri->pos_fragment = 0; !iteri->gend_fragment(); iteri->next_fragment() ) // screen residues
		{
			i=iteri->pos_fragment; // shorter
			// hessian matrix rank
			sizei += propsi[i].nan; // Each residue may have different number of dihedral-angles!
			// indices of first pseudo-atoms of each residue
			if ( i != 0 )
				propsi[i].k1 = propsi[i - 1].k1 + propsi[i].nat;
		}
		printf( "imove> Number of dihedral angles: %d\n", sizei );
		delete iteri;

		int n_seg=0;
		int old_size;
		old_size = sizei; // temp
		pdbIter *iter_seg = new pdbIter( moli ); // Iterator to screen segments
		n_seg = iter_seg->num_segment();
		n_seg--; // (n_seg-1)
		sizei += 6*n_seg; // (n_seg-1)*6 aditional degrees of freedom must be added
		delete(iter_seg);
		printf( "imove> Input CG-model Inter-segment coords: %d (Rot+Trans)\n", sizei-old_size );
		printf( "imove> Input CG-model Internal Coordinates: %d (sizei)\n", sizei );

		// FIXATION of VARIABLES
		// Mobile Internal Coordinates selection "i"-model
		// Fixation array is related to "i"-model, not to "current" one!
		// input: fixedi ---> output: fixedx
		bool *fixedi=NULL;
		if(fixmodel != 0)
		{
			if( !(fixedi = (bool *)malloc(sizeof(bool)*sizei) ) ) // Allocate memory
			{
				printf("Unable to allocate memory. Forcing exit!\n");
				exit(1);
			}
			old_size = sizei;
			switch(fixmodel)
			{
			case 2:
				sizei = read_fixIC(fix_file,moli,propsi,fixedi);
				break;
			case 3:
				sizei = read_fixDH(fix_file,moli,propsi,fixedi,typei,modeli);
				break;
			default:
				printf("Sorry, unknowk fixation method!\nForcing exit\n");
				exit(1);
				break;
			}
			printf("imove> Input CG-model Fixed Internal Coordinates: %d\n", old_size-sizei );
			printf("imove> Input CG-model Final Internal Coordinates (sizei) = %d\n\n",sizei);
		}

		// Creates two auxiliar arrays with segment properties (needed due to fixing):
		//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
		//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
		// (Allocates memory itself, if it's needed)
		bool *addrot; // Should be added ROTATIONs? (with fixing)
		int *effseg; // Effective segment indices
		sizei = seg_props(moli, propsi, fixedi, modeli, typei, &addrot, &effseg);
		printf( "imove> Number of ICs predicted by seg_props(): %d\n", sizei );

		// Some checking: Number of ptraj vector components must match internal coordinates found in PDB (size)
		if( ncompi != sizei )
		{
			printf("imove> Sorry!\nimove> Mismatch between input ptraj file (%d) and input coarse-graining model (%d).\nimove> Please, check input: pdb, ptraj or --type\n",ncompi,sizei);
			exit(1);
		}

		if(fixedi != NULL)
		{
			if( !(fixedx = (bool *)malloc(sizeof(bool)*size) ) )
			{
				printf("Unable to allocate memory. Forcing exit!\n");
				exit(1);
			}

			// Change fix-arrays (internal coordinates) from a "i" model into a "f" model.
			// WARNING: "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
			//          (no memory allocation is performed)
			// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
			// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
			// "i"-model should be allways "finer" than "f"-model
			if(model > modeli) // x-model > i-model  -->  "inverse" ("x" filled with "i")
				changefixIC(molx, propsx, propsi, fixedx, fixedi, model, type, modeli, typei, true);
			else if (model == modeli) // current-model == i-model --> lets see...
			{
				if(type >= typei) // x-type >= i-type --> "inverse" ("x" filled with "i")
					changefixIC(molx, propsx, propsi, fixedx, fixedi, model, type, modeli, typei, true);
				else // x-type < i-type --> "normal" ("x" filled with "i")
					changefixIC(moli, propsi, propsx, fixedi, fixedx, modeli, typei, model, type, false);
			}
			else // x-model < i-model  --> "normal" ("x" filled with "i")
				changefixIC(moli, propsi, propsx, fixedi, fixedx, modeli, typei, model, type, false);

			// Counting new DoFs (for output model)
			old_size=0;
			for(i=0; i<size; i++)
				if(!fixedx[i]) // if "f" DoF is fixedx
					old_size++;
			printf( "imove> Output CG-model fixed coordinates: %d \n", old_size );
			printf( "imove> Output CG-model previous mobile coordinates: %d \n", size );
			size -= old_size;
			printf( "imove> Output CG-model mobile coordinates (size): %d \n", size );
		}

		// "f" FLOATING precission eigenvectors
		printf("imove> Eigenvector matrix memory-size (floating)= %.3f Mb\n",((float)sizeof(floating)*size*n_vects)/1e6);
		if( !( evec = (floating *) malloc(size*n_vects * sizeof(floating)) ) )
		{
			printf("Sorry, unable to allocate Eigenvectors memory!!! Try a lower \"nevec\".\nForcing exit!!!\n\n");
			exit(1);
		}
		// "eigvec" initialization
		for(i=0; i<size*n_vects; i++)
			evec[i]=0.0;

		// Converting Input IC-eigenvectors model/type into current one
		// Change normal modes (internal coordinates) from a "i" model into a "f" model.
		// WARNING: "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
		// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
		// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
		// "i"-model should be allways "finer" than "f"-model
		if(model > modeli) // "inverse" we want: current-model <-- i-model
			changemodelIC(molr, propsx, propsi, evec, eveci, fixedx, fixedi, n_vects, size, sizei, model, type, modeli, typei, true);
		else if (model == modeli) // i-model == f-model  --> normal
		{
			if(type >= typei) // "inverse" we want: current-model <-- i-model
				changemodelIC(molr, propsx, propsi, evec, eveci, fixedx, fixedi, n_vects, size, sizei, model, type, modeli, typei, true);
			else // "normal" we want: i-model --> current-models
				changemodelIC(moli, propsi, propsx, eveci, evec, fixedi, fixedx, n_vects, sizei, size, modeli, typei, model, type, false);
		}
		else // current-model < i-model -->  // "normal" we want: i-model --> current-models
			changemodelIC(moli, propsi, propsx, eveci, evec, fixedi, fixedx, n_vects, sizei, size, modeli, typei, model, type, false);
	}
	else
	{
		// Eigenvectors will be the readed ones
		evec = eveci;

		if(fixmodel != 0)
		{
			if( !(fixedx = (bool *)malloc(sizeof(bool)*size) ) ) // Allocate memory
			{
				printf("Unable to allocate memory. Forcing exit!\n");
				exit(1);
			}
			old_size = size;
			switch(fixmodel)
			{
			case 2:
				size = read_fixIC(fix_file,molx,propsx,fixedx);
				break;
			case 3:
				size = read_fixDH(fix_file,molx,propsx,fixedx,type,model);
				break;
			default:
				printf("Sorry, unknown fixation method!\nForcing exit\n");
				exit(1);
				break;
			}
			printf("imove> Number of Fixed Internal Coordinates: %d\n", old_size-size );
			printf("imove> Final Internal Coordinates (size) = %d\n\n",size);
		}

		// Creates two auxiliar arrays with segment properties (needed due to fixing):
		//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
		//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
		// (Allocates memory itself, if it's needed)
		int *effsegx=NULL; // Effective segment indices
		if(mov_type != 3) // if not-cartesian motion
		{
			size = seg_props(molx, propsx, fixedx, model, type, &addrotx, &effsegx);
			printf( "imove> Number of ICs predicted by seg_props(): %d\n", size );
			free(effsegx);
		}

		// Some checking: Number of ptraj vector components must match internal coordinates found in PDB (size)
		if( ncompi != size )
		{
			printf("imove> Sorry!\n"
					"imove> Mismatch between input ptraj file (%d) and input CG-model (%d).\n"
					"imove> Please, check input: pdb, ptraj or --type\n",ncompi,size);
			exit(1);
		}
	}

	Macromolecule *dummy; // dummy model pdb
	Macromolecule *dummyCA; // dummy model pdb
	double *uu;
//	float *coord;
	trd *der = NULL; // forces memory allocation
	Tcoor pos;
	double ***V = NULL, ***W = NULL;
//	bool **body1 = NULL;
	int **body1 = NULL; // (num_atoms x 4) sized integer matrix to store body 1/2 boundary information
	int nframe = 0;
	int written_frames = 0;

	// SELECTING MOTION TYPE:
	//
	switch(mov_type)
	{
	case 3: // LINEAL MOTION
	{
		if(proj_switch) // ED-trajectory regeneration enabled...
		{
			FILE *f_proj;
			if( !(f_proj = fopen(proj_file, "r")) )
			{
				printf("Msg(imove): Unable open file to be read %s\nForcing exit!\n\n",proj_file);
				exit(1);
			}

			double *projs;
			if( !(projs = (double *) malloc( sizeof(double) * (nev+1) )) )
			{
				fprintf(stderr,"Error, memory allocation failed \nForcing exit!\n");
				exit(1);
			}

			double *vector;
			if( !(vector = (double *)malloc( sizeof(double) * size)) )
			{
				fprintf(stderr,"Unable to allocate vector memory...\nForcing exit!\n");
				exit(1);
			}

			char c;
			int ok = 1;
			while( ok > 0 ) //
			{
				// Jumps to the next line... (and skips first line)
				do
				  c = fgetc(f_proj);
				while (c != '\n');

				if( fscanf(f_proj,"%d",&nframe) == EOF ) // get frame index
					ok = 0;
				for(int i=0; i <= nev;i++)
					if( fscanf(f_proj,"%lf",projs+i) == EOF ) // get PC projection
						ok = 0;

				if( ((initial <= 0 && final <= 0) || (nframe >= initial && nframe <= final)) && ok > 0)
				{
					printf("\r");
					printf("imove> Projecting frame= %5d...",nframe);
//					printf("imove> Projecting frame= %5d, projections=",nframe);
//					for(int i=0; i <= nev;i++)
//						printf(" %f",projs[i]);

					// Reset Cartesian vector...
					for(int j=0; j<size; j++)
						vector[j] = 0.0;

					// Converting from PC to Cartesian..
					for(int i=0; i <= nev;i++)
						for(int j=0; j<size; j++)
							vector[j] += evec[i*size + j] * projs[i];

					// Moving in Cartesian coords. from input structure...
					dummy = new Macromolecule(molx); // dummy model pdb
					move_cart(dummy, vector, 1);
					// printf("imove> Added frame %5d to Multi-PDB\n",nframe);
					dummy->writeMPDB(file_out,nframe);
					written_frames++; // counting written frames
					delete dummy; // needed to free memory
				}
			}
			printf("\n");
			fclose(f_proj);
		}
		else // Standard mode animation...
		{
			index = nev * size; // Eigenvector index
//			molx->coordMatrix( &coord );

			double energy2=0.0,energy_old;
			for(int k = -ndiv/2; k <= ndiv/2; k++) // screen frames
			{
				energy_old = energy;
				// Moving (translating internal coordinates into cartesian space)
				if(k!=0)
				{
					// Each frame is now generated from the first one! (minimizes error)
					dummy = new Macromolecule(molx); // dummy model pdb
					move_cart(dummy, evec+index, k*increment);
					if(harmonic_energy)
					{
						// Computes the harmonic energy of "mol" by using the provided force constants matrix "decint"
						energy = molenergy(dummy,coord,decint,nipa);
						printf("imove> Adding frame %5d to Multi-PDB E= %10.3lf E2= %10.3f Eh= %8.3lf\n",k,energy,energy_old-energy,evali[nev]*pow(k*increment,2));
					}
					else
						printf("imove> Adding frame %5d to Multi-PDB\n",k);
					dummy->writeMPDB(file_out,k);
					delete dummy; // needed to free memory
				}
				else
				{
					if(harmonic_energy)
					{
						// Computes the harmonic energy of "mol" by using the provided force constants matrix "decint"
						energy = molenergy(molx,coord,decint,nipa);
						printf("imove> Adding frame %5d to Multi-PDB E= %10.3lf E2= %10.3f Eh= %8.3lf <-- 1st model (%s)\n",k,energy,energy_old-energy,evali[nev]*pow(k*increment,2),file_pdb);
					}
					else
						printf("imove> Adding frame %5d to Multi-PDB <-- 1st model (%s)\n",k,file_pdb);
					// Each frame is now generated from the first one! (minimizes error)
					molx->writeMPDB(file_out,k);
				}
			}
			free(coord);
		}
		break;
	}

	case 0:
	{
		printf("imove> Moving with Jacobian (K-matrix) (iterative linear approach)\n");
		fflush(stdout);

		Macromolecule **movie;
		movie = (Macromolecule **) malloc( sizeof(Macromolecule *) * ndiv); // frames array (Macromolecules)
		index = nev*size; // Eigenvector index

		// First, Centering Macromolecule
		// Computing the PDB's Center of Mass (CoM)
		double mtot,mta;
		double r[3];
		r[0] = r[1] = r[2] = 0.0f;
		mtot = 0.0;

		if(model == 0 || model == 3)
			iter_atoms = new pdbIter( molNCACx ); // current frame iterator (N,CA,C selection)
		else
			iter_atoms = new pdbIter( molx ); // current frame iterator

		for( iter_atoms->pos_atom = 0; !iter_atoms->gend_atom(); iter_atoms->next_atom() )  // screens all-atoms
		{
			atom = ( iter_atoms->get_atom() );
			mta = atom->getPdbocc(); // Load mass...
			mtot += mta;
			atom->getPosition(pos);
			/* Sum(mass*coord) before putting the CoM at 0 */
			r[0] += mta * pos[0];
			r[1] += mta * pos[1];
			r[2] += mta * pos[2];
		}
		r[0] /= mtot;
		r[1] /= mtot;
		r[2] /= mtot;
		printf( "Msg(imove): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, r[0], r[1], r[2] );
		// shift CM of pdb_model to 0,0,0
		for( iter_atoms->pos_atom = 0; !iter_atoms->gend_atom(); iter_atoms->next_atom() )  // screens all-atoms
		{
			atom = ( iter_atoms->get_atom() );
			atom->getPosition(pos);
			pos[0] -= r[0];
			pos[1] -= r[1];
			pos[2] -= r[2];
			atom->setPosition(pos);
		}
		delete iter_atoms;
		// mol->writePDB( "mol.pdb" );

		// Backwards motion
		iframe_old = ndiv/2;
		if(model == 0 || model == 3)
			// N,CA,C selection
			movie[iframe_old] = new Macromolecule(molNCACx); // current frame pdb
		else
			movie[iframe_old] = new Macromolecule(molx); // current frame pdb

		if(!lineal_switch) // "sinusoidal" bounce
			increment = 0.0;

		for(int k = -1; k >= -(ndiv-1)/2; k--) // screen "backwards" frames
		{								       // this minimizes error...
			if(!lineal_switch) // "sinusoidal" bounce
			{
				// Updating increment (sinusoidal motion)
				increment_old = increment;
				increment = 0.5 * sin( (k*PI)/(ndiv-1) ) * max_factor; // 0.5 factor needed in sinusoidal motion
				step = increment - increment_old;
//				fprintf(stderr,"k= %d  sin( %f )= %f  increment = %f  step= %f\n",k,(k*PI)/(ndiv-1),sin((k*PI)/(ndiv-1)), increment, step);
			}
			else
				step = -increment;

			iframe = k+ndiv/2; // frame index
			movie[iframe] = new Macromolecule(movie[iframe_old]); // current frame pdb
			// (iframe_old is updated at the loop end)

			// Getting coordinates single row (pseudo-atom model)
			movie[iframe_old]->coordMatrix( &coord );

			ht_timer.restart(); // timer
			switch(model)
			{
			case 0:
			case 3:
				// (N,CA,C-model)
				printf("imove>  ?) Computing CA-only K-matrix with dydqMCA3x() time=");
				dydqMCA3x(movie[iframe_old],coord, &der, size, fixedx);
				break;
			case 1:
			case 2:
				printf("imove>  ?) Computing Full-Atom/3BB2R K-matrix with dydqMFAx() time=");
				dydqMFAx(movie[iframe_old],coord, propsx, &der, type, model, size, fixedx, addrotx);
				break;
			}
			printf("%s\n",ht_timer.print_time_sec());

			// MOVING
			iter_atoms = new pdbIter( movie[iframe] ); // current frame iterator
			for(iter_atoms->pos_atom=0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
			{   // screen atoms
				i_atom = size * iter_atoms->pos_atom; // a bit shorter...
				atom = iter_atoms->get_atom();
				atom->getPosition(pos); // get position
				for(int l=0; l< size; l++) // screen Dihedrals
				{
					pos[0] += der[i_atom + l].x * evec[index + l] * step;
					pos[1] += der[i_atom + l].y * evec[index + l] * step;
					pos[2] += der[i_atom + l].z * evec[index + l] * step;
				}
				atom->setPosition(pos); // set new position
			}

			iframe_old = iframe;
			delete iter_atoms;
		}

		// Forward motion
		iframe_old = ndiv/2;

		if(!lineal_switch) // "sinusoidal" bounce
			increment = 0.0;
		for(int k = 1; k <= ndiv/2; k++) // screen "backwards" frames
		{								   // this minimizes error...
			if(!lineal_switch) // "sinusoidal" bounce
			{
				// Updating increment (sinusoidal motion)
				increment_old = increment;
				increment = 0.5 * sin( (k*PI)/(ndiv-1) ) * max_factor; // 0.5 factor needed in sinusoidal motion
				step = increment - increment_old;
//				fprintf(stderr,"k= %d  sin( %f )= %f  increment = %f  step= %f\n",k,(k*PI)/(ndiv-1),sin((k*PI)/(ndiv-1)), increment, step);
			}
			else
				step = increment;

			iframe = k+ndiv/2; // frame index
			movie[iframe] = new Macromolecule(movie[iframe_old]); // current frame pdb

			// Getting coordinates single row (pseudo-atom model)
			movie[iframe_old]->coordMatrix( &coord );

			ht_timer.restart(); // timer
			switch(model)
			{
			case 0:
			case 3:
				// (N,CA,C-model)
				printf("imove>  ?) Computing CA-only K-matrix with dydqMCA3x() time=");
				dydqMCA3x(movie[iframe_old],coord, &der, size, fixedx);
				break;
			case 1:
			case 2:
				printf("imove>  ?) Computing Full-Atom/3BB2R K-matrix with dydqMFAx() time=");
				dydqMFAx(movie[iframe_old],coord, propsx, &der, type, model, size, fixedx, addrotx);
				break;
			}
			printf("%s\n",ht_timer.print_time_sec());

			// MOVING
			iter_atoms = new pdbIter( movie[iframe] ); // current frame iterator
			for(iter_atoms->pos_atom=0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
			{ // screen atoms
				i_atom = size * iter_atoms->pos_atom; // a bit shorter...
				atom = iter_atoms->get_atom();
				atom->getPosition(pos); // get position
				for(int l=0; l< size; l++) // screen Dihedrals
				{
					pos[0] += der[i_atom + l].x * evec[index + l] * step;
					pos[1] += der[i_atom + l].y * evec[index + l] * step;
					pos[2] += der[i_atom + l].z * evec[index + l] * step;
				}
				atom->setPosition(pos); // set new position
			}

			iframe_old = iframe;
			delete iter_atoms;
		}

		// Saving Movie
		for(int k = 0; k < ndiv; k++) // screen frames
		{
			// If CA-model
			if(model == 0)
			{
				dummy = movie[k];
				movie[k] = dummy->select_cpy( calpha2 ); // CA selection-copy
				delete dummy; // needed to free memory
			}

			// Placing back each frame into the original position
			iter_atoms = new pdbIter( movie[k] ); // current frame iterator
			for(iter_atoms->pos_atom=0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
			{ // screen atoms
				atom = iter_atoms->get_atom();
				atom->getPosition(pos); // get position
				pos[0] += r[0];
				pos[1] += r[1];
				pos[2] += r[2];
				atom->setPosition(pos); // set new position
			}
			delete iter_atoms;

			// Saving multi-pdb
			if(k!=ndiv/2)
				printf("imove> Added frame %5d to Multi-PDB\n",k);
			else
				printf("imove> Added frame %5d to Multi-PDB <-- 1st model (%s)\n",k,file_pdb);
			movie[k]->writeMPDB(file_out,k);
			delete movie[k];
		}
		break;
	}
	case 1: // Moving with V/W arrays
	{
		printf("imove> Moving with Jacobian computed \"in-situ\" from V/W-arrays (memory efficient)\n");
		fflush(stdout);

		Macromolecule **movie;
		movie = (Macromolecule **) malloc( sizeof(Macromolecule *) * ndiv); // frames array (Macromolecules)
		index = nev*size; // Eigenvector index

		double mtot,mta;
		double r[3],deriv[3],rk[3];
		r[0] = r[1] = r[2] = 0.0;
		mtot = 0.0;

		if(model == 0 || model == 3)
			iter_atoms = new pdbIter( molNCACx ); // current frame iterator (N,CA,C selection)
		else
			iter_atoms = new pdbIter( molx ); // current frame iterator

		// First, Centering Macromolecule
		// Computing the PDB's Center of Mass (CoM)
		for( iter_atoms->pos_atom = 0; !iter_atoms->gend_atom(); iter_atoms->next_atom() )  // screens all-atoms
		{
			atom = ( iter_atoms->get_atom() );
			mta = atom->getPdbocc(); // Load mass...
			mtot += mta;
			atom->getPosition(pos);
			/* Sum(mass*coord) before putting the CoM at 0 */
			r[0] += mta * pos[0];
			r[1] += mta * pos[1];
			r[2] += mta * pos[2];
		}
		r[0] /= mtot;
		r[1] /= mtot;
		r[2] /= mtot;
		printf( "imove> Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, r[0], r[1], r[2] );
		// shift CM of pdb_model to 0,0,0
		for( iter_atoms->pos_atom = 0; !iter_atoms->gend_atom(); iter_atoms->next_atom() )  // screens all-atoms
		{
			atom = ( iter_atoms->get_atom() );
			atom->getPosition(pos);
			pos[0] -= r[0];
			pos[1] -= r[1];
			pos[2] -= r[2];
			atom->setPosition(pos);
		}
		delete iter_atoms;

		if(nthreads != 0) // NOT parallel
			printf("imove> Parallelization enabled, using: %d threads\n",nthreads);
		else
			printf("imove> Parallelization disabled.\n");

		// Backwards motion
		iframe_old = ndiv/2;
		if(model == 0 || model == 3)
		{
			// N,CA,C selection
			movie[iframe_old] = new Macromolecule(molNCACx); // current frame pdb
		}
		else
			movie[iframe_old] = new Macromolecule(molx); // current frame pdb

		if(!lineal_switch) // "sinusoidal" bounce
			increment = 0.0;

		for(int k = -1; k >= -ndiv/2; k--) // screen "backwards" frames
		{								   // this minimizes error...
			if(!lineal_switch) // "sinusoidal" bounce
			{
				// Updating increment (sinusoidal motion)
				increment_old = increment;
				increment = 0.5 * sin( (k*PI)/(ndiv-1) ) * max_factor; // 0.5 factor needed in sinusoidal motion
				step = increment - increment_old;
//				fprintf(stderr,"k= %d  sin( %f )= %f  increment = %f  step= %f\n",k,(k*PI)/(ndiv-1),sin((k*PI)/(ndiv-1)), increment, step);
			}
			else
				step = -increment;

			iframe = k+ndiv/2; // frame index
			movie[iframe] = new Macromolecule(movie[iframe_old]); // current frame pdb

			// Getting coordinates single row (pseudo-atom model)
			movie[iframe_old]->coordMatrix( &coord );

			ht_timer.restart(); // timer
			switch(model)
			{
			case 0:
			case 3:
				printf("imove>  %2d Computing CA-only V/W-arrays with vwMCA3x() t=",iframe);
				fflush(stdout);
				vwMCA3x(movie[iframe_old],coord,&V,&W,&body1,size,fixedx);
				break;
			case 1:
			case 2:
				printf("imove>  %2d Computing HA/C5 V/W-arrays with vwMFAx() t=",iframe);
				fflush(stdout);
				vwMFAx(movie[iframe_old],coord,propsx,&V,&W,&body1,type,model,size,fixedx,addrotx);
				break;
			}
			printf("%s  mov=",ht_timer.print_time_sec());
			fflush(stdout);

			// MOVING
			ht_timer.restart(); // timer
			t1.startTimer();

			if(nthreads == 0) // NOT parallel
				move_VW(movie[iframe], step, evec, index, size, V, W, body1, model );
			else
				move_VWpar(movie[iframe], step, evec, index, size, V, W, body1, nthreads, model);

			printf("%s ",ht_timer.print_time_sec());
			t1.stopTimer();
			printf("(realtime= %lf s)\n",t1.getElapsedTime()); // Real elapsed time in seconds
			fflush(stdout);

			iframe_old = iframe;
		}

		// Forward motion
		if(!lineal_switch) // "sinusoidal" bounce
			increment = 0.0;

		iframe_old = ndiv/2;
		for(int k = 1; k <= ndiv/2; k++) // screen "backwards" frames
		{								   // this minimizes error...
			if(!lineal_switch) // "sinusoidal" bounce
			{
				// Updating increment (sinusoidal motion)
				increment_old = increment;
				increment = 0.5 * sin( (k*PI)/(ndiv-1) ) * max_factor; // 0.5 factor needed in sinusoidal motion
				step = increment - increment_old;
//				fprintf(stderr,"k= %d  sin( %f )= %f  increment = %f  step= %f\n",k,(k*PI)/(ndiv-1),sin((k*PI)/(ndiv-1)), increment, step);
			}
			else
				step = increment;

			iframe = k+ndiv/2; // frame index
			movie[iframe] = new Macromolecule(movie[iframe_old]); // current frame pdb

			// Getting coordinates single row (pseudo-atom model)
			movie[iframe_old]->coordMatrix( &coord );

			ht_timer.restart(); // timer
			switch(model)
			{
			case 0:
			case 3:
				printf("imove>  %2d Computing CA-only V/W-arrays with vwMCA3x() t=",iframe);
				fflush(stdout);
				vwMCA3x(movie[iframe_old],coord,&V,&W,&body1,size,fixedx);
				break;
			case 1:
			case 2:
				printf("imove>  %2d Computing HA/C5 V/W-arrays with vwMFAx() t=",iframe);
				fflush(stdout);
				vwMFAx(movie[iframe_old],coord,propsx,&V,&W,&body1,type,model,size,fixedx,addrotx);
				break;
			}
			printf("%s  mov=",ht_timer.print_time_sec());
			fflush(stdout);

			// MOVING
			ht_timer.restart(); // timer
			t1.startTimer();

			if(nthreads == 0) // NOT parallel
				move_VW(movie[iframe], step, evec, index, size, V, W, body1, model );
			else
				move_VWpar(movie[iframe], step, evec, index, size, V, W, body1, nthreads, model);

			printf("%s ",ht_timer.print_time_sec());
			t1.stopTimer();
			printf("(realtime= %lf s)\n",t1.getElapsedTime()); // Real elapsed time in seconds
			fflush(stdout);

			iframe_old = iframe;
		}

		// Saving Movie
		for(int k = 0; k < ndiv; k++) // screen frames
		{
			// If CA-model
			if(model == 0)
			{
				dummy = movie[k];
				movie[k] = dummy->select_cpy( calpha2 ); // CA selection-copy
				delete dummy; // needed to free memory
			}

			// Placing back each frame into the original position (removing CoM centering)
			iter_atoms = new pdbIter( movie[k] ); // current frame iterator
			for(iter_atoms->pos_atom=0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
			{ // screen atoms
				atom = iter_atoms->get_atom();
				atom->getPosition(pos); // get position
				pos[0] += r[0];
				pos[1] += r[1];
				pos[2] += r[2];
				atom->setPosition(pos); // set new position
			}
			delete iter_atoms;

			// Saving multi-pdb
			if(k!=ndiv/2)
				printf("imove> Added frame %5d to Multi-PDB\n",k);
			else
				printf("imove> Added frame %5d to Multi-PDB <-- 1st model (%s)\n",k,file_pdb);
			movie[k]->writeMPDB(file_out,k);
			delete movie[k];
		}
		break;
	}
	case 2:
	{
		printf("imove> Moving analytically (moving internal coordinates)\n");
		fflush(stdout);

		// Selected eigenvector memory allocation
		if( !(uu = (double *) malloc( sizeof(double) * size )) )
		{
			printf("imove> Memory allocation error!\n");
			exit(1);
		}

		// Getting selected mode eigenvector components
		if(model == 3 || model == 0)
			//			dihedral_comps(molNCACx,evec,nev,size,propsx,&uu);
			dihedral_comps(evec,nev,size,&uu);
		else
			//			dihedral_comps(molx,evec,nev,size,propsx,&uu);
			dihedral_comps(evec,nev,size,&uu);

		M4Rot *matrix4_op;
		for(int k = -ndiv/2; k <= ndiv/2; k++) // screen frames
		{
			// Moving (translating internal coordinates into cartesian space)
			if(k!=0)
			{
				// "sinusoidal" motion
				if(!lineal_switch)
					step = 0.5 * sin( (k*PI)/(ndiv-1) ) * max_factor; // 0.5 factor needed in sinusoidal motion
				else
					step = k*increment;

				if(model == 3 || model == 0)
				{
					// Each frame is now generated from the first one! (minimizes error)
					dummy = new Macromolecule(molNCACx); // dummy model pdb
					move_dihedralMCAx(dummy,uu, propsx, step, type, model, fixedx);
					molNCACx->minRmsd(dummy,matrix4);
				}
				else
				{
					// Each frame is now generated from the first one! (minimizes error)
					dummy = new Macromolecule(molx); // dummy model pdb
					move_dihedralMFAx(dummy, uu, propsx, step, type, model, fixedx, addrotx);
					molx->minRmsd(dummy,matrix4);
				}

				// The protein is re-built from N-terminal every time we move it,
				// so we need to align it respecting the readed PDB model. (It's fast and works nice)
				matrix4_op = new M4Rot(matrix4);
				dummy->applyAtoms(matrix4_op);
				delete(matrix4_op);
				printf("imove> Added frame %5d to Multi-PDB\n",k);
			}
			else
			{
				printf("imove> Added frame %5d to Multi-PDB <-- 1st model (%s)\n",k,file_pdb);
				if(model == 3 || model == 0)
					// Each frame is now generated from the first one! (minimizes error)
					dummy = new Macromolecule(molNCACx); // dummy model pdb
				else
					dummy = new Macromolecule(molx); // dummy model pdb
			}

			// If CA-model
			if(model == 0)
			{
				dummyCA = dummy->select_cpy( calpha2 ); // CA selection-copy
				dummyCA->writeMPDB(file_out,k);
			}
			else
				dummy->writeMPDB(file_out,k);
			delete dummy; // needed to free memory
		}
		break;
	}
	}

	if(proj_switch)
		printf("imove>\nimove> Success!\nimove> Written movie: %s (%d frames)\nimove> Bye!\n",file_out,written_frames);
	else
		printf("imove>\nimove> Success!\nimove> Written movie: %s (%d frames)\nimove> Bye!\n",file_out,ndiv);
}

/*==============================================================================================*/
void parseOptions(int argc, char** argv)
{
	std::string temp;
	CmdLine cmd("imove","Internal coordinates MOVEment tool.", version );

	try {
		// Labeled arguments definition
//		SwitchArg Verb("", "verb","Enables Hessian visualization and comparation with Naive version", true);
//		cmd.add( Verb );

        // DEVELOPER INPUT
        ValueArg<std::string> HarmonicEnergy("","energy", "Enables harmonic energy calculations from the provided force constants file (Kfile) (DEVELOPER's).",false,"fixstring","string");
        cmd.add( HarmonicEnergy );
        ValueArg<std::string> FixIC("","fixIC", "Plain-text file defining the fixed Internal Coordinates. "
	    		"Each line will contain the index (0,1,...) of the ICs to be removed "
	    		"(DEVELOPER's).",false,"fixstring","string");
        cmd.add( FixIC );
		SwitchArg ChiOut("","chi_out", "Considers first CHI dihedral angle in output models (default=disabled) (DEVELOPER's).", true);
		cmd.add( ChiOut );
        // END DEVELOPER INPUT

		SwitchArg NoModel("", "nomodel","Disables PDB model building. "
				"Warning: introduced PDB model must match the CG selected with the -m option (default=disabled).", false);
		cmd.add( NoModel );

        ValueArg<int> NThreads("","nthreads", "Number of threads for parallel processing (only used if --mov 1)",false,1,"int");
		cmd.add( NThreads );

		ValueArg<int> FrameF("","final", "Final frame index as defined in proj file (default= last) (only active if --proj option is enabled)",false,99999999,"int");
		cmd.add( FrameF );
		ValueArg<int> FrameI("","initial", "Initial frame index as defined in proj file (default= first)",false,0,"int");
		cmd.add( FrameI );

        ValueArg<std::string> Proj("p","proj", "ED-trajectory projection file for filtration. The projection file with the components for each frame should be provided here.",false,"proj-file","string");
        cmd.add( Proj );

		ValueArg<int> MovType("","mov", "Motion Type (default=2): "
				"0=K-matrix, 1=V/W-arrays, 2=Simple-Rotations, 3=Linear (if Cartesian modes).",false,2,"int");
		cmd.add( MovType );
		SwitchArg LinealMotion("","linear", "Enables linear motion instead of sinusoidal bounce (default=disabled).", true);
		cmd.add( LinealMotion );
        ValueArg<int> ModelOut("","model_out", "Output Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom. "
        		"The default output model will be that selected with the -m option.",false,-1,"int");
        cmd.add( ModelOut );
		SwitchArg Cart("","cart", "Mandatory if Cartesian modes are supplied, otherwise moving in Internal Coordinates (default=disabled).", true);
		cmd.add( Cart );

	    ValueArg<std::string> FixFile("f","fixFile", "Input ASCII file defining the ICs that were fixed during NMA (default=disabled). "
	    		"If modes were computed removing arbitrary ICs, the user must introduce here the file generated by iMode's --save_fixfile option.",false,"fixstring","string");
        cmd.add( FixFile );

//		ValueArg<float> MaxF("a","amp", "Amplitude factor applied (default=1)",false,1,"float");
//		cmd.add( MaxF );
		ValueArg<double> Temp("T","temperature", "Temperature [K] (default=300).",false,300,"double");
		cmd.add( Temp );
//		ValueArg<double> AmpF("a","amp", "Amplitude energy factor. Mode energy = a*Kb*T (default=5)",false,5,"float");
//		ValueArg<double> AmpF("a","amp", "Amplitude linear factor. Mode amplitude =  factor * sqrt( 5*Kb*T/lambda,i ) (default=1)",false,1,"float");
		ValueArg<double> AmpF("a","amp", "Amplitude linear factor to scale motion (default=2).",false,2,"float");
		cmd.add( AmpF );
		ValueArg<int> Fram("c","frames", "Number of conformations generated (default=11). It should be an odd number!",false,11,"int");
		cmd.add( Fram );
//		ValueArg<int> Nev("n","nev", "Mode number excited (1,2,...,size) (default=1)",false,1,"int");
//		cmd.add( Nev );
		SwitchArg Chi("x","chi", "Considers first CHI dihedral angle (default=disabled).", true);
		cmd.add( Chi );
//		ValueArg<int> Model("m","model", "Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom, 3=NCAC, 4=CA-only. Note input model should match input data! (default=2)",false,2,"int");
		ValueArg<int> Model("m","model", "Input modes Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom, 3=NCAC, 4=CA-only (default=2).",false,2,"int");
        cmd.add( Model );

		// Defining required arguments (not labeled)
		UnlabeledValueArg<std::string> Pdb("pdb1","PDB input file.","default","in_pdb");
		cmd.add( Pdb );
		UnlabeledValueArg<std::string> Ptraj("ptraj1","Normal modes input file (.evec).","default","ptraj");
		cmd.add( Ptraj );
		UnlabeledValueArg<std::string> Out("out1","Output Multi-PDB file.","default","out_pdb");
		cmd.add( Out );
		UnlabeledValueArg<int> Nev("nev", "Mode number to be moved (1,2,...,size). If --proj option enabled, then it sets the number of components to be used in ED-trajectory filtration.",1,"int");
		cmd.add( Nev );

		// Parsing the command line.
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		strcpy(file_pdb,((temp=Pdb.getValue()).c_str())); // Gets PDB file name
		strcpy(file_ptraj,((temp=Ptraj.getValue()).c_str())); // Gets PDB file name
		strcpy(file_out,((temp=Out.getValue()).c_str())); // Gets PDB file name
//		amp_factor = AmpF.getValue();
		factor = AmpF.getValue();
        Ta = Temp.getValue();
		nev = Nev.getValue();
//		verb = Verb.isSet();
		ndiv = Fram.getValue();

		if(NoModel.isSet())
	    	nomodel = true;

	    if(Model.isSet())
		{
			// "normal" run (Input model/chi equal to output)
			if( !ModelOut.isSet() && !ChiOut.isSet() ) // (no ModelOut and no ChiOut)
			{
				model = Model.getValue();
				modeli = -1; // forces "normal" run
				if(Chi.isSet())
					type = 2;
				else
					type = 0;
			}
			// if input model/chi different from output
			if( ModelOut.isSet() || ChiOut.isSet() )
			{
				modeli = Model.getValue();
				if(Chi.isSet())
					typei = 2;
				else
					typei = 0;
				model = ModelOut.getValue();
				if(ChiOut.isSet())
					type = 2;
				else
					type = 0;
			}
			if(model == 4) // if CA-only, then cartesian motion...
			{
				nomodel = true; // disabling formatting (reduces dirty warning messages when CA-only model...)
				mov_type = 3;

			}
		}
		else
			model = 2;

		if(Chi.isSet())
			type = 2;
		else
			type = 0;

		if(MovType.isSet())
			mov_type = MovType.getValue();
		else
			if(Cart.isSet())
				mov_type = 3; // lineal motion
			else
				mov_type = 2; // simple rotations motion

		lineal_switch = LinealMotion.isSet();

        if(Proj.isSet())
    	{
    		strcpy(proj_file,((temp=Proj.getValue()).c_str())); // Gets ED-Trajectory-Projection-file
    		proj_switch = true;
    	}

		if(FrameI.isSet())
		{
			initial = FrameI.getValue();
			printf("pcatool> Parser: Initial frame= %d\n",initial);
			// initial--;
		}
		if(FrameF.isSet())
		{
			final = FrameF.getValue();
			printf("pcatool> Parser: Final frame= %d\n",final);
			// final--;
		}

        if(FixIC.isSet())
    	{
    		strcpy(fix_file,((temp=FixIC.getValue()).c_str())); // Gets Fix-file
    		fixmodel = 2;
    	}

        if(HarmonicEnergy.isSet())
    	{
    		strcpy(k_file,((temp=HarmonicEnergy.getValue()).c_str())); // Gets Fix-file
    		harmonic_energy = true;
    	}

        if(FixFile.isSet())
    	{
    		strcpy(fix_file,((temp=FixFile.getValue()).c_str())); // Gets Fix-file
    		fixmodel = 3;
    	}


		if(ndiv%2 == 0) // checking input
		{
			printf("Parser> Sorry, \"--frames\" should be an odd number!\n");
			exit(1);
		}

		if(NThreads.isSet())
		{
			nthreads = NThreads.getValue();
			printf("Parser> Using %d threads.\n",nthreads);
		}

	}
	catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }

}

// Computes the harmonic energy of "mol" by using the provided force constants matrix "decint"
double molenergy(Macromolecule *mol,float *coord,twid *decint,int nipa)
{
	float *coord2;
	double energy=0.0;
	int k,l;
	mol->coordMatrix(&coord2); // get single row coordinates from current Macromolecule

	for(int i=0; i<nipa; i++)
	{
		k = 3*decint[i].k;
		l = 3*decint[i].l;
		energy += decint[i].C * ( pow(coord[k]-coord2[k],2) + pow(coord[k+1]-coord2[k+1],2) + pow(coord[k+2]-coord2[k+2],2) );
	}

	free(coord2);
	return(energy);
}
