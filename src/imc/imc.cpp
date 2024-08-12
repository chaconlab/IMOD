/*************************************************************************
 *                           iMC                                         *
 *************************************************************************
 * This program is part of iMOD: http://chaconlab.org/imod/index.html    *
 * (c) Jose Ramon Lopez-Blanco, Jose Ignacio Garzon and Pablo Chacon.    *
 * IQFR-CSIC's Structural Bioinformatics Group (2004-2011).              *
 *************************************************************************
 *                                                                       *
 *   It applies either Cartesian or Internal Coordinates modes to a PDB  *
 *   model to generate a Monte-Carlo trajectory (Harmonic Energy).       *
 *   Motion can be carried out different ways:                           *
 *     -lineally (CCS Modes),                                            *
 *     -K-matrix (= Jacobian, ICS modes)                                 *
 *     -V/W-arrays (= Jacobian, ICS modes, memory efficient              *
 *     -Simple Rotations Scheme (ICS modes, very fast).                  *
 *   (It takes into account Multiple-Chains and different CG-models)      *
 *                                                                       *
 *************************************************************************
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <stdio.h> // needed by some compilers
#include "cmdl/CmdLine.h"
//#include "libnma/nma.h" // Mon's NMA library (Dihedral included)
#include "libnma/include/libnma_misc.h" // Mon's NMA Internal Coordinates related library
#include "libnma/include/libnma_io.h" // Mon's NMA Input-Output library
#include "libnma/include/libnma_cg.h" // Mon's NMA Coarse-Graining library
#include "libnma/include/libnma_move.h" // Mon's NMA IC-Motion library
#include "libnma/include/libnma_deriv.h" // Mon's NMA Derivatives library

// double KbT = 0.00198717 * 300; // Kboltz * Tempereature(K)
double Kb = 0.00198717; // Kboltz (in Angstroms)
// Boltzmann constant = 1.3806503 � 10-23 m2 kg s-2 K-1
// SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM   1.9872065e-3L
// 	This is the gas constant R in (kcal/mol)/K.
// SimTK_BOLTZMANN_CONSTANT_KCAL_ANGSTROM   SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM
// 	Boltzmann's constant in Kcal-Angstrom units of (kcal/mol)/kelvin; same as R.

// double rfactor= 2.75 * sqrt((double)2); // "agressivity"-related factor (amplitude-related)
//                                // sqrt(2) ==> because the eigenvalues are /2 in "m-mc-eigen_mon.pl"
double rfactor = 3.889087297*2; // 7.778174594

char version[]="1.11"; // version code
char file_pdb[FILE_NAME];
char file_ptraj[FILE_NAME];
char fix_file[FILE_NAME];
char file_out[FILE_NAME];
char name[FILE_NAME];
char dummy_string[FILE_NAME];
int nmodes=0; // number of modes to be considered
int nframes=0; // number of frames
int nmovs=0; // number of iterations (movements) per frame
int traj_format=2; // Trajectory output format
int mov_method=0; // Movement type
int type=0; // Dihedral-model type
int model=0; // Coarse Graining model
int ndiv=5; // number of iterations per K-matrix model generated
int fixmodel = 0; // =0 --> no ICs fixation
int valid=0; // valid biased structures
int novalid=0; // not valid biased structures
int maxframes=100000; // maximum number of trial models
int maxtest=20000; // maximum number of test models
int wintest=50; // window size for testing stage
unsigned int seed;
double KbT;  // Kboltz * Tempereature
double Ta; // Temperature
double Ef; // Energy factor
double linear_factor = 1.0; // linear amplitude factor
float tscore = 0.0; // Target score
float cscore = 0.0; // Current score
float scorethr = 1.0; // Default score tolerance (A)
float scorerat = 0.01; // Default radius of gyration tolerance ratio
float score = 0.0; // Current score
float score0 = 0.0; // Initial score
float optfactor = 0.05; // Optimization stage factor
float optscale = 1.0; // Optimization final scale factor
float optrmsd = 0.0; // Target RMSD for optimization
float nevec_fact = -1; // % of eigenvectors to be computed (set by parser)
bool weight_switch=false; // takes into account mass weighting
bool verb_switch=false; // verbose output
bool var_switch=false; // =false --> force / =true --> variance
bool firstwrite_switch=false; // = true --> writes the initial model as first frame
bool fixFrag_switch = false;
bool fixIC_switch = false;
bool read_fixDH_switch = false;
bool saveformat_switch = false; // = true --> save formated input PDB
bool savemodel_switch = false; // = true --> save model PDB
bool nomodel = false; // true = CG-model building and formating disabled (initial model)
bool validconf = true; // Sometimes false for Rg biased conformational sampling.
bool filter_switch = false; // True enable biased conformational sampling.
bool opt_switch = false; // True enable energy/stiffness optimization
bool rg_switch = false; //
bool rmsd_switch = false; //
//bool carmsd_switch = true; // During Ef optimization and Filtration stages, CA-RMSD wil be used

extern CRandomMersenne *rg; // Mersenne Twister global object

using namespace TCLAP;
void parseOptions(int argc, char** argv);

int main( int argc, char * argv[] )
{
	printf("imc>\nimc> Monte-Carlo IC-NMA based program v%s\nimc>\n",version);

	// Parsing Input
	parseOptions(argc,argv);

	// CA  Conditions
	Condition *calpha = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	calpha->add(" CA "); // For Proteins
	calpha->add(" P  "); // For Nucleic acids
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

	// Mersenne Twister Seed Initialization
	// (It outputs a random double number in the interval 0 <= x < 1, with rg->Random() )
	rg = new CRandomMersenne( seed );

	Htimer ht_timer,ht_timer2; // timers
	KbT = Kb * Ta; // setting thermal-energy
	int natoms;

	printf( "imc> Reading Input PDB file\n" );
	Macromolecule *molr, *mol, *molNCAC;
	molr = new Macromolecule( file_pdb );
	molr->readPDB(file_pdb);
	molr->info(stdout);
//	natoms = molr->get_num_atoms();

	// FORMATS THE PDB FIRST !!!
	if(!nomodel)
	{
		printf( "imc> Formatting pdb names\n" );
		molr->format_residues(); // Forces predefined atom sorting
		if(saveformat_switch)
		{
			sprintf(dummy_string,"%s_format.pdb",name);
			molr->writePDB( dummy_string );
		}
	}

	// Setting Coarse-Graining model
	switch(model)
	{
	case 0: // CA-only + (NH and CO)-terminal model
		printf( "imc> Coarse-Graining model: CA-model\n");

		// N,CA,C selection
		molNCAC = molr->select( ncac2 );

		// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
		// Warning, atoms are not copied, they're just pointers to the original atoms.
		// setmass = true --> adds masses to Occupancy and Bfactor, otherwise left unchanged.
		// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
		//                   (unit mass is divided between NH or CO and their CAs)
		// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
		//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
		if(!nomodel) // Makes model (if it's needed)
			mol = cg_CA( molNCAC, true, false ); // masses will be computed
		else
			mol = cg_CA( molNCAC, false, false ); // input pdb masses will be used

		// Saving "NCAC" model
		if(savemodel_switch)
		{
			sprintf(dummy_string,"%s_ncac.pdb",name);
			molNCAC->writePDB( dummy_string );
		}
		break;

	case 4: // CA-only
		printf( "imc> Coarse-Graining model: CA-only\n");
		// if CA selection
		mol = molr->select( calpha2 );
		break;

	case 3: // N,CA,C-model
		printf( "imc> Coarse-Graining model: N,CA,C-model\n");

		// N,CA,C selection
		mol = molr->select( ncac2 );
		if(!nomodel) // Set masses (if it's needed)
			mass_NCAC( mol, false );
		break;

	case 1: // 3BB2R model
		printf( "imc> Coarse-Graining model: 3BB2R\n");
		mol = molr;
		if(!nomodel) // Makes 3BB2R model (if it's needed)
		{
			// CREATES a 3BB2R reduced model
			//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
			//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
			//     There are a few exceptions: for Ala, Gly and Pro,
			//     and for the 1st and last residues.
			printf("imc> Creating 3BB2R reduced model:\n");
			cg_3BBR2(mol);
		}
		break;

	case 2: // Full-Atom
		printf( "imc> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
		mol = molr;
		if(!nomodel) // Sets masses (if it's needed)
			mass_FA(mol, false); // sets masses
		break;
	}
	natoms = mol->get_num_atoms();
	printf( "imc> Selected model number of (pseudo)atoms: %d\n", natoms );

	// Output file name (different extensions according to output file type)
	FILE *f_out;
	switch(traj_format)
	{
	case 0:
		sprintf(file_out,"%s.mc",name);
		break;
	case 1:
		sprintf(file_out,"%s.pdb",name);
		break;
	case 2:
		sprintf(file_out,"%s.x",name);
		break;
	}
	printf("imc> Opening output trajectory file: %s\n",file_out);
	f_out = fopen(file_out,"w"); // opening output file

	tri *props;
	int size=0;
	bool *fixed = NULL;
	int old_size;
	pdbIter *iter_atom = new pdbIter( molr );
	bool *addrotx=NULL; // Should be added ROTATIONs? (with fixing)

	if(mov_method != 3) // internal coordiantes normal modes
	{
		// Initializing some internal coords. related residue properties
		// (# atoms, # units, # internal coordinates, etc...)
		int *unat;

		switch(model)
		{
		case 0:
			properCA(molNCAC,&props,&unat);
			break;
		case 3:
			properCA(mol,&props,&unat);
			break;
		case 1:
		case 2:
			properMFA(mol,&props,&unat,type,model);
			break;
		}

		// Theoric number of DoFs= 3*T + 3*R -6 + Dihedral
		// (Before fixing, i.e. the fix-file format DoFs...)
		int old_size,seg_atoms=0;
		pdbIter *iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
		size = 0;
		old_size = size; // temp (Dihedral ICs)
		for( iter_seg->pos_fragment = 0; !iter_seg->gend_fragment(); iter_seg->next_fragment() ) // screen residues
			size += props[iter_seg->pos_fragment].nan; // Each residue may have different number of dihedral-angles!
		printf( "imc> Number of Dihedral angles: %d\n", size );
		old_size = size;
		for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // screen segments
		{
			seg_atoms = (iter_seg->get_segment())->num_atoms();
			if(seg_atoms > 1)
				size += 6;
			else
				size += 3;
		}
		printf( "imc> Rotational/Translational ICs (Non-Eckart): %d\n", size-old_size );
		size -= 6; // Eckart conditions
		printf( "imc> Predicted number of ICs (Eckart): %d\n", size );

		// Mobile Internal Coordinates selection
		if(fixmodel != 0)
		{
			if( !(fixed = (bool *)malloc(sizeof(bool)*size) ) ) // Allocate memory
			{
				printf("Unable to allocate memory. Forcing exit!\n");
				exit(1);
			}
			old_size = size;
			switch(fixmodel)
			{
			case 2:
				size = read_fixIC(fix_file,mol,props,fixed);
				break;
			case 3:
				size = read_fixDH(fix_file,mol,props,fixed,type,model);
				break;
			default:
				printf("Sorry, unknowk fixation method!\nForcing exit\n");
				exit(1);
				break;
			}
			printf("imc> Input CG-model Fixed Internal Coordinates: %d\n", old_size-size );
			printf("imc> Input CG-model Final Internal Coordinates (sizei) = %d\n",size);
		}

		// Creates two auxiliar arrays with segment properties (needed due to fixing):
		//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
		//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
		// (Allocates memory itself, if it's needed)
		int *effsegx=NULL; // Effective segment indices
		size = seg_props(mol, props, fixed, model, type, &addrotx, &effsegx);
		printf( "imc> Number of ICs predicted by seg_props(): %d\n", size );
		free(effsegx);
	}
	else // cartesian normal modes
		size = 3*natoms;

	// Reads ptraj file (allocating memory)
	double *invec,*inval,*evec_CA;
	int natoms_ptraj, nvects, ncomps;
	read_ptraj(file_ptraj, &invec, &inval, &natoms_ptraj, &nvects, &ncomps);
	// Showing some "ptraj" info...
//	printf("File: %s\n\tatoms: %5d\n\tvectors: %5d\n\tcomponents: %5d\n",file_ptraj,natoms_ptraj,nvects,ncomps);
	printf("imc> Ptraj info: %s  vectors=%d  components=%d\n",file_ptraj,nvects,ncomps);

	// Some checking: Number of ptraj vector components must match internal coordinates found in PDB (size)
	if( ncomps != size )
	{
		printf("imc> Sorry!\n"
				"imc> Mismatch between ptraj file (%d) and coarse-graining model (%d).\n"
				"imc> Please, check input: pdb, ptraj or --type\n",ncomps,size);
		exit(1);
	}

	// Number of eigenvectors to be computed (we need to know "size" first)
    if(nevec_fact >= 1.0) // number of modes
    	nmodes = (int) nevec_fact;
    else if(nevec_fact > 0)
    	nmodes = (int) (nevec_fact * size);
    else // if nevec_fact < 0 (default)
    	nmodes = nvects; // setting to the maximum available
	if(nvects < nmodes)
	{
		nmodes = nvects; // setting to the maximum available
		printf( "imc> Warning: More modes requested than available; setting \"nmodes\" to the maximum (%d)\n",nmodes);
	}
	printf( "imc> Number of used modes: %d\n", nmodes);

	// Un-mass-weighting eigenvectors (if input vectors are mass-weighed)
	if(weight_switch)
	{
		float *rmasses;
		rmasses = (float *) malloc( sizeof(float) * 3 * natoms );

		// Getting reduced-masses
		for(iter_atom->pos_atom=0; !iter_atom->gend_atom(); iter_atom->next_atom() )
		{
			rmasses[3*iter_atom->pos_atom] =
				rmasses[3*iter_atom->pos_atom+1] =
					rmasses[3*iter_atom->pos_atom+2] = sqrt( ( iter_atom->get_atom() )->getPdbocc() );
		}

		// Un-mass-weighting eigenvectors
		for(int i=0; i< nmodes; i++)
		{
			for(int j=0; j< ncomps; j++)
				invec[ncomps*i + j] /= rmasses[j];
		}
	}

	// Converting eigenvalues to Force Constants (if they are variances from PCA)
	// see: ref. from Karplus & McCammon "Quasi-harmonic analisys..." elsewere...
	if(var_switch)
		for(int i=0; i< nmodes; i++)
			inval[i] = KbT/inval[i];

	// Multi-PDB trajectory (First frame)
	if(traj_format == 1)
	{
		fclose(f_out); // deleting previous file...
		if(firstwrite_switch)
			molr->writeMPDB(file_out,1); // writing first frame (initial model)
			// In VMD, if the first frame is distorted, then the
			// wire representation looks very ugly.
	}

	// Center of Mass positioning
	// K-matrix computation needs the PDB model placed in its CoM
	// (each snapshot will be placed back to its initial position below)
	double rd[3];
	if(mov_method==0 || mov_method==1)
	{
		// Shifting to the CoM
		// Computing the PDB's Center of Mass (CoM)
		double mtot = 0.0;
		double mta = 0.0;
		rd[0] = rd[1] = rd[2] = 0.0;
		Atom *atom0;
		Tcoor pos0;
		for ( iter_atom->pos_atom = 0; !iter_atom->gend_atom(); iter_atom->next_atom() )  // screens all-atoms
		{
			atom0 = iter_atom->get_atom();
			atom0->getPosition(pos0);
			mta = atom0->getPdbocc(); // Load mass...
			mtot += mta;
			// Sum(mass*coord) before putting the CoM at 0
			rd[0] += mta * pos0[0];
			rd[1] += mta * pos0[1];
			rd[2] += mta * pos0[2];
		}
		rd[0] /= mtot;
		rd[1] /= mtot;
		rd[2] /= mtot;

		// shift CM of pdb_model to 0,0,0
		for ( iter_atom->pos_atom = 0; !iter_atom->gend_atom(); iter_atom->next_atom() )  // screens all-atoms
		{
			atom0 = iter_atom->get_atom();
			atom0->getPosition(pos0);
			pos0[0] -= rd[0];
			pos0[1] -= rd[1];
			pos0[2] -= rd[2];
			atom0->setPosition(pos0);
		}
	}

	// MAIN LOOP's variables
	int niters; // maximum number of MC iterations
	int mode; // selected mode
	int iter; // current iteration index
	int accepted = 0; // number of accepted movements
	int cont,index,k,l,i_atom;
	double ds,dx,delta,increment;
	double enertot = 0; // total energy
	double enertry = 0;
	double dpos[3];
	double deriv[3];
	Tcoor pos;
	Macromolecule *molx,*molxCA,*molCA; // allocating copy of the readed macromolecule
	pdbIter *iter_mol,*iter_molCA,*iter_temp;
	Atom *atom;
	float *coord;
	trd *der = NULL; // K-matrix pointer (NULL forces memory allocation)
	M4Rot *matrix4_op;
	double ***V = NULL, ***W = NULL;
//	bool **body1 = NULL;
	int **body1 = NULL; // (num_atoms x 4) sized integer matrix to store body 1/2 boundary information

	// Allocating "moving" Macromolecule only once!
	if(model == 0 ) // CA-model
	{
		if(mov_method == 3) // CA-model && lineal motion
			molx = new Macromolecule( mol ); // creating copy of the model molecule
		else
			molx = new Macromolecule( molNCAC ); // creating copy of the model molecule
	}
	else // remaining models: NCAC, 3BB2R, Full-Atom, and CA-only
		molx = new Macromolecule( mol ); // creating copy of the model molecule
	iter_mol = new pdbIter( molx ); // Just creating the iterator once!

	if(model == 0 && mov_method != 3) // CA-model & Not lineal motion
	{
		// Creates (allocate memory) a CA-model with first NH and last CO pseudo-atoms of each segment.
		// Warning, atoms are not copied, they're just pointers to the original atoms.
		molxCA = cg_CA( molx, false, false ); // already set masses will be used
		molCA = cg_CA( molNCAC, false, false ); // already set masses will be used
		if(traj_format == 2)
			iter_molCA = new pdbIter( molxCA );
	}

	if(filter_switch)
	{
		if(rg_switch)
		{
			if(model == 0 && mov_method != 3) // CA-model & Not lineal motion
				score = radius_gyration(molxCA);
			else
				score = radius_gyration(molx);
			printf( "imc> Initial model radius of gyration (A): %f\n", score);
			printf( "imc> Target Rg (A): %f +- %f\n",tscore,scorethr);
		}
		else if(rmsd_switch)
		{
			score = 0.0; // RMSD to inital model
			printf( "imc> Initial model RMSD (A): %f\n", score);
			printf( "imc> Target Rg (A): %f +- %f\n",tscore,scorethr);
		}
		else
		{
			printf("imc> Sorry, you must choose a valid filtering method\n");
			exit(2);
		}
	}

	// Displacement vector allocation ("uu")
	double *uu; // current motion
	if(traj_format != 0)
	{
		uu = (double *) malloc( sizeof(double) * size );
		if( !(uu = (double *) malloc( sizeof(double) * size )) )
		{
			printf("imc> Memory allocation error (uu-vector)!\n");
			exit(1);
		}
	}

	// Allocating scaling factor array (normal-mode space)
	double *scaling = (double *) malloc( sizeof(double) * nmodes );
	// Allocating absolute conformation (normal-mode space)
	double *xconf = (double *) malloc( sizeof(double) * nmodes );
	// Allocating Energy (normal-mode space)
	double *ener = (double *) malloc( sizeof(double) * nmodes );
	// Initialization
	for(int i=0; i< nmodes; i++)
	{
		xconf[i] = 0.0;
		ener[i] = 0.0;
	}
	// Computing scaling factors to have an "isotropic step"
	// 		see: Noguti & Go. Biopolymers, Vol 24, 527-46 (1985).
	// also see Manu's script: "m-mc-eigen_mon.pl"
	// Note that "sqrt( KbT/inval[i] )" is proportional to the RMSD due to mode "i";
	// 		see: Go, N. "A theorem on amplitudes of thermal atomic fluctuations in large molecules
	// 		assuming specific conformations calculated by NMA". Biophys. Chem., 53 (1990) 105-112
	for(int i=0; i< nmodes; i++)
	{
		inval[i] *= Ef; // Applying first the user-introduced Energy-factor (another scaling term)
		scaling[i] = rfactor * sqrt( KbT/inval[i] ); // Computing scale factor
	}
	//	scaling[i] = rfactor * ( KbT/inval[i] );
	// "rfactor" is related to the step amplitude, and is included into the scaling factor

	// STIFFNESS TEST: It determines optimal E-factor to obtain the desired average RMSD.
	int iavg = 0;
	double mavg=0.0; // current moving average
	valid = 0;
	if(opt_switch)
	{
		Ef = 1.0; // reset the energy factor
		printf("imc> Energy (stiffness) optimization stage. Target RMSD= %f\n",optrmsd);

		// Allocating Energy (normal-mode space)
		float *scores = (float *) malloc( sizeof(float) * wintest );

		for(iter=1; valid < maxtest; iter++)
		{
			// Choosing random mode
			mode = (int) ( nmodes * rg->Random() ); // [0:1) Playing dice with Mersenne!
			// Choosing random displacement
			dx = ( (double) rg->Random() ) - 0.5; // random displacement (-0.5:0.5) (in normal-mode space)
			// Scaling displacement increment (anisotropic step)
			ds = scaling[mode] * dx; // "ds" is the scaled amplitude of the motion (in normal-mode space)
			// Energy increment due to "ds" in "mode"
			enertry = 0.5 * inval[mode] * Ef * pow(xconf[mode] + ds,2) - ener[mode];

			// Metropolis et al. acceptance test
			if( enertry < 0 || exp(-enertry/KbT) > (double) rg->Random() )
			{
//				accepted++;
				enertot += enertry; // updating total energy
				ener[mode] += enertry; // updating mode energy
				xconf[mode] += ds; // updating conformation (in normal-mode space)
			}

			// Saving (and generating) frame each "nmovs" iterations
			if(iter % nmovs == 0)
			{
				valid++; // all are valid ones

				// Cartesian Coordinate Space (CCS) trajectory
				// Always moving from initial structure (minimize errors)
				// Copy trial macromolecule coordinates (faster than allocation...)
				if(model == 0) // CA-model
				{
					if(mov_method == 3) // CA-model && lineal motion
						molx->copy_coordinates(mol); // copy initial coords
					else
						molx->copy_coordinates(molNCAC); // copy initial NCAC coords.
				}
				else // remaining models: NCAC, 3BB2R, Full-Atom, and CA-only
					molx->copy_coordinates(mol); // copy initial coords.

				for(int l=0; l< size; l++) // screen Dihedrals
					uu[l] = 0.0; // initialization
				for(int i=0; i< nmodes; i++) // screen used modes
					for(int l=0,index=ncomps*i; l< size; l++,index++) // screen Dihedrals
						uu[l] += invec[index] * xconf[i];
				for(int l=0; l< size; l++) // screen Dihedrals
					uu[l] *= linear_factor;

				printf("imc> Testing frame %5d  Ef= %10.4e ",valid,Ef);

				// Translating from Normal Mode space into Cartesian space
				// (updating atom positions)
				switch(mov_method)
				{
				case 3: // LINEAL MOTION
					// (from cartesian Normal Modes)
					// It applies many consecutive CCS Normal Modes with their amplitudes
					move_cart(molx, uu, 1.0);
					break;

				case 0: // MOVING LINEARLY, BUT UPDATING K-MATRIX EACH TIME
					// (ICS --> CCS, via Jacobian)

					// Divide motion (uu-vector)
					for(int l=0; l< size; l++) // screen Dihedrals
						uu[l] /= ndiv; // K-matrix motion increment

					for(int k = 1; k <= ndiv; k++) // ndiv=10 seems to be a good aprox.
					{
						// Getting coordinates single row (pseudo-atom model)
						molx->coordMatrix( &coord );

						// Compute derivatives of cartesian coordinates respect dihedrals
						// (Warning, "mol" should be placed on the Center of Mass)
						switch(model)
						{
						case 0:
						case 3:
							dydqMCA3x(molx,coord, &der, size, fixed);
							break;
						case 1:
						case 2:
							dydqMFAx(molx,coord, props, &der, type, model, size, fixed, addrotx);
							break;
						}

						// MOVING
						for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
						{   // screen atoms
							i_atom = size * iter_mol->pos_atom; // a bit faster...
							atom = iter_mol->get_atom();
							atom->getPosition(pos); // get current position
							dpos[0] = dpos[1] = dpos[2] = 0.0;

							for(int l=0; l< size; l++) // screen Dihedrals
							{
								dpos[0] += der[i_atom + l].x * uu[l];
								dpos[1] += der[i_atom + l].y * uu[l];
								dpos[2] += der[i_atom + l].z * uu[l];
							}

							pos[0] += (float) dpos[0];
							pos[1] += (float) dpos[1];
							pos[2] += (float) dpos[2];
							atom->setPosition(pos); // set new position
						}
						free(coord);
					}

					// Placing frame back to the initial position
					for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
					{   // screen atoms
						atom = iter_mol->get_atom();
						atom->getPosition(pos); // get current position
						pos[0] += (float) rd[0]; // rd[] places it back to initial origin
						pos[1] += (float) rd[1];
						pos[2] += (float) rd[2];
						atom->setPosition(pos); // set new position
					}
					break;

				case 1: // MOVING LINEARLY WITH V/W-ARRAYS, BUT UPDATING EACH TIME
					// (ICS --> CCS, via Jacobian)

					// Divide motion (uu-vector)
					for(int l=0; l< size; l++) // screen Dihedrals
						uu[l] /= ndiv; // K-matrix motion increment

					for(int k = 1; k <= ndiv; k++) // ndiv=10 seems to be a good aprox.
					{
						// Getting coordinates single row (pseudo-atom model)
						molx->coordMatrix( &coord );

						// Compute derivatives of cartesian coordinates respect dihedrals
						// (Warning, "mol" should be placed on the Center of Mass)
						switch(model)
						{
						case 0:
						case 3:
							vwMCA3x(molx,coord,&V,&W,&body1,size,fixed);
							break;
						case 1:
						case 2:
							vwMFAx(molx,coord,props,&V,&W,&body1,type,model,size,fixed,addrotx);
							break;
						}

						//
						// MON: FIX THIS, USE V/W MOTION ROUTINE (SINGLE AND MULTI THREADED) !!!!
						//
						// MOVING
						for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
						{   // screen atoms
							atom = iter_mol->get_atom();
							atom->getPosition(pos); // get current position
							dpos[0] = dpos[1] = dpos[2] = 0.0;

							for(int l=0; l< size; l++) // screen Dihedrals
							{
								if( body1[l][iter_mol->pos_atom] ) // if body 1 atom
								{
									// i-der for k-atom --> = v + (w x r) (vectorial product)
									deriv[0] = V[l][0][0] + W[l][0][1] * pos[2] - W[l][0][2] * pos[1];
									deriv[1] = V[l][0][1] + W[l][0][2] * pos[0] - W[l][0][0] * pos[2];
									deriv[2] = V[l][0][2] + W[l][0][0] * pos[1] - W[l][0][1] * pos[0];
								}
								else // if body 2 atom
								{
									// i-der for k-atom  --> = v + (w x r) (vectorial product)
									deriv[0] = V[l][1][0] + W[l][1][1] * pos[2] - W[l][1][2] * pos[1];
									deriv[1] = V[l][1][1] + W[l][1][2] * pos[0] - W[l][1][0] * pos[2];
									deriv[2] = V[l][1][2] + W[l][1][0] * pos[1] - W[l][1][1] * pos[0];
								}

								dpos[0] += deriv[0] * uu[l];
								dpos[1] += deriv[1] * uu[l];
								dpos[2] += deriv[2] * uu[l];
							}

							pos[0] += (float) dpos[0];
							pos[1] += (float) dpos[1];
							pos[2] += (float) dpos[2];
							atom->setPosition(pos); // set new position
						}
						free(coord);
					}

					// Placing frame back to the initial position
					for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
					{   // screen atoms
						atom = iter_mol->get_atom();
						atom->getPosition(pos); // get current position
						pos[0] += (float) rd[0]; // rd[] places it back to initial origin
						pos[1] += (float) rd[1];
						pos[2] += (float) rd[2];
						atom->setPosition(pos); // set new position
					}
					break;

				case 2:	// MOVING IN INTERNAL COORDINATES
					// (from internal coordinates normal modes)
					float matrix4[4][4];

					if(model == 3 || model == 0) // NCAC- or CA-model
					{
						move_dihedralMCAx(molx,uu, props,1.0,type,model,fixed);
						molNCAC->minRmsd(molx,matrix4);
					}
					else
					{
						move_dihedralMFAx(molx,uu, props,1.0,type,model,fixed,addrotx);
						mol->minRmsd(molx,matrix4);
					}

					// Alignment respect to the first PDB model
					// The protein is re-built from N-terminal every time we move it,
					// so we need to align it respecting the initial PDB model. (It's fast and works nice)
					matrix4_op = new M4Rot(matrix4);
					molx->applyAtoms(matrix4_op);
					delete(matrix4_op);
					break;
				}

				// Computing Moving average
				if(iavg >= wintest) // reset index if outside window...
					iavg = 0;
				// Score computation (RMSD to inital model)
				if(model == 0 && mov_method != 3) // CA-model & Not lineal motion
					score = molCA->rmsd(molxCA);
				else
					score = mol->rmsd(molx);
				scores[iavg] = score;
				iavg++; // next window element

				if(valid >= wintest) // If there are enough score values to compute window
				{
					mavg = 0.0;
					for(index=0; index<wintest; index++)
						mavg += scores[index];
					mavg /= wintest;

					// Updating Energy/Stiffness factor ("Ef")
					if(valid % wintest == 0 && valid != maxtest)
					{
						cscore = mavg-optrmsd*optscale;
						Ef *= 1.0 + optfactor * cscore; // it should be lower/higher depending on "cscore" sign
						// Updating the scaling factors for next iters
						for(int i=0; i< nmodes; i++)
							scaling[i] = rfactor * sqrt( KbT/(inval[i]*Ef) ); // Computing scale factor
					}
				}

				printf("%f %f %s\r",mavg,score,ht_timer.print_time_sec());
				fflush(stdout);
				ht_timer.restart(); // timer
			}
		}
		printf("\nimc> Scores:");
		for(index=0; index<wintest; index++)
		{
			printf(" %.3f",scores[index]);
		}
		free(scores);
		printf("\nimc>\n");

		// Apply the final "Ef"
		for(int i=0; i< nmodes; i++)
		{
			inval[i] *= Ef;
			scaling[i] = rfactor * sqrt( KbT/inval[i] ); // Computing scale factor
		}
	}

	// MONTE CARLO MAIN LOOP
	ht_timer.restart(); // timer
	ht_timer2.restart(); // global timer
	valid = 0;
	printf("imc> Monte-Carlo stage using Ef= %f\n",Ef);
	for(iter=1; valid < nframes && novalid < maxframes; iter++)
	{
		// Choosing random mode
		mode = (int) ( nmodes * rg->Random() ); // [0:1) Playing dice with Mersenne!
		// Choosing random displacement
		dx = ( (double) rg->Random() ) - 0.5; // random displacement (-0.5:0.5) (in normal-mode space)
		// Scaling displacement increment (anisotropic step)
		ds = scaling[mode] * dx; // "ds" is the scaled amplitude of the motion (in normal-mode space)
		// Energy increment due to "ds" in "mode"
		enertry = 0.5 * inval[mode] * pow(xconf[mode] + ds,2) - ener[mode];

		// Metropolis et al. acceptance test
		if( enertry < 0 || exp(-enertry/KbT) > (double) rg->Random() )
		{
			accepted++;
			enertot += enertry; // updating total energy
			ener[mode] += enertry; // updating mode energy
			xconf[mode] += ds; // updating conformation (in normal-mode space)
		}
//printf("mon> enertry= %f  enertot= %f\n",enertry, enertot);

		// Saving (and generating) frame each "nmovs" iterations
		if(iter % nmovs == 0)
		{
			if(traj_format == 0) // Normal-Mode Space (NMS) trajectory (just plain-text mode amplitudes)
			{
				// Output frame header
				fprintf(f_out, "****\n");
				fprintf(f_out, "ITERATION: %d, MODE: %d, ETOT:%12.4e, DX: %12.4e, CUM-DX(%d): %12.4e\n",iter,mode,enertot,ds,mode,0);
				fprintf(f_out, "DX for %d Modes\n",nmodes);
				fprintf(f_out, "-------------------------------------------------------------------\n");

				// Output Normal Mode coordinates
				cont=1;
				for(int j=0; j< nmodes; j++)
				{
					if ( cont % 10 != 0 )
					{
						fprintf(f_out, "%10.4f",xconf[j]);
						cont++;
					}
					else
					{
						fprintf(f_out, "%10.4f\n",xconf[j]);
						cont=1;
					}
				}
				if(cont==1)
					fprintf(f_out,"\n");
				else
					fprintf(f_out,"\n\n");
			}
			else // Cartesian Coordinate Space (CCS) trajectory
			{
				// Always moving from initial structure (minimize errors)
				// Copy trial macromolecule coordinates (faster than allocation...)
				if(model == 0) // CA-model
				{
					if(mov_method == 3) // CA-model && lineal motion
						molx->copy_coordinates(mol); // copy initial coords
					else
						molx->copy_coordinates(molNCAC); // copy initial NCAC coords.
				}
				else // remaining models: NCAC, 3BB2R, Full-Atom, and CA-only
					molx->copy_coordinates(mol); // copy initial coords.

				for(int l=0; l< size; l++) // screen Dihedrals
					uu[l] = 0.0; // initialization
				for(int i=0; i< nmodes; i++) // screen used modes
					for(int l=0,index=ncomps*i; l< size; l++,index++) // screen Dihedrals
						uu[l] += invec[index] * xconf[i];
				for(int l=0; l< size; l++) // screen Dihedrals
					uu[l] *= linear_factor;

				printf("imc> Frame %5d  Etot= %f ",iter/nmovs,enertot);

				// Translating from Normal Mode space into Cartesian space
				// (updating atom positions)
				switch(mov_method)
				{
				case 3: // LINEAL MOTION
						// (from cartesian Normal Modes)
					// It applies many consecutive CCS Normal Modes with their amplitudes
					move_cart(molx, uu, 1.0);
					break;

				case 0: // MOVING LINEARLY, BUT UPDATING K-MATRIX EACH TIME
						// (ICS --> CCS, via Jacobian)

					// Divide motion (uu-vector)
					for(int l=0; l< size; l++) // screen Dihedrals
						uu[l] /= ndiv; // K-matrix motion increment

					for(int k = 1; k <= ndiv; k++) // ndiv=10 seems to be a good aprox.
					{
						// Getting coordinates single row (pseudo-atom model)
						molx->coordMatrix( &coord );

						// Compute derivatives of cartesian coordinates respect dihedrals
						// (Warning, "mol" should be placed on the Center of Mass)
						switch(model)
						{
						case 0:
						case 3:
							dydqMCA3x(molx,coord, &der, size, fixed);
							break;
						case 1:
						case 2:
							dydqMFAx(molx,coord, props, &der, type, model, size, fixed, addrotx);
							break;
						}

						// MOVING
						for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
						{   // screen atoms
							i_atom = size * iter_mol->pos_atom; // a bit faster...
							atom = iter_mol->get_atom();
							atom->getPosition(pos); // get current position
							dpos[0] = dpos[1] = dpos[2] = 0.0;

							for(int l=0; l< size; l++) // screen Dihedrals
							{
								dpos[0] += der[i_atom + l].x * uu[l];
								dpos[1] += der[i_atom + l].y * uu[l];
								dpos[2] += der[i_atom + l].z * uu[l];
							}

							pos[0] += (float) dpos[0];
							pos[1] += (float) dpos[1];
							pos[2] += (float) dpos[2];
							atom->setPosition(pos); // set new position
						}
						free(coord);
					}

					// Placing frame back to the initial position
					for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
					{   // screen atoms
						atom = iter_mol->get_atom();
						atom->getPosition(pos); // get current position
						pos[0] += (float) rd[0]; // rd[] places it back to initial origin
						pos[1] += (float) rd[1];
						pos[2] += (float) rd[2];
						atom->setPosition(pos); // set new position
					}
					break;

				case 1: // MOVING LINEARLY WITH V/W-ARRAYS, BUT UPDATING EACH TIME
					    // (ICS --> CCS, via Jacobian)

					// Divide motion (uu-vector)
					for(int l=0; l< size; l++) // screen Dihedrals
						uu[l] /= ndiv; // K-matrix motion increment

					for(int k = 1; k <= ndiv; k++) // ndiv=10 seems to be a good aprox.
					{
						// Getting coordinates single row (pseudo-atom model)
						molx->coordMatrix( &coord );

						// Compute derivatives of cartesian coordinates respect dihedrals
						// (Warning, "mol" should be placed on the Center of Mass)
						switch(model)
						{
						case 0:
						case 3:
							vwMCA3x(molx,coord,&V,&W,&body1,size,fixed);
							break;
						case 1:
						case 2:
							vwMFAx(molx,coord,props,&V,&W,&body1,type,model,size,fixed,addrotx);
							break;
						}

						// MOVING
						for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
						{   // screen atoms
							atom = iter_mol->get_atom();
							atom->getPosition(pos); // get current position
							dpos[0] = dpos[1] = dpos[2] = 0.0;

							for(int l=0; l< size; l++) // screen Dihedrals
							{
								if( body1[l][iter_mol->pos_atom] ) // if body 1 atom
								{
									// i-der for k-atom --> = v + (w x r) (vectorial product)
									deriv[0] = V[l][0][0] + W[l][0][1] * pos[2] - W[l][0][2] * pos[1];
									deriv[1] = V[l][0][1] + W[l][0][2] * pos[0] - W[l][0][0] * pos[2];
									deriv[2] = V[l][0][2] + W[l][0][0] * pos[1] - W[l][0][1] * pos[0];
								}
								else // if body 2 atom
								{
									// i-der for k-atom  --> = v + (w x r) (vectorial product)
									deriv[0] = V[l][1][0] + W[l][1][1] * pos[2] - W[l][1][2] * pos[1];
									deriv[1] = V[l][1][1] + W[l][1][2] * pos[0] - W[l][1][0] * pos[2];
									deriv[2] = V[l][1][2] + W[l][1][0] * pos[1] - W[l][1][1] * pos[0];
								}

								dpos[0] += deriv[0] * uu[l];
								dpos[1] += deriv[1] * uu[l];
								dpos[2] += deriv[2] * uu[l];
							}

							pos[0] += (float) dpos[0];
							pos[1] += (float) dpos[1];
							pos[2] += (float) dpos[2];
							atom->setPosition(pos); // set new position
						}
						free(coord);
					}

					// Placing frame back to the initial position
					for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom())
					{   // screen atoms
						atom = iter_mol->get_atom();
						atom->getPosition(pos); // get current position
						pos[0] += (float) rd[0]; // rd[] places it back to initial origin
						pos[1] += (float) rd[1];
						pos[2] += (float) rd[2];
						atom->setPosition(pos); // set new position
					}
					break;

				case 2:	// MOVING IN INTERNAL COORDINATES
						// (from internal coordinates normal modes)
					float matrix4[4][4];

					if(model == 3 || model == 0) // NCAC- or CA-model
					{
						move_dihedralMCAx(molx,uu, props,1.0,type,model,fixed);
						molNCAC->minRmsd(molx,matrix4);
					}
					else
					{
						move_dihedralMFAx(molx,uu, props,1.0,type,model,fixed,addrotx);
						mol->minRmsd(molx,matrix4);
					}

					// Alignment respect to the first PDB model
					// The protein is re-built from N-terminal every time we move it,
					// so we need to align it respecting the initial PDB model. (It's fast and works nice)
					matrix4_op = new M4Rot(matrix4);
					molx->applyAtoms(matrix4_op);
					delete(matrix4_op);
					break;
				}

				if(filter_switch)
				{
					validconf = false;

					if(filter_switch)
					{
						if(rg_switch) // current Rg
						{
							if(model == 0 && mov_method != 3) // CA-model & Not lineal motion
								score = radius_gyration(molxCA);
							else
								score = radius_gyration(molx);
						}
						else if(rmsd_switch) // RMSD to inital model
						{
							if(model == 0 && mov_method != 3) // CA-model & Not lineal motion
								score = molCA->rmsd(molxCA);
							else
								score = mol->rmsd(molx);
						}
					}

					// Checking Score
//					if(fabsf(score-tscore) < scorethr) // if "in range"
					if(fabs(score-tscore) < scorethr) // if "in range"
						validconf = true; // enables model saving
					else
						novalid++;

					printf("%5d valid   %5d invalid  ",valid,novalid);
				}

				if(validconf) // If it is a valid conformation (if not biased, always true)
				{
					valid++;
					// Saving fame as Multi-PDB trajectory (.pdb)
					if(traj_format == 1)
						if(model == 0 && mov_method != 3) // CA-model & Not lineal motion
							molxCA->writeMPDB(file_out,(int)iter/nmovs);
						else
							molx->writeMPDB(file_out,(int)iter/nmovs);

					// Saving fame as AMBER's trajectory (.x)
					if(traj_format == 2)
					{
						// CA-model & Not lineal motion & Amber's traj_format
						if(model == 0 && mov_method != 3)
						{
							// Swaping iterators to save in AMBER's trajectory (.x) with CA-model
							iter_temp = iter_mol;
							iter_mol = iter_molCA;
						}

						cont = 1;
						for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom() )
						{
							( iter_mol->get_atom() )->getPosition(pos);
							for(int l=0; l< 3; l++) // x, y, z
							{
								// writing coord.
								if(cont % 10 != 0)
								{
									fprintf(f_out,"%8.3f",pos[l]);
									cont++;
								}
								else
								{
									fprintf(f_out,"%8.3f\n",pos[l]);
									cont = 1;
								}
							}
						}
						if(cont != 1)
							fprintf(f_out,"\n");

						iter_mol = iter_temp; // restoring swapped iterator after saving
					}
				}

				printf("%s\r",ht_timer.print_time_sec());
				fflush(stdout);
				ht_timer.restart(); // timer
			}
		}
	}

	if(mov_method == 0)
		free(der); // freeing K-matrix
	if(mov_method == 1)
	{
		// Lets Free V/W arrays (not used anymore)
		for(int i=0; i<size; i++)
		{
			free(V[i][0]);
			free(V[i][1]);
			free(V[i]);
			free(W[i][0]);
			free(W[i][1]);
			free(W[i]);
			free(body1[i]);
		}
		free(V);
		free(W);
		free(body1);
	}
	printf("\nimc>           Total        %s\n",ht_timer2.print_time_sec());

	// Computing aceptance ratio
	float aratio = (float) accepted/iter;
	printf( "imc>\nimc> ACEPTANCE RATIO(%%): %f\n",aratio*100);

	if(filter_switch)
		printf( "imc> Total frames: %d  ==>  Valid: %d %.2f%%   Invalid: %d %.2f%%\n",valid+novalid,valid,(float)(valid*100)/(valid+novalid),novalid,(float)(novalid*100)/(valid+novalid));

	if(traj_format == 0)
	{
		fprintf(f_out,"\n");
		fprintf(f_out, "-------------------\n");
		fprintf(f_out, "RATIO(%%): %f\n",aratio*100);
		fprintf(f_out, "-------------------\n");
	}

	if( traj_format != 1 )
		fclose(f_out); // closing output file

	printf("imc>\nimc> Success!\nimc> Written MC-trajectory: %s (%d frames)\nimc> Bye!\n",file_out,valid);

	return 0;
}

/*==============================================================================================*/
void parseOptions(int argc, char** argv)
{
	std::string temp;
	CmdLine cmd("imc","Monte-Carlo program based on normal modes.", version );

	try {
		// Define labeled arguments

		SwitchArg Verb("", "verb","Enables verbose.", true);
		cmd.add( Verb );

		// DEVELOPER INPUT
        ValueArg<std::string> FixIC("","fixIC", "Plain-text file defining the fixed Internal Coordinates. "
	    		"Each line will contain the index (0,1,...) of the ICs to be removed "
	    		"(DEVELOPER's).",false,"fixstring","string");
        cmd.add( FixIC );
        ValueArg<int> MovType("","mov", "Motion Type (default=2): "
				"0=K-matrix, 1=V/W-arrays, 2=Simple-Rotations, 3=Linear (if Cartesian modes) (DEVELOPER's).",false,2,"int");
        cmd.add( MovType );
        // END DEVELOPER INPUT

		ValueArg<double> Rfactor("","rfact", "Agressivity factor (default=7.778) (DEVELOPER's).",false,7.778174594,"double");
		cmd.add( Rfactor );
		ValueArg<float> OptF("","optf", "Factor to scale the optimal energy/stiffness factor (default=disabled) (DEVELOPER's).",false,1.0,"float");
		cmd.add( OptF );
		SwitchArg Weight("", "unweight","Un-mass-weights the input vectors (default=false). The Mass-weighted modes (wcart) will be converted into Cartesian coordiantes (cart) (DEVELOPER's).", false);
		cmd.add( Weight );

		SwitchArg Var("", "var","Input eigenvalues will be considered as variance (pca), otherwise the will be force constants (nma) (default=false).", false);
		cmd.add( Var );
		SwitchArg FirstWrite("", "include_first","Includes input model as first frame in the Multi-PDB trajectory (default=disabled).", false);
		cmd.add( FirstWrite );
        ValueArg<unsigned int> Seed("","seed", "Set the random number generator SEED (Mersenne Twister) (default=random-seed from /dev/urandom).",false,386,"unsigned int");
		cmd.add( Seed );
		SwitchArg Cart("","cart", "Mandatory if Cartesian modes are supplied, otherwise moving in Internal Coordinates (default=disabled).", true);
		cmd.add( Cart );
		ValueArg<int> TrajFormat("","otraj", "Output trajectory format: 0-Normal-Mode, 1-Multi-PDB, 2-AMBER (default=1).",false,1,"int");
		cmd.add( TrajFormat );
		ValueArg<float> Rmsd("","Rmsd", "Filter models by target RMSD (default=disabled).",false,0.0,"float");
		cmd.add( Rmsd );
		ValueArg<float> Thr("","thr", "Enable filtering by absolute tolerance (default=disabled).",false,0.1,"float");
		cmd.add( Thr );
		ValueArg<float> Ratio("","perc", "Enable filtering by tolerance percentage (default=1, i.e. 1%).",false,1,"float");
		cmd.add( Ratio );
		ValueArg<double> AmpF("a","amp", "Amplitude linear factor to scale motion (default=1).",false,1,"float");
		cmd.add( AmpF );
		ValueArg<double> Temp("T","temperature", "Temperature [K] (default=300).",false,300,"double");
		cmd.add( Temp );
		SwitchArg Chi("x","chi", "Considers first CHI dihedral angle (default=disabled).", true);
		cmd.add( Chi );
		ValueArg<int> Iter("i","iter", "Number of MC iterations per output structure (default=1000).",false,1000,"int");
		cmd.add( Iter );
		ValueArg<double> Efact("E","efact", "Energy/Stiffness scaling factor (mode energy will be multiplied by this value) (default=1.0).",false,1,"double");
		cmd.add( Efact );
		ValueArg<int> NumFrames("c","frames", "Number of output conformations (default= 100).",false,100,"int");
		cmd.add( NumFrames );
        ValueArg<float> Nevs("n","nevs", "Number of eigenvectors to be employed, either number [1,N] <integer>, or ratio from maximum available [0,1) <float> "
        		"(default=5).",false,5,"int");
        cmd.add( Nevs );
		ValueArg<float> Rg("g","Rg", "Filter models by target radius of gyration (default=disabled).",false,0.0,"float");
		cmd.add( Rg );
		ValueArg<float> OptRmsd("p","optRmsd", "Finds the optimal energy/stiffness scaling factor to obtain the desired average RMSD (A) from the initial model (default=disabled).",false,0.0,"float");
		cmd.add( OptRmsd );
	    ValueArg<std::string> FixFile("f","fixFile", "Input ASCII file defining the ICs that were fixed during NMA (default=disabled). "
	    		"If modes were computed removing arbitrary ICs, the user must introduce here the file generated by iMode's --save_fixfile option.",false,"fixstring","string");
        cmd.add( FixFile );
        ValueArg<std::string> Name("o","name", "Output files basename (default=imc).",false,"imc","string");
        cmd.add( Name );
//		ValueArg<int> Model("m","model", "Input Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom, 3=NCAC, 4=CA-only (default=2).",false,2,"int");
		ValueArg<int> Model("m","model", "Input Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom (default=2).",false,2,"int");
        cmd.add( Model );

        if(Verb.isSet())
		{
			verb_switch = true;
			printf("imc> Verbose enabled\n");
		}

		// Define required arguments not labeled
		UnlabeledValueArg<std::string> Pdb("pdb1","PDB input file. ","default","pdb");
		cmd.add( Pdb );
		UnlabeledValueArg<std::string> Evec("ptraj1","Normal modes input file (.evec), either from NMA or PCA.","default","ptraj");
		cmd.add( Evec );

		// Parse the command line.
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		strcpy(file_pdb,((temp=Pdb.getValue()).c_str())); // Gets PDB file name
		strcpy(file_ptraj,((temp=Evec.getValue()).c_str())); // Gets Ptraj file name
		strcpy(name,((temp=Name.getValue()).c_str())); // Gets Base-Name

    	// Number of eigenvectors to be computed
		if(Nevs.isSet())
		{
			nevec_fact = Nevs.getValue();
			if(nevec_fact <= 0) // checking
			{
				printf("imc> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
				exit(1);
			}
		}

		nframes= NumFrames.getValue();
		nmovs= Iter.getValue();
		traj_format= TrajFormat.getValue();
		if(MovType.isSet())
			mov_method = MovType.getValue();
		else
			if(Cart.isSet())
				mov_method = 3; // lineal motion
			else
				mov_method = 2; // simple rotations motion

		// Setting model and chi
    	model = Model.getValue();
		if(Chi.isSet())
			type = 2; // phi,chi,psi
		else
			type = 0; // phi,psi
        Ta = Temp.getValue();
        Ef = Efact.getValue();
        rfactor = Rfactor.getValue();
		linear_factor = AmpF.getValue();

        if(FixIC.isSet())
    	{
    		strcpy(fix_file,((temp=FixIC.getValue()).c_str())); // Gets Fix-file
    		fixmodel = 2;
    	}
        if(FixFile.isSet())
    	{
    		strcpy(fix_file,((temp=FixFile.getValue()).c_str())); // Gets Fix-file
    		fixmodel = 3;
    	}

        if(FirstWrite.isSet())
        {
    		printf("imc> Input PDB will be part of the Multi-PDB.\n");
        	firstwrite_switch = true;
        }

        if(OptRmsd.isSet())
        {
        	opt_switch=true;
        	optrmsd = OptRmsd.getValue();
        	printf("Parser input: Enabled Energy/Stiffness optimization to match RMSD= %f\n",optrmsd);

        	if(OptF.isSet())
        	{
        		optscale = OptF.getValue();
        		printf("Parser input: Factor to scale the optimal energy/stiffness factor: %f\n",optscale);
        	}
        }

        if(Rg.isSet())
        {
			filter_switch = true;
        	rg_switch = true;
			tscore = Rg.getValue();
			printf("Parser input: Filtering by Target Radius of gyration = %fA\n",tscore);
        }

        if(Rmsd.isSet())
        {
			filter_switch = true;
        	rmsd_switch = true;
			tscore = Rmsd.getValue();
			printf("Parser input: Filtering by Target RMSD = %fA\n",tscore);
        }

    	scorethr = Thr.getValue(); // default (always it may be needed)
        if(Thr.isSet()) // threshold filtering defined by absolute tolerance
        {
            filter_switch = true;
        	printf("Parser input: Absolute tolerance = %fA\n",scorethr);
        }

        if(Ratio.isSet()) // threshold filtering defined by tolerance ratio
        {
            filter_switch = true;
        	scorerat = Ratio.getValue();
        	scorethr = tscore * scorerat / 100;
        	printf("Parser input: Tolerance percentage = %.2f (%fA)\n",scorerat,scorethr);
        }

		if(Seed.isSet()) // Fixed seed
			seed = (unsigned int) Seed.getValue(); // Gets seed
		else // Random seed (time initialization)
//			seed = (unsigned int) time(0); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)
//			seed = (unsigned int) clock(); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)
		{     // Needed to avoid seed repetition between different runs.
			  FILE *fp;
			  unsigned char b[4];
			  int l=0;
			  if ((fp = fopen("/dev/urandom", "r")) == NULL) {
			    fprintf(stderr, "Error! Could not open /dev/urandom for read\n");
			    exit(2);
			  }
			  fread(b,1,4,fp);
			  l |= b[0] & 0xFF;
			  l <<= 8;
			  l |= b[1] & 0xFF;
			  l <<= 8;
			  l |= b[2] & 0xFF;
			  l <<= 8;
			  l |= b[3] & 0xFF;
			  seed = (unsigned int) l;
			  fclose(fp);
		}
		if(verb_switch)
			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);

		if(Weight.isSet())
		{
			weight_switch = true;
			printf("imc> Input eigenvectors/values will be considered as mass weighted\n");
		}
		if(Var.isSet())
		{
			var_switch = true;
			printf("imc> Input eigenvalues will be considered variance (from PCA)\n");
		}

	} catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}
