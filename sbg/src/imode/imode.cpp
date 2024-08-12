/*************************************************************************
 *                            iMODE                                      *
 *************************************************************************
 * This program is part of iMOD: http://chaconlab.org/imod/index.html    *
 * (c) Jose Ramon Lopez-Blanco, Jose Ignacio Garzon and Pablo Chacon.    *
 * IQFR-CSIC's Structural Bioinformatics Group (2004-2011).              *
 *************************************************************************
 *                                                                       *
 *   Tool for NMA computation in Internal Coordinate Space (ICS)         *
 *   using different Coarse-Graining models.                             *
 *   It takes into account Multiple-Chains, and outputs both Internal     *
 *   and Cartesian Coordinates Normal Modes.                             *
 *                                                                       *
 *************************************************************************
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <stdio.h> // needed in 48Gb-ubuntuserver
#include "cmdl/CmdLine.h"
#include "libnma/include/libnma_misc.h" // Mon's NMA Internal Coordinates related library
#include "libnma/include/libnma_io.h" // Mon's NMA Input-Output library
#include "libnma/include/libnma_cg.h" // Mon's NMA Coarse-Graining library
#include "libnma/include/libnma_deriv.h" // Mon's NMA Derivatives library
#include "libnma/include/libnma_kinetic.h" // Mon's NMA Kinetic Energy Matrix library
#include "libnma/include/libnma_hessian.h" // Mon's NMA Hessian Matrix library
#include "libnma/include/libnma_diag.h" // Mon's LAPACK Matrix Diagonalization calls
#include "libnma/include/libnma_def.h" // Pablo's Deformability library
#include "libnma/include/libnma_time.h" // Real-timer (Santi's)


/*==============================================================================================*/
char version[]="1.20"; // version code
// Input variables
char file_pdb[FILE_NAME];
char Kfile[FILE_NAME]; // Kfile input
char name[FILE_NAME];
char text[FILE_NAME];
char file_ptraj[FILE_NAME];
char file_ptraj_cart[FILE_NAME];
char dummy_string[FILE_NAME];
char fix_file[FILE_NAME];
char file_ss[FILE_NAME];
char file_func[FILE_NAME];
char fix_ss[30]; // string which allocates single character SS identifiers
char saved_files[2000]; // string to buffer screen output until the program end.
char file_inevec[FILE_NAME]; // Input eigenvectors file name

// INPUT PARAMETERS (CUSTOMIZABLE by parser)
int contacts = 0; // Contacting method
int model = 0; // Model Coarse-Graining
int type = 0; // Internal Coordinates Coarse-Graining type
int nthreads = 0; // Number of threads used by move_vwMFA_par()
int eigensolver = 1; // Eigensolver choice: 0= dspgvx, 1= dsdrv1_AP_BP_W_mon
bool SSout_switch = false; // outputs secondary structure (if computed internally)
bool Kout_switch = false; // outputs contacts
bool nomass_switch = false; // true --> Constant mass for every atom (=1.0)
bool notors_switch = false; // true --> Disables TORSional springs addition (hessianMDHx)
bool norm_modes = false; // true = Normalizes modes (norm=1)
bool save_cart = false; // true --> compute cartesian modes
bool save_cart2 = false; // true --> save cartesian modes
bool save_wcart = false; // true --> compute mass-weighted cartesian modes
bool save_wcart2 = false; // true --> save mass-weighted cartesian modes
bool save_covar = false; // true --> save covariance matrix
bool save_covarf = false; // true --> save output model covariance matrix
bool save_dcovar = false; // true --> save CA-based distance-covariance matrix
bool save_dh = true; // true --> save Internal Coordinates modes
bool save_matrices = false; // true --> save hessian and kinetic energy matrices
bool just_matrices = false; // true --> disables matrices memory allocation and forces exit after calculations.
bool save_matrices_text = false; // false --> saving matrices in text format (otherwise binary)
bool save_fixfile = false; // true --> save fixation mask file
bool deform_switch=false; // = true --> enable deformability computations
bool intermolec_switch = false; // Enables inter-molecular force constant factor
bool inevec_switch = false; // Disables Hess and Kinetic energy matrices calculation and diagonalization. Some input eigenvectors must be provided!
bool squared_matrices = false; // true --> Using squared matrices instead of the memory-efficient packed triangular storage.
bool delHydrogens_switch = true; // Delete hydrogens
bool delHeteros_switch = false; // Delete heteroatoms
bool delWaters_switch = true; // Delete waters

float cte_k0 = 1.0;
float cte_k1 = 1.0;
float x0 = 3.8;
float power = 6.0;
float cutoff_k0 = 10;
float cutoff_k1 = 10;
float cutoff_k2 = 10;
float swapmatrix_mem = 1; // RAM memory used during Hessian matrix swapping... (in GB)
float intermolec_factor = 1; // Inter-molecular force constant factor
double Ta; // Temperature [T]
double dc; // Characteristic distance factor for Deformability computations (see: hardy and compute_def functions in libnma_def.cpp)
int typef = 0; // Internal coordinates type (output type of normal modes)
int modelf = -1; // Internal Coordinates model (output model of normal modes)

// INPUT PARAMETERS (PRE-DEFINED)
int debug_code=0; // Debugging code
int conv_method=1; // Conversion method from IC into Cartesian coordinates (0=K-matrix, 1=VW-arrays)
int hessian_method=2; // Hessian matrix computation method (0=K-matrix, 1=VW-arrays, 2=Fast)
int kinetic_method=2; // Kinetic energy matrix computation method (0=K-matrix, 1=VW-arrays, 2=Fast)
bool saveformat_switch = false; // = true --> save formated input PDB
bool savemodel_switch = true; // = true --> save 3BB2R-model input PDB
bool ss_switch = false; // the NMA force constants will be set according to SS rules
bool func_switch = false; // input file with function coefficients (Topology and SS)
bool parse_verb = false;

// Some constants
float cte = 1.0;
int fixmodel = 0; // =0 --> no ICs fixation
float nevec_fact = -1; // % of eigenvectors to be computed (set by parser)
int nevec = -1; // number of eigenvectors to be computed (set by parser)
int modes_saved = 20; // now controled by "nevs" (if there are enougth modes available)

bool save_CA = false; // true --> compute cartesian CA modes
bool save_CA2 = false; // true --> save cartesian CA modes
bool nomodel = false; // true --> Avoids CG-model building from input model.
float fixRand_prob=0.5; // random IC fixation probability
bool debug = false;
int verb = 0;
int saved_files_len; // store the "saved_files" array current length
unsigned int seed; // random number generator seed
float change_timer = 300; // Time threshold (seconds) to use: timer, instead of Htimer

/*==============================================================================================*/
using namespace TCLAP;
void parseOptions(int argc, char** argv);
/*==============================================================================================*/
int main( int argc, char * argv[] )
{
	printf("imode>\nimode> Welcome to the NMA in Internal Coordinates tool v%s\nimode>\n",version);

	// Parsing Input
	parseOptions(argc,argv);

	// Mersenne Twister Seed Initialization
	// Output random float number in the interval 0 <= x < 1, with Random()
	rg = new CRandomMersenne( seed );

	double Kb = 0.00198717; // Kboltz (in Angstroms)
	double KbT = Kb * Ta;  // Kboltz * Tempereature (thermal-energy)

	if(verb > 0)
		printf("imode> Mersenne's twister seed = %d\n",seed);

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

	// Saving Input Command-Log File
	FILE *f_com;
	sprintf(text,"%s.log",name);
	if( !(f_com=(FILE *)fopen(text,"w") ) )
	{
		printf("Sorry, unable to open LOG FILE: %s\n",text);
		exit(1);
	}
	sprintf(saved_files,"imode> Log file:                            %35s\n",text);

	bool already_seed=false;
	for(int i=0; i<argc; i++)
	{
		fprintf(f_com,"%s ",argv[i]);
		if(strcmp(argv[i],"--seed") == 0 )
			already_seed = true;
	}
	if(!already_seed)
		fprintf(f_com,"--seed %u\n",seed); // This allows user to carry out again the same run.
	fclose(f_com);

	twid *decint; // contacts = 1
	int nipa; // contacts = 1
	double temp;
	Residue *res;
	Atom *at;
	timerReal timer; // Real timer (independent of parallelization)

	// Initialize aminoacids and nucleotids
	init_aminoacids();

	// READING INPUT PDB
	printf( "imode> Reading PDB file: %s\n",file_pdb );
	Macromolecule * molr = new Macromolecule( file_pdb );
	molr->readPDB(file_pdb);
	if(delHydrogens_switch)
	{
		printf( "imode> Deleting Hydrogen atoms (if any)...\n" );
		molr->deleteHYDS();
	}
	if(delHeteros_switch || model == 0) // CA model can't deal with HETATM's... (TO DO)
	{
		printf( "imode> Deleting Hetero-atoms (if any)...\n" );
		molr->delete_heteros();
	}
	if(delWaters_switch)
	{
		printf( "imode> Deleting Water molecules (if any)...\n" );
		molr->delete_waters();
	}
	molr->info(stdout);

	// Formating PDB first (if necessary)
	if(!nomodel)
	{
		if(verb > 1)
			printf( "imode> Formatting residues order and checking for missing atoms\n" );
		if(molr->format_residues(false,model) > 0)
		{
			printf( "imode> Error, missing atom(s) found! Forcing exit!\n" );
			exit(1);
		}
		if(saveformat_switch)
		{
			sprintf(dummy_string,"%s_format.pdb",name);
			molr->writePDB( dummy_string );
			sprintf(text,"imode> Formated model PDB:                  %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
		}
	}

	Macromolecule *mol,*molNCAC,*molini;
	pdbIter *iter;
	pdbIter *iterA,*iter2;
	float imass;
	int num_res = 0, num_atoms = 0, num_atomsNCAC = 0;

	// Needed further to convert normal modes into the appropriate atomic model
	molini = new Macromolecule(molr); // Initial (readed and formatted) macromolecule copy

	// Setting Coarse-Graining model ("mol" will hold current CG'ed Macromol.)
	switch(model)
	{
	case 0: // CA-IC model: CA + (NH and CO)-terminal model
	{
		printf( "imode> Coarse-Graining model: CA-model\n");

		// N,CA,C selection
		molNCAC = molr->select( ncac2 );
		num_atomsNCAC = molNCAC->get_num_atoms();


		// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
		// Warning, atoms are not copied, they're just pointers to the original atoms.
		// setmass = true --> adds masses to Occupancy and Bfactor, otherwise left unchanged.
		// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
		//                   (unit mass is divided between NH or CO and their CAs)
		// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
		//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)

		// Makes model (if it's needed)
		mol = cg_CA( molNCAC, !nomodel, nomass_switch ); // masses will be computed

		// Saving "NCAC" model
		if(savemodel_switch)
		{
			sprintf(dummy_string,"%s_ncac.pdb",name);
			molNCAC->writePDB( dummy_string );
			sprintf(text,"imode> NCAC model PDB:                      %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
		}
		break;
	}
	case 3: // N,CA,C-model
	{
		printf( "imode> Coarse-Graining model: N,CA,C-model (Experimental)\n");

		// N,CA,C selection
		mol = molr->select( ncac2 );

		if(!nomodel) // Makes model (if it's needed)
			mass_NCAC( mol, nomass_switch );
		break;
	}
	case 1: // 3BB2R model
	{
		printf( "imode> Coarse-Graining model: 3BB2R\n");
		mol = molr;
		if(!nomodel) // Makes 3BB2R model (if it's needed)
		{
			// CREATES a 3BB2R reduced model
			//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
			//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
			//     There are a few exceptions: for Ala, Gly and Pro,
			//     and for the 1st and last residues.
			printf("imode> Creating 3BB2R reduced model:\n");
			cg_3BBR2(mol, nomass_switch, nomass_switch);
		}
		break;
	}
	case 2: // Full-Atom
	{
		printf( "imode> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
		mol = molr;
		if(!nomodel)
			mass_FA(mol, nomass_switch); // sets masses
		break;
	}
	}

	// Saving current model
	if(savemodel_switch)
	{
		sprintf(dummy_string,"%s_model.pdb",name);
//		mol->writePDB( dummy_string, false ); // renumbers the PDB
		mol->writePDB( dummy_string ); // NO-renumber the PDB
		sprintf(text,"imode> Model PDB:                           %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
	}

	int resn, j = 0, i = 0;
	// Initializing some internal coords. related residue properties
	// (# atoms, # units, # internal coordinates, etc...)
	tri *props;
	int *unat;
	switch(model)
	{
	case 0:
	case 3:
		properCA(mol,&props,&unat);
		break;
	case 1:
	case 2:
		properMFA(mol,&props,&unat,type,model);
		break;
	}

	// Computing Total number of degrees of freedom (hessian matrix rank)
	int size = 0;
	iter = new pdbIter( mol, true, true, true, true ); // iter to screen fragments (residues)

//	mol->info(stderr);

	num_atoms = iter->num_atom();
	num_res = iter->num_fragment();
	printf( "imode> Selected model residues: %d\n", num_res );
	printf( "imode> Selected model (pseudo)atoms: %d\n", num_atoms );

	// Theoric number of DoFs= 3*T + 3*R -6 + Dihedral
	// (Before fixing, i.e. the fix-file format DoFs...)
	int old_size,seg_atoms=0;
	size = 0;
	old_size = size; // temp (Dihedral ICs)
	for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		size += props[iter->pos_fragment].nan; // Each residue may have different number of dihedral-angles!
	printf( "imode> Number of Dihedral angles: %d\n", size );
	old_size = size;
	for( iter->pos_segment = 0; !iter->gend_segment(); iter->next_segment() ) // screen segments
	{
		seg_atoms = (iter->get_segment())->num_atoms();
		if(seg_atoms > 1)
			size += 6;
		else
		{
//			fprintf(stderr, "Single atom segment found! seg= %d\n", iter->pos_segment);
			size += 3;
		}
	}
	printf( "imode> Rotational/Translational ICs (Non-Eckart): %d\n", size-old_size );
	size -= 6; // Eckart conditions
	printf( "imode> Predicted number of ICs (Eckart): %d\n", size );

	char *ss_table; // SS-table
	TSfunc *funcs; // SS and Topology functions
	int nfunc=0; // number of functions
	int nss=0;

	if(ss_switch) // NMA with SS (optional)
	{
		// Reading ss-file (allocating table memory)
		read_ss(file_ss,&ss_table,&nss);
		// Some checking...
		if( nss != num_res)
		{
			printf("imode> Sorry... the SS-file should have the same number of entries (%d) as the number of residues (%d) in the input PDB\n",nss,num_res);
			exit(1);
		}
	}
	else if(contacts == 3 || fixmodel == 6)
		ss_table = mol->secondary_structure(true); // "make_ipasTS" should work as well...
	else
		ss_table = NULL; // disables SS-checking in "make_ipasTS" (Topology only)

	if(SSout_switch && (ss_switch || contacts == 3 || fixmodel == 6))
	{
		sprintf(dummy_string,"%s.ss",name);
		write_SSfile(ss_table,num_res,dummy_string); // Write Secondary Structure file
		sprintf(text,"imode> Secondary Structure file:            %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		if(verb > 1)
		{
			printf("imode> Computed Secondary Structure:\n");
			for(int i=0;i<num_res;i++)
				printf("%c ",ss_table[i]);
			printf("\n");
		}
	}

	// Reading Secondary structure file and T/SS Functions.
	if(contacts == 3)
	{
		// Reads input TS functions (mandatory)
		if(func_switch)
			read_TSfunc(file_func, &funcs, &nfunc);
		else
		{
			printf("imode> Sorry... At this moment you should include always a connection functions file\n");
			exit(1);
		}

		if(verb > 0)
			for(int i=0; i<nfunc; i++)
				printf("imode> %c%c %d %f %f %f\n", funcs[i].i, funcs[i].j, funcs[i].t, funcs[i].a, funcs[i].b, funcs[i].c );
	}

	// FREEZING VARIABLES
	// Mobile Internal Coordinates selection "i"-model
	// Fixation array is related to "i"-model, not to "current" one!
	// input: fixedi ---> output: fixedx
	bool *fixed=NULL;
	if(fixmodel != 0)
	{
		save_fixfile = true; // By default, it saves the fixing mask... (otherwise user will not be able to animate modes)
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
		case 4:
			size = fixRand(fixed,size,fixRand_prob);
			break;
		case 5:
			size = fixRandDH(mol,props,fixed,type,model,fixRand_prob);
			break;
		case 6:
			if(!ss_switch) // if not specified a SS file, it will be computed
				ss_table = mol->secondary_structure(true);
			size = fixSS(mol,props,fixed,ss_table,fix_ss,type,model);
			break;
		default:
			printf("Sorry, unknown fixation method!\nForcing exit\n");
			exit(1);
			break;
		}
		printf("imode> Input CG-model Fixed Internal Coordinates: %d\n", old_size-size );
		printf("imode> Input CG-model Final Internal Coordinates (sizei) = %d\n",size);
	}

	// CHECKING SINGLE RESIDUE SEGMENTS...
	if(model == 0 || model == 3)
	{
		pdbIter *iter_frag;
		for( iter->pos_segment = 0; !iter->gend_segment(); iter->next_segment() ) // screen segments
		{
			iter_frag = new pdbIter( iter->get_segment() ); // iters current segment
			if(iter_frag->num_fragment() == 1) // If single residue segment...
			{
				printf( "Error! Found segment(s) with only one residue!\n");
				exit(1);
			}
			delete iter_frag;
		}
	}

	if(save_fixfile)
	{
		sprintf(dummy_string,"%s.fix",name);
		write_fixDH(dummy_string,mol,props,fixed,type,model); // Output fixation mask
		if(verb > 1)
			printf("imode> Fixation mask written: %s\n",dummy_string);
		sprintf(text,"imode> Fixation mask file:                  %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
	}

	// Creates two auxiliar arrays with segment properties (needed due to fixing):
	//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
	//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
	// (Allocates memory itself, if it's needed)
	bool *addrot=NULL; // Should be added ROTATIONs? (with fixing)
	int *effseg=NULL; // Effective segment indices
	size = seg_props(mol, props, fixed, model, type, &addrot, &effseg);
	printf( "imode> Number of ICs predicted by seg_props(): %d\n", size );

	if(size<=0)
	{
		printf("imode> Sorry invalid size (%d), forcing exit!\n",size);
		exit(1);
	}

	// Number of eigenvectors to be computed (we need to know "size" first)
	if(nevec_fact >= 1.0) // number of modes
	{
		nevec = (int) nevec_fact;
		if(!inevec_switch)
			printf( "imode> Range of computed modes: 1-%d\n", nevec);
		else
			printf( "imode> Range of processed modes: 1-%d\n", nevec);
	}
	else
	{
		nevec = (int) (nevec_fact * size);
		if(!inevec_switch)
			printf( "imode> Range of computed modes: 1-%d (%.0f%)\n", nevec, nevec_fact*100);
		else
			printf( "imode> Range of processed modes: 1-%d\n", nevec);
	}
	// Checking
	if(nevec > size)
	{
		printf("imode> Sorry, more eigenvectors requested (%d) than available (%d), forcing maximum.\n",nevec,size);
		nevec = size;
	}
	else if(nevec <= 0) // checking
	{
		printf("imode> Error, invalid number of eigenvectors requested %d (%f)!\nForcing exit!\n",nevec,nevec_fact);
		exit(1);
	}

	// You only save what you compute!
	modes_saved = nevec;
	if(verb > 1)
		printf( "imode> Number of Saved modes: %d\n", modes_saved);

	// Getting coordinates single row (pseudo-atom model)
	float *coord,*coordNCAC;
	if(verb > 1)
		printf ("imode> Getting coordinates single row (pseudo-atom model)\n");
	mol->coordMatrix( &coord );
	if(model==0)
		molNCAC->coordMatrix( &coordNCAC );

	// Shifting coordinates to their Center of Mass (Computing the CoM)
	// The structure should be centered on its CoM in order
	// to compute Kinetic-energy and Hessian Matrices.
	iter2 = new pdbIter( mol ); // Iterator to screen atoms
	double mtot = 0.0;
	double mta = 0.0;
	double rd[3];
	rd[0] = rd[1] = rd[2] = 0.0;
	for ( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )  // screens all-atoms
	{
		mta = ( iter2->get_atom() )->getPdbocc(); // Load mass from occupancies...
		mtot += mta;
		// Sum(mass*coord) before putting the CoM at 0
		rd[0] += mta * coord[iter2->pos_atom * 3];
		rd[1] += mta * coord[iter2->pos_atom * 3 + 1];
		rd[2] += mta * coord[iter2->pos_atom * 3 + 2];
	}
	delete(iter2);
	rd[0] /= mtot;
	rd[1] /= mtot;
	rd[2] /= mtot;
	if(verb > 1)
		printf( "imode> Mass %12.8f Center %.8f %.8f %.8f --> to 0,0,0\n", mtot, rd[0], rd[1], rd[2] );

	// Shifting CoM into origin
	for(int k = 0; k < num_atoms; k++)
	{
		coord[k * 3] -= rd[0];
		coord[k * 3 + 1] -= rd[1];
		coord[k * 3 + 2] -= rd[2];
	}

	// Shifting CoM of NCAC-model into origin (needed by naive-derivatives computation)
	if(model==0)
		for(int k = 0; k < num_atomsNCAC; k++)
		{
			coordNCAC[k * 3] -= rd[0];
			coordNCAC[k * 3 + 1] -= rd[1];
			coordNCAC[k * 3 + 2] -= rd[2];
		}

	// ******************************************
	// * CONTACTING
	// ******************************************
	FILE *f_Kout;
	if(!inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
	{
		printf("imode> Creating pairwise interaction potential:\n");
		// "decint" stores the ipas-list
		if( !(decint = ( twid * ) malloc( 1 * sizeof( twid ) ) ) )
		{
			printf("Sorry, \"decint\" memory allocation failed!\n");
			exit(1);
		}
		timer.startTimer(); // timer
		switch(contacts)
		{
		case 0: // INVERSE EXPONENTIAL (power of distance for contact matrix)
		{
			// Making Interacting Pair of (non-virtual) Atoms (ipas)
			make_ipas_new0(mol, &decint, &nipa, (float) cutoff_k0 );
			printf("imode> Inverse Exponential (%d nipas) cutoff= %.1f, k= %f, x0= %.1f ",nipa,cutoff_k0,cte_k0,x0);
			for(int i=0; i<nipa; i++)
				decint[i].C = inv_exp(cte_k0,decint[i].d,x0,power); // setting Force Constants
			break;
		}

		case 1: // DISTANCE CUTOFF METHOD
		{
			// Making Interacting Pair of (non-virtual) Atoms (ipas)
			make_ipas_new0( mol, &decint, &nipa, (float) cutoff_k1 );
			printf("imode> Cutoff Distance (%d nipas) cutoff= %.1f, k= %f ",nipa,cutoff_k1,cte_k1);
			for(int i=0; i<nipa; i++)
				decint[i].C = cte_k1; // setting Force Constants
			break;
		}

		case 2:	// HINSEN DISTANCE METHOD
		{
			printf("imode> Hinsen's distance criterion ");
			for(int i=0; i<nipa; i++)
			{
				if(decint[i].d <= 4.0)
					decint[i].C = 86000*decint[i].d-23900; // setting Force Constants
				else
					decint[i].C = 128/ (pow( decint[i].d/10, 6 ) );
			}
			break;
		}

		case 3: // Setting Force Constants according to Secondary Structure and Topology
		{
			// Reads input TS functions (mandatory)
			printf("imode> Topology and/or Secondary Structure (make_ipasTS) ");
			if(func_switch)
				read_TSfunc(file_func, &funcs, &nfunc);
			else
			{
				printf("imode> Sorry... At this moment you should include always a connection functions file\n");
				exit(1);
			}

			if(verb > 0)
			{
				for(int i=0; i<nfunc; i++)
					printf("imode> %c%c %d %f %f %f\n", funcs[i].i, funcs[i].j, funcs[i].t, funcs[i].a, funcs[i].b, funcs[i].c );
			}
			make_ipasTS(mol, &decint, &nipa, cutoff_k2, funcs, nfunc, ss_table);
			break;
		}

		case 4: // Laura's "Mixed" model
			// Making Interacting Pair of (non-virtual) Atoms (ipas)
			make_ipas_mix0(mol, &decint, &nipa);
			printf("imode> \"Mixed\" (%d nipas) ",nipa);
			break;

		case 5: // Taking contacts from a KLIST-file
		{
			// Reading contacts form a force constants file (Klist)
			printf("imode> Reading force constants file: %s",Kfile);
			// exit(0);
			// Read Force constants file (Kfile) allocating memory
			// Warning: if "coord" not provided, distances ".d" will not be updated!
			read_Kfile(&decint,&nipa, Kfile, coord);
			printf(" (%d nipas) ",nipa);
			break;
		}

		default:
			printf("Please, introduce a valid Connection method to continue!!!\n\nForcing exit!\n\n");
			exit(1);
			break;
		}
		timer.stopTimer();
		printf("(%lf s)\n",timer.getElapsedTime());
	}

	if(intermolec_switch)
	{
		printf("\nimode> Multiply each inter-molecule force constant by factor= %f\n",intermolec_factor);
		// Multiplies by "factor" every ipa force constant ("C") if both atoms belong to different molecules.
		modify_intermolec_ipas(mol, decint, nipa, intermolec_factor);
	}

//	// Sorting "ipas" to facilitate Hessian computation (mandatory for huge systems!)
//	sort_ipas(decint,nipa,unat);

	// Write Force constants file (Kfile)
	if(Kout_switch && !inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
	{
		sprintf(dummy_string,"%s_Kfile.dat",name);
		write_Kfile(decint,nipa,dummy_string);
		if(verb > 1)
			printf( "imode> Contacts found: %d --> Written to %s\n", nipa, dummy_string);
		sprintf(text,"imode> Force constants file:                %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
	}

	// IPAs checking
	if(verb > 10 && !inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
		for(int i=0;i<nipa;i++)
			printf("ipa %4d: k= %d  l= %d  d= %f  C= %f\n",i,decint[i].k,decint[i].l,decint[i].d,decint[i].C);

	// K-matrix computation (K-matrix)
	trd *der = NULL; // = NULL enables memory allocation!
	trd *derNCAC = NULL; // = NULL enables memory allocation!
	double ***V=NULL, ***W=NULL;
	int **body1=NULL;

	// Derivatives (dy/dq) array computation
	if(!inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
		if( hessian_method == 0 || hessian_method == 10 || kinetic_method == 0 || kinetic_method == 10 )
		{
			timer.startTimer(); // timer
			switch(model)
			{
			case 0: // (first-N, CAs, last-C model)
				printf("imode> Computing CA-only K-matrix with dydqMCAx() ");
				dydqMCAx(mol,coordNCAC, &der, size, fixed);
				break;
			case 3: // (N,CA,C-model)
				printf("imode> Computing CA-only K-matrix with dydqMCA3x() ");
				dydqMCA3x(mol,coord, &der, size, fixed);
				break;
			case 1: // 3BB2R model
			case 2: // Full-Atom model
				printf("imode> Computing Full-Atom/3BB2R K-matrix with dydqMFAx() ");
				dydqMFAx(mol,coord, props, &der, type, model, size, fixed, addrot);
				break;
			}
			timer.stopTimer();
			printf("%lf s\n",timer.getElapsedTime());
		}

	if( hessian_method == 1 || hessian_method == 11 || hessian_method == 12 || kinetic_method == 1 || kinetic_method == 11 )
	{
		timer.startTimer(); // timer
		switch(model)
		{
		case 0: // (first-N, CAs, last-C model)
			//			printf("imode> Computing CA-only K-matrix with vwMCAx() ");
			//			vwMCAx(mol,coordNCAC,&V,&W,&body1,size,fixed);
			printf("imode> Computing CA-only K-matrix with vwMCA3x() ");
			vwMCA3x(molNCAC,coordNCAC,&V,&W,&body1,size,fixed,model);
			break;
		case 3: // (N,CA,C-model)
			printf("imode> Computing CA-only K-matrix with vwMCA3x() ");
			vwMCA3x(mol,coord,&V,&W,&body1,size,fixed);
			break;
		case 1: // 3BB2R model
		case 2: // Full-Atom model
			printf("imode> Computing Full-Atom/3BB2R K-matrix with vwMFAx() ");
			vwMFAx(mol,coord,props,&V,&W,&body1,type,model,size,fixed,addrot);
			break;
		}
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());
	}

	// BUILDING HESSIAN MATRIX
	floating *hess_tr = NULL; // Hessian generation routines will allocate memory!
	floating *hess_sq = NULL; // To store the squared matrices instead of packed ones...

	char *matrixfile = NULL;

	if(!inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
	{
		timer.startTimer(); // timer
		printf("imode> Packed-Hessian/Kinetic matrices mem.= %.3f Mb (rank= %d)\n",((float)sizeof(double)*size*(size+1)/2)/1000000,size);

		// Saving matrix
		if(save_matrices)
		{
			sprintf(dummy_string,"%s_hess.bin",name);
			matrixfile = dummy_string;
		}

		switch( hessian_method )
		{
		case 0: // "naive" Kovacs' method (PRE-DEFINED PRECISION)
		case 3: // "naive" Kovacs' method (PRE-DEFINED PRECISION)
			printf("imode> Slow Hessian Matrix Building O(n^4) [hess_naive0()] ");
			fflush(stdout);
			hess_naive0(decint,der,nipa,size,coord,&hess_tr); // avoids zeros (x2 faster)
			// hess_naive(decint,der,nipa,size,coord,&hess_tr);
			break;
		case 10: // "naive" Kovacs' method (FORCING DOUBLE PRECISION)
		case 30: // "naive" Kovacs' method (FORCING DOUBLE PRECISION)
			printf("imode> Slow Hessian Matrix Building O(n^4) [hess_naive0_double()] ");
			fflush(stdout);
			hess_naive0_double(decint,der,nipa,size,coord,&hess_tr); // avoids zeros (x2 faster)
			break;
		case 1: // Memory efficient "naive" method (via V/W arrays) (PRE-DEFINED PRECISION)
			//		printf("imode> Slow Hessian Matrix Building O(n^4) [hess_naiveVW()]\n");
			printf("imode> Slow Hessian Matrix Building O(n^4) [hess_naiveVW0()] ");
			fflush(stdout);
			//		hess_naiveVW(coord,decint,V,W,body1,nipa,size,&hess_tr);
			hess_naiveVW0(coord,decint,V,W,body1,nipa,size,&hess_tr,model); // avoids zeros cool (x20)faster!!!
			break;
		case 11: // Memory efficient "naive" method (via V/W arrays) (FORCING DOUBLE PRECISION)
			printf("imode> Slow Hessian Matrix Building O(n^4) [hess_naiveVW()] ");
			fflush(stdout);
			hess_naiveVW(coord,decint,V,W,body1,nipa,size,&hess_tr,model);
			break;
		case 12: // Memory efficient "naive" method (via V/W arrays)
			printf("imode> Slow Hessian Matrix Building O(n^4) \"row by row\" [hess_naiveVW0()] ");
			fflush(stdout);
			//		hess_naiveVW(coord,decint,V,W,body1,nipa,size,&hess_tr);
			for(int row=0;row<size;row++)
				hess_naiveVW0(coord,decint,V,W,body1,nipa,size,&hess_tr,model,row); // avoids zeros cool (x20)faster!!!
			break;
		case 22:
			switch(model)
			{
			case 0:
				printf("imode> Low-mem & fast CA-model Hessian Matrix Building O(n^2) [hessianMCAxHD_par()] ");
				fflush(stdout);
				hessianMCAxHD_par(mol,decint,nipa,(long int)size,coordNCAC,coord,&hess_tr,props,unat,fixed,nthreads,matrixfile,just_matrices);
				break;
			default:
				exit(0);
			}
			break;

		case 2:
			switch(model)
			{
			case 0:
				//			printf("imode> Fast CA-model Hessian Matrix Building O(n^2) [hessianMCAx()] ");
				if(debug_code == 1)
				{
					printf("\ndebug> Pause before computations... (press any key to continue)\n");
					fflush(stdout);
					getchar();
				}
				//			hessianMCAx(mol,decint,nipa,size,coordNCAC,coord,&hess_tr,props,unat,fixed);
				printf("imode> Low-mem & fast CA-model Hessian Matrix Building O(n^2) [hessianMCAxHD()] ");
				fflush(stdout);
				hessianMCAxHD(mol,decint,nipa,(long int)size,coordNCAC,coord,&hess_tr,props,unat,fixed,matrixfile,just_matrices);
//				fprintf(stderr,"Saving matrix into hess_ref.bin\n");
//				save_matrixB(hess_tr,size,"hess_ref.bin");
//				show_matrix(hess_tr, size, "Hessian Matrix (F) hess_tr:");
//				exit(0);


				//			double *zero=NULL;
				//			int zero_size;
				//			read_matrixB2(&zero, &zero_size, matrixfile);
				//			show_matrix(zero, zero_size, matrixfile);
				//
				//			sprintf(dummy_string,"%s_hessOK.bin",name);
				//			save_matrixB(hess_tr,size,dummy_string);
				//			double *zero2=NULL;
				//			int zero2_size;
				//			read_matrixB(&zero2, &zero2_size, dummy_string);
				//			show_matrix(zero2, zero2_size, dummy_string);
				//
				//			sprintf(dummy_string,"%s_hess.bin",name);

				if(save_matrices)
				{
					timer.stopTimer();
					printf("%lf s\n",timer.getElapsedTime());

					// Remove this:
					printf("imode> Swapping Hessian matrix: ");
					fflush(stdout);
					timer.startTimer(); // timer

					swapSPmatrix(matrixfile, (long int) (pow(2,30)/8 * swapmatrix_mem)); // Swap matrix using stripes... (2GB)
				}

				if(debug_code == 1)
				{
					printf("\ndebug> Pause after computations... (press any key to continue)\n");
					fflush(stdout);
					getchar();
				}
				break;
			case 1:
			case 2:
				printf("imode> Fast Hessian Matrix Building O(n^2) [hessianMFAx()] ");
				fflush(stdout);
				hessianMFAx(mol,decint,nipa,size,coord,&hess_tr,props,unat,type,model,fixed,addrot);
				break;
			case 3:
				printf("imode> N,CA,C-model K-matrix with hessianMCA3x() NOT IMPLEMENTED YET\n"
						"Try \"naive method via K-matrix or V/W-arrays\" instead.\n");
				exit(1);
				break;
			}
			break;
		}

		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());
		free(decint); // free contact list

		// Adds torsional springs to a previously built hessian matrix (triangular packing)
		// (ec.5) from:   Lu, Poon and Ma. J.Chem. Theory Comput. 2006, 2, 464-471.
		if(!notors_switch)
		{
			if(!save_matrices)
				hessianMDHx(mol,props,hess_tr,1.0,size,fixed,addrot);
			else // if(save_matrices)
			{
				if(model==0)
				{
					fprintf(stderr,"\nUSING hessianMDHxHD(), copón!\n");
					hessianMDHxHD(mol,props,matrixfile,1.0,fixed,addrot);
					if(hess_tr != NULL && !just_matrices) // Hessian "update" should be also carried out in memory...
					{
						hessianMDHx(mol,props,hess_tr,1.0,size,fixed,addrot);
					}
				}
				else // Given "hessianMFAxHD" is not yet developed, there is no a hessian-matrix .bin file... (TO DO LIST)
					hessianMDHx(mol,props,hess_tr,1.0,size,fixed,addrot);
			}
		}

		// Saving matrix
		if( save_matrices && !(hessian_method == 2 && model == 0) )
		{
			timer.startTimer();
			save_matrixB(hess_tr,size,dummy_string);
			if(verb > 1)
			{
				timer.stopTimer();
				printf("imode> Saving Hessian (A) matrix in %s: %lf s\n",dummy_string,timer.getElapsedTime());
				fflush(stdout);
			}
		}
		if( save_matrices ) // Storing the file saving summary
		{
			sprintf(text,"imode> Hessian matrix binary file:          %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
		}

		if(verb > 2)
			show_matrix(hess_tr, size, "Hessian Matrix (F) hess_tr:");

		if(squared_matrices) // If symmetric squared matrices instead of packed triangular storage...
		{
			// Allocating single precision memory
			if( !(hess_sq = (double *) malloc( sizeof( double ) * size*size ) ) )
			{
				printf("Msg(diag): I'm sorry, unable to allocate %d bytes\nForcing exit\n",sizeof( double )*size*size);
				exit(1);
			}
			for ( i = 0; i < size; i++ )
				for ( j = i; j < size; j++ ) // diagonal included
				{
					hess_sq[i+size*j] = hess_tr[i + j*(j+1)/2]; // upper triangle + diagonal
					if(i != j) // avoids writing diagonal twice...
						hess_sq[j+size*i] = hess_tr[i + j*(j+1)/2]; // lower triangle
				}
			free(hess_tr); // not needed any more
		}
	}

	// BUILDING KINETIC ENERGY MATRIX
	floating *mass_tr = NULL; // triangular matrix input arrays
	floating *mass_sq = NULL; // To store the squared matrices instead of packed ones...

	if(!inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
	{
		timer.startTimer();

		if(save_matrices)
		{
			sprintf(dummy_string,"%s_mass.bin",name);
			matrixfile = dummy_string;
		}

		switch( kinetic_method )
		{
		case 0: // Slow Kinetic Energy matrix
		case 3: // Slow Kinetic Energy matrix
			printf ("imode> Slow Kinetic-Energy matrix Building O(n^3) [kinetic_naive()] ");
			fflush(stdout);
			kinetic_naive( mol, &mass_tr, der, size );
			break;
		case 10: // Slow Kinetic Energy matrix
		case 30: // Slow Kinetic Energy matrix
			printf ("imode> Slow Kinetic-Energy matrix Building O(n^3) [kinetic_naive_double()] ");
			fflush(stdout);
			kinetic_naive_double( mol, &mass_tr, der, size );
			break;
			//	case 1: // Slow (memory efficient) Kinetic Energy matrix
			//		printf ("imode> Slow Kinetic-Energy matrix Building O(n^3) [kinetic_naiveVW()] ");
			//		fflush(stdout);
			//		kinetic_naiveVW( mol,coord,&mass_tr,V,W,body1,size);
			//		break;
			//	case 11: // Slow (memory efficient) Kinetic Energy matrix (Forcing DOUBLE)
			//		printf ("imode> Slow Kinetic-Energy matrix Building O(n^3) [kinetic_naiveVW_double()] ");
			//		fflush(stdout);
			//		kinetic_naiveVW_double( mol,coord,&dmass_tr,V,W,body1,size);
			//		mass_tr = (floating *) malloc(sizeof(floating) * size*(size+1)/2);
			//		for(int i=0; i<size*(size+1)/2; i++)
			//			mass_tr[i] = (floating) dmass_tr[i];
			//		free(dmass_tr);
			//		break;
		case 2: // Fast Kinetic Energy matrix
			switch(model)
			{
			case 0:
				//			printf("imode> Fast CA-model Kinetic Matrix Building O(n^2) [kineticMCAx()] ");
				printf("imode> Fast CA-model Kinetic Matrix Building O(n^2) [kineticMCAxHD()] ");
				fflush(stdout);
				if(debug_code == 1)
				{
					printf("\ndebug> Pause before computations... (press any key to continue)\n");
					fflush(stdout);
					getchar();
				}
				//			kineticMCAx( mol, coordNCAC, props, size, &mass_tr, model,fixed ); // H-matrix (triangular)
				kineticMCAxHD( mol, coordNCAC, props, size, &mass_tr, model,fixed,matrixfile,just_matrices ); // H-matrix (triangular)
				if(debug_code == 1)
				{
					printf("\ndebug> Pause after computations... (press any key to continue)\n");
					fflush(stdout);
					getchar();
				}
				break;
			case 1:
			case 2:
				printf ("imode> Fast Kinetic-Energy matrix Building O(n^2) [kineticMFAx()] ");
				fflush(stdout);
				kineticMFAx( mol, coord, props, size, &mass_tr, type, model, fixed, addrot ); // H-matrix (triangular)
				break;
			case 3:
				printf("imode>  ?) N,CA,C-model K-matrix with kineticMCA3x() NOT IMPLEMENTED YET\n"
						"Try \"naive method via K-matrix or V/W-arrays\" instead.\n");
				exit(1);
				break;
			}
			break;
		}
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());
		fflush(stdout);

		// Saving matrix
		if( save_matrices && !(kinetic_method == 2 && model == 0) )
		{
			timer.startTimer();
			save_matrixB(mass_tr,size,dummy_string);
			if(verb > 1)
			{
				timer.stopTimer();
				printf("imode> Saving Kinetic Energy (B) matrix in %s: %lf s\n",dummy_string,timer.getElapsedTime());
				fflush(stdout);
			}
		}
		if(save_matrices) // Storing the file saving summary
		{
			sprintf(text,"imode> Kinetic energy matrix binary file:   %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
		}

		if(verb > 2)
			show_matrix(mass_tr, size, "Kinetic Energy Matrix (H) mass_tr:");

		if(squared_matrices) // If symmetric matrices instead of packed storage...
		{
			// Allocating single precision memory
			if( !(mass_sq = (double *) malloc( sizeof( double ) * size*size ) ) )
			{
				printf("Msg(diag): I'm sorry, unable to allocate %d bytes\nForcing exit\n",sizeof( double )*size*size);
				exit(1);
			}
			for ( i = 0; i < size; i++ )
				for ( j = i; j < size; j++ ) // diagonal included
				{
					mass_sq[i+size*j] = mass_tr[i + j*(j+1)/2]; // upper triangle + diagonal
					if(i != j) // avoids writing diagonal twice...
						mass_sq[j+size*i] = mass_tr[i + j*(j+1)/2]; // lower triangle
				}
			free(mass_tr); // not needed any more
		}

		// Forces exit if "just_matrices"
		if(just_matrices)
		{
			printf("imode> Just matrices! Matrices were computed and saved!\nExiting!\n");
			exit(0);
		}
	}

	floating *eigvec = NULL;
	float *seigvec = NULL;
	floating *eigval; // DOUBLE PRECISION
	float *seigval; // SINGLE PRECISION
	if(!inevec_switch) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
	{
		// Allocating eigenvector storage memory "hess_matrix"
		// FLOATING precission eigenvectors
		printf("imode> Eigenvector matrix mem. = %.3f Mb\n",((float)sizeof(floating)*size*nevec)/1e6);

		if( !( eigvec = (floating *) malloc(size*nevec * sizeof(floating)) ) )
		{
			printf("Sorry, unable to allocate Eigenvectors memory!!! Try a lower \"nevec\".\nForcing exit!!!\n\n");
			exit(1);
		}
		// "eigvec" initialization
		for(i=0; i<size*nevec; i++)
			eigvec[i]=0.0;

		// Eigenvalue memory allocation & initialization
		if(verb > 1)
			printf("imode> Allocating eigenvalues memory.\n");
		if( !( eigval  = (floating *) malloc(size * sizeof(floating)) ) )
		{
			printf("Sorry, unable to allocate Eigenvalue memory!!!\nForcing exit!!!\n\n");
			exit(1);
		}
		for(i=0; i<size; i++)
			eigval[i]=0.0; // initialization

		// DIAGONALIZATION
		timer.startTimer();

		switch(eigensolver)
		{
		case 0:
			printf("imode> Diagonalization with LAPACK/BLAS-based XSPGVX()... ");
			fflush(stdout);
			diag_xspgvx(hess_tr, mass_tr, eigval, eigvec, size, nevec); // detects with sizeof() the floating point precision
			break;
		case 1:
			printf("imode> Diagonalization with ARPACK-based dsdrv1_AP_BP_W_mon()... ");
			fflush(stdout);
			dsdrv1_AP_BP_W_mon(hess_tr, mass_tr, eigval, eigvec, size, nevec);
			break;
		case 2:
			printf("imode> Diagonalization with ARPACK-based dsdrv1_A_B_W_mon()... ");
			fflush(stdout);
			dsdrv1_A_B_W_mon(hess_sq, mass_sq, eigval, eigvec, size, nevec);
			break;
		}

		//	printf("imode>  Diagonalizing with LAPACK's DSPGVX()...\n");
		//	fflush(stdout);
		//	diag_dspgvx(shess_tr, smass_tr, seigval, seigvec, size, nevec);
		//	printf("imode>  Diagonalizing with LAPACK's SSPGVX()...\n");
		//	fflush(stdout);
		//	diag_sspgvx(shess_tr, smass_tr, seigval, seigvec, size, nevec);
		//	printf("imode>  Diagonalizing with LAPACK's DSPGVD()...\n");
		//	fflush(stdout);
		//	diag_dspgvd(hess_tr, mass_tr, eigval, eigvec, size);
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());

		if(squared_matrices)
		{
			free(hess_sq);
			free(mass_sq);
		}
		else
		{
			free(hess_tr);
			free(mass_tr);
		}
	}
	else // If the eigenvalues/vectors were already provided by the user
	{
		// Reads ptraj file (allocating memory)
		int natoms_ptraj, nvecs, ncomps;
		printf("imode> Reading eigenvectors/values from: %s",file_inevec);
		read_ptraj(file_inevec, &eigvec, &eigval, &natoms_ptraj, &nvecs, &ncomps);
		if(nevec < nvecs)
			modes_saved = nevec; // Saving only the required number of modes (nevec)
		else // If maximum number of available modes was exceeded by required number of modes (nevec)
			nevec = nvecs; // then all read eigenvectors will be translated into Cartesian.
		printf(" (%d vectors and %d components)\n",nevec,ncomps);
	}

	// Showing the first eigenvalues
	int max=10;
	if(max > nevec)
		max = nevec;
	printf("imode> Showing the first %d eigenvalues:\n",max);
	printf("imode>\nimode> %4s %12s\n","MODE","EIGENVALUE");
	// FLOATING precision diagonalization
	for(int i=0; i<max; i++)
		printf("imode> %4d %12.5e\n",i+1, eigval[i]);
	printf("imode>\n");

// Save eigenvectors in text-format (why?...)
//
//	if(save_matrices)
//	{
//		// Saving a rectangular matrix (into file)
//		sprintf(dummy_string,"%s_evec.txt",name);
//		if(verb > 1)
//			printf("imode> Saving Eigenvector matrix in %s\n ",dummy_string);
//		sprintf(text,"imode> Eigenvector matrix text file:        %35s\n",dummy_string);
//		saved_files_len = strlen(saved_files);
//		strcpy(saved_files+saved_files_len,text);
//		save_matrix_rec(eigvec, nevec, size, dummy_string);
//	}

	// Saving & Normalizing(or not) Internal Coordinates Eigenvectors in "ptraj" format
	if(save_dh)
	{
		sprintf(file_ptraj,"%s_ic.evec",name); // ptraj name
		if(verb > 1)
			printf("imode> Saving the %d modes (Internal Coordinates): %s\n",modes_saved,file_ptraj);
		sprintf(text,"imode> ICS eigenvectors:                    %35s\n",file_ptraj);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		save_ptraj_modes(file_ptraj, size, 0, modes_saved, eigval, eigvec, norm_modes); // Saves & Normalizes
	}

	// Conversion from IC into CC is needed
	// K-matrix computation
	if( conv_method == 0 && (save_cart || save_wcart || deform_switch) && der == NULL )
		// The K-matrix (derivatives that convert internal coords. into cartesian coords)
		// is needed to compute Cartesian and Mass-Weighted Cartesian modes.
		//  Compute derivatives of cartesian coordinates w/r/t dihedrals
		switch(model)
		{
		case 0: // (first-N, CAs, last-C model)
			if(verb > 1)
				printf ("imode> Computing K-matrix with dydqMCA()\n");
			dydqMCAx(mol,coordNCAC, &der, size, fixed);
			break;
		case 3: // (N,CA,C-model)
			if(verb > 1)
				printf("imode> Computing N,CA,C-model K-matrix with dydqMCA3x()\n");
			dydqMCA3x(mol,coord, &der, size, fixed);
			break;
		case 1: // 3BB2R model
		case 2: // Full-Atom model
			if(verb > 1)
				printf ("imode> Computing K-matrix with dydqMFAx()\n");
			dydqMFAx(mol,coord, props, &der, type, model, size, fixed, addrot);
			break;
		}

	// K-matrix computation (using VW-arrays)
	if( conv_method == 1 && (save_cart || save_wcart || deform_switch) && ( V == NULL || W == NULL ) ) // then V/W-arrays needed
		switch(model)
		{
		case 0: // (first-N, CAs, last-C model)
		case 3:
			printf("imode> Computing CA-only K-matrix with vwMCA3x() \n");
			vwMCA3x(molNCAC,coordNCAC,&V,&W,&body1,size,fixed,model);
			break;
		case 1: // 3BB2R model
		case 2: // Full-Atom model
			printf("imode> Computing HA or C5 K-matrix with vwMFAx() \n");
			vwMFAx(mol,coord,props,&V,&W,&body1,type,model,size,fixed,addrot);
			break;
		}

	// Allocating Cartesian modes memory
	floating *evec, *evec_CA;
	int size_cart, size_CA;
	if( save_cart || save_wcart || deform_switch )
	{
		size_cart = num_atoms * 3;
		if( !(evec=(floating *) malloc( size_cart * sizeof(floating) * modes_saved ) ) )
		{
			printf("Error, unable to allocate Cartesian NM's memory!\n");
			exit(1);
		}
	}

	// K-matrix IC to CC conversion
	if( conv_method == 0 && (save_cart || save_wcart || deform_switch) )
	{
		if(save_cart || save_CA || save_wcart || deform_switch) // then Compute CARTESIAN MODES
		{
			// Translating Internal Coords. Eigenvectors into Cartesian Coords.
			printf("imode> Converting ICS modes into CCS modes [di2cart()]\n");
			di2cart(eigvec,der,evec,size,num_atoms,modes_saved);
		}
//		if(save_cart2)
//		{
//			// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format (3BB2R-model)
//			sprintf(file_ptraj_cart,"%s_cart.evec",name); // ptraj name
//			if(verb > 1)
//				printf("imode> Saving %d modes with \"ptraj\" format (cartesian): %s\n",modes_saved,file_ptraj_cart);
//			sprintf(text,"imode> CCS eigenvectors:                    %35s\n",file_ptraj_cart);
//			saved_files_len = strlen(saved_files);
//			strcpy(saved_files+saved_files_len,text);
//			save_ptraj_modes(file_ptraj_cart, size_cart, 0, modes_saved, eigval, evec, norm_modes);
//		}
	}

	// V/W-arrays for IC to CC conversion
	if( conv_method == 1 && (save_cart || save_wcart || deform_switch) ) // then V/W-arrays needed
	{
		if(save_cart || save_wcart || deform_switch) // Compute CARTESIAN MODES
		{
			// Translating Internal Coords. Eigenvectors into Cartesian Coords.
			if(nthreads > 0) // Parallel
			{
				printf("imode> Converting ICS modes into CCS modes in parallel (%d threads) [di2cartVW_par()]\n",nthreads);
				if(model == 0)
					di2cartVW_par(eigvec,coord,V,W,body1,evec,size,num_atoms,modes_saved,model,nthreads);
				else if( model == 3)
					di2cartVW_par(eigvec,coordNCAC,V,W,body1,evec,size,num_atoms,modes_saved,model,nthreads);
				else
					di2cartVW_par(eigvec,coord,V,W,body1,evec,size,num_atoms,modes_saved,model,nthreads);
			}
			else // Not Parallel
			{
				printf("imode> Converting ICS modes into CCS modes [di2cartVW()]\n");
				if(model == 0)
					di2cartVW(eigvec,coord,V,W,body1,evec,size,num_atoms,modes_saved,model);
				else if(model == 3)
					di2cartVW(eigvec,coordNCAC,V,W,body1,evec,size,num_atoms,modes_saved,model);
				else
					di2cartVW(eigvec,coord,V,W,body1,evec,size,num_atoms,modes_saved,model);
			}
		}

//		if( save_cart2 )
//		{
//			// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format (3BB2R-model)
//			sprintf(file_ptraj_cart,"%s_cart.evec",name); // ptraj name
//			if(verb > 1)
//				printf("imode> Saving %d modes \"ptraj\" format (cartesian): %s\n",modes_saved,file_ptraj_cart);
//			sprintf(text,"imode> CCS eigenvectors:                    %35s\n",file_ptraj_cart);
//			saved_files_len = strlen(saved_files);
//			strcpy(saved_files+saved_files_len,text);
//			save_ptraj_modes(file_ptraj_cart, size_cart, 0, modes_saved, eigval, evec, norm_modes);
//		}
		// Should we include CA eigenvectors, or just deleting it ????? (see above)

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

	if( save_cart2 )
	{
		// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format (3BB2R-model)
		sprintf(file_ptraj_cart,"%s_cart.evec",name); // ptraj name
		if(verb > 1)
			printf("imode> Saving %d modes \"ptraj\" format (cartesian): %s\n",modes_saved,file_ptraj_cart);
		sprintf(text,"imode> CCS eigenvectors:                    %35s\n",file_ptraj_cart);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		save_ptraj_modes(file_ptraj_cart, size_cart, 0, modes_saved, eigval, evec, norm_modes);
	}

	if(save_CA)
	{
		// Cartesian Eigenvector ---> CA-cartesian eigenvector
		size_CA = num_res * 3;
		if( !(evec_CA=(floating *) malloc( size_CA * sizeof(floating) * modes_saved ) ) )
		{
			printf("Error, unable to allocate CA Cartesian NM's memory!\n");
			exit(1);
		}

		// Extracting CA eigenvector components for CA coarse-graining model comparison
//			pdbIter *iter1 = new pdbIter( molr ); // <-- CHECK THIS ???
		pdbIter *iter1;
		Atom *at;
		int atom_index=0;
		int res_index=0;

		for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screen atoms
		{
			at = iter->get_atom();
			if ( !strncmp( at->getName(), " CA ", 4 ) || !strncmp( at->getName(), " P  ", 4 ) )
			{
				for(int n=0; n<modes_saved; n++)
				{
					evec_CA[n*size_CA+res_index*3]   = evec[n*size_cart+atom_index*3];
					evec_CA[n*size_CA+res_index*3+1] = evec[n*size_cart+atom_index*3+1];
					evec_CA[n*size_CA+res_index*3+2] = evec[n*size_cart+atom_index*3+2];
//					fprintf(stderr,"ires= %d  x= %f  n= %d\n",res_index,evec_CA[n*size_CA+res_index*3],n);
				}
				res_index++;
			}
			atom_index++;
		}
//		for ( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			res_index = iter->pos_fragment; // shorter
//			res = (Residue *) iter->get_fragment();
//			iter1 = new pdbIter( res );
//			for ( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() ) // screen atoms
//			{
//				if(	iter1->pos_atom == 1 ) // CA position index
//				{						   // (needed to compute "atom_index" properly!)
//					for(int n=0; n<modes_saved; n++)
//					{
//						evec_CA[n*size_CA+res_index*3]   = evec[n*size_cart+atom_index*3];
//						evec_CA[n*size_CA+res_index*3+1] = evec[n*size_cart+atom_index*3+1];
//						evec_CA[n*size_CA+res_index*3+2] = evec[n*size_cart+atom_index*3+2];
//						fprintf(stderr,"ires= %d  x= %f  n= %d\n",res_index,evec_CA[n*size_CA+res_index*3],n);
//					}
//				}
//				atom_index++;
//			}
//			delete iter1;
//		}

		// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format (CA-model)
		if(save_CA2)
		{
			sprintf(file_ptraj_cart,"%s_ca.evec",name); // ptraj name
			if(verb > 1)
				printf("imode> Saving %d modes with \"ptraj\" format: %s\n",modes_saved,file_ptraj_cart);
			sprintf(text,"imode> CCS eigenvectors (CA):               %35s\n",file_ptraj_cart);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			save_ptraj_modes(file_ptraj_cart, size_CA, 0, modes_saved, eigval, evec_CA, norm_modes);
		}
	}
	iter->clean_virtual();
	delete iter;

	// COMPUTING DEFORMABILITY
	if(deform_switch)
	{
		double *deform, *mob;
		double *bf;
		char *list;
		bf = (double *) malloc(num_atoms * sizeof(double));

//		molr->writePDB( "before.pdb" );

		mol->get_Pdbfact(bf); // ein?
//		for(int i=0; i<num_atoms; i++)
//			fprintf(stderr,"i= %4d  bf= %10f\n",i,bf[i]);

//		molr->writePDB( "after.pdb" );

		// compute pairwise distance matrix
		double *dist_matrix;
		mol->distanceMatrix( &dist_matrix );

		printf("imode> Deformability computations started...\n");
		timer.startTimer();
		compute_def (num_atoms, coord, dist_matrix, eigval, evec, nevec, 0, &deform, &mob, dc);
		timer.stopTimer();
		printf("imode> Time computing deformability %lf s\n",timer.getElapsedTime());

		// saving deformability/mobility
		list = NULL; // no res-id because of added N- and C- terminals...
		sprintf(dummy_string,"%s_defmob.txt",name);
		sprintf(text,"imode> Deformability/Mobility table:        %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		save_ascii_defmob (num_atoms, deform, mob, bf, list, dummy_string);

		//saving mobility in PDB
		norm_defmob (num_atoms, deform);
		mol->exchange_Pdbfact (deform);
		sprintf(dummy_string,"%s_def.pdb",name);
		if(verb > 1)
			printf("imode> Saving deformability in: %s\n",dummy_string);
		sprintf(text,"imode> Deformability pdb:                   %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		mol->writePDB(dummy_string);

		//saving mobilitiy in PDB
		norm_defmob(num_atoms,mob);
		mol->exchange_Pdbfact(mob);
		sprintf(dummy_string,"%s_mob.pdb",name);
		if(verb > 1)
			printf("imode> Saving mobility in: %s\n",dummy_string);
		sprintf(text,"imode> Mobility pdb:                        %35s\n",dummy_string);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		mol->writePDB(dummy_string);

		free( dist_matrix );
		free(deform);
		free(mob);
	}

	double *mass = NULL;
	if(save_wcart)
	{
		// Masses are needed for mass-weighting
		mass = ( double * ) malloc( num_atoms * sizeof( double ) );
		iter2 = new pdbIter( mol );
		Atom *atom;
		for( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
		{
			atom = ( Atom * ) iter2->get_atom();
			mass[iter2->pos_atom] = sqrt(atom->getPdbocc()); // get mass from occupancy
		}

		// Mass-weighting eigenvectors
		cart2wcart(evec,num_atoms,modes_saved,mass);

		// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format
		if(save_wcart2)
		{
			sprintf(file_ptraj_cart,"%s_wcart.evec",name); // ptraj name
			if(verb > 1)
				printf("imode> Saving %d modes \"ptraj\" format (mass-weighted cartesian): %s\n",modes_saved,file_ptraj_cart);
			sprintf(text,"imode> Mass-Weighted CCS eigenvectors:      %35s\n",file_ptraj_cart);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			save_ptraj_modes(file_ptraj_cart, size_cart, 0, modes_saved, eigval, evec, norm_modes);
		}
	}

	// COVARIANCE MATRIX COMPUTATION
	// See: eq.28 from: Noguti & Go. 1983. Journal of the Physical Society of Japan. pp.3283-88
	//      eq.13 from: Nishikawa & Go. 1987. Proteins: Structure, Function and Genetics. 2:308-29.
	//      eqs. 1 and 9 from: Tomimoto, Kitao and Go. NMA of furanose ring in DAS. 1996.

	// Computing & Saving the NMA predicted Covariance matrix (C)
	double *covar = NULL;
	if(save_covar)
	{
		// Allocating/Initializing Covariance matrix (C)
		printf("pcatool> Allocating the Covariance matrix (C) bytes=%d\n",sizeof(double) * size_cart*(size_cart+1)/2);
		if( !(covar = (double *) malloc(sizeof(double) * size_cart*(size_cart+1)/2)) )
		{
			printf("pcatool> Sorry, unable to allocate Covariance matrix (C): %d bytes\nExiting!\n",sizeof(double) * size_cart*(size_cart+1)/2);
			exit(1);
		}
		for(int i=0; i< size_cart*(size_cart+1)/2; i++)
			covar[i] = 0.0; // initialization

		// Computing Covariance matrix (C)
		timer.startTimer(); // timer
		printf("imode> Computing the Cartesian Covariance matrix at %.1fK ",Ta);
		fflush(stdout);
		for(int i=0; i<size_cart; i++) // iter rows. (coordinate)
			for(int j=i; j<size_cart; j++) // iter cols. (coordinates, including diagonal)
			{
				// Sum up each mode contribution
				for(int k=0; k<nevec; k++) // iter computed modes
					covar[i + j*(j+1)/2] += (evec[k*size_cart+i] * evec[k*size_cart+j])/eigval[k];
				// Boltzmann Temperature factor
				covar[i + j*(j+1)/2] *= KbT;
			}
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());

		// Saving a triangular packed Covariance matrix
		timer.startTimer(); // timer
		if(save_matrices_text) // Text format
		{
			sprintf(dummy_string,"%s_covar.txt",name);
			save_matrix(covar, size_cart, dummy_string, 2, true); // save indices and symmetric square covariance matrix
			printf("imode> Saving covariance matrix (C=%d^2) (text): %s ",size_cart,dummy_string);
			sprintf(text,"imode> Covariance matrix (text):            %35s\n",dummy_string);
		}
		else // Binary format
		{
			sprintf(dummy_string,"%s_covar.bin",name); // ptraj name
			printf("imode> Saving covariance matrix (C=%d^2) (binary): %s ",size_cart,dummy_string);
			save_matrixB(covar, size_cart, dummy_string);
			sprintf(text,"imode> Covariance matrix (binary):          %35s\n",dummy_string);
		}
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());
		fflush(stdout);
		free(covar);
	}

	// See: Ichiye T, Karplus M. Collective Motions in Proteins: a Covariance Analysis of Atomic Fluctuations in Molecular-Dynamics and Normal Mode Simulations. Proteins 1991 11; 205-217 (Eq. 2)
	if(save_dcovar) // Save CA-based Distance-Covariance Matrix
	{
		// Allocating/Initializing Covariance matrix (C)
		printf("pcatool> Allocating the Distance-Covariance matrix (Cd) bytes=%d\n",sizeof(double) * num_res*(num_res+1)/2);
		if( !(covar = (double *) malloc(sizeof(double) * num_res*(num_res+1)/2)) )
		{
			printf("pcatool> Sorry, unable to allocate Covariance matrix (C): %d bytes\nExiting!\n",sizeof(double) * num_res*(num_res+1)/2 );
			exit(1);
		}
		for(int i=0; i < num_res*(num_res+1)/2; i++)
			covar[i] = 0.0; // initialization

		// Computing Covariance matrix (C)
		timer.startTimer(); // timer
		printf("imode> Computing the Distance-Covariance matrix... ");
		fflush(stdout);

		double *dr;
		if( !(dr = (double *) malloc(sizeof(double) * num_res) ) )
		{
			printf("pcatool> Sorry, unable to allocate dr memory: %d bytes\nExiting!\n",sizeof(double) * num_res );
			exit(1);
		}
		for(int i=0; i<num_res; i++) // iter rows. (coordinate)
		{
			dr[i] = 0.0;
			for(int k=0; k<nevec; k++) // iter computed modes
				dr[i] += dotp(evec_CA+k*size_CA+3*i, evec_CA+k*size_CA+3*i) / eigval[k];
		}

		double drij;
		for(int i=0; i<num_res; i++) // iter rows. (coordinate)
			for(int j=i; j<num_res; j++) // iter cols. (coordinates, including diagonal)
			{
				// Sum up each mode contribution
				drij = 0.0;
				for(int k=0; k<nevec; k++) // iter computed modes
					drij += dotp(evec_CA+k*size_CA+3*i, evec_CA+k*size_CA+3*j) / eigval[k];
				covar[i + j*(j+1)/2] = drij / sqrt(dr[i]*dr[j]);
			}
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());
		free(dr);

		// Saving a triangular packed Covariance matrix
		timer.startTimer(); // timer
		if(save_matrices_text) // Text format
		{
			sprintf(dummy_string,"%s_dcovar.txt",name);
			printf("imode> Saving covariance matrix (C=%d^2) (text): %s ",num_res,dummy_string);
			save_matrix(covar, num_res, dummy_string, 2, true); // save indices and symmetric square covariance matrix
			sprintf(text,"imode> CA-Distance-Cov. matrix (text):      %35s\n",dummy_string);
		}
		else // Binary format
		{
			sprintf(dummy_string,"%s_dcovar.bin",name); // ptraj name
			printf("imode> Saving covariance matrix (C=%d^2) (binary): %s ",num_res,dummy_string);
			save_matrixB(covar, num_res, dummy_string);
			sprintf(text,"imode> CA-Distance-Cov. matrix (binary):    %35s\n",dummy_string);
		}
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		timer.stopTimer();
		printf("%lf s\n",timer.getElapsedTime());
		fflush(stdout);
		free(covar);
	}

	// OUTPUT DIFFERENT MODEL
	Macromolecule *molf;
	if(modelf != -1 )
	{
		// Setting IC-Coarse-Graining model
		switch(modelf)
		{
		case 0: // CA-only + (NH and CO)-terminal model
			printf( "imode> Output Coarse-Graining model: CA-only\n");

			// N,CA,C selection
			molNCAC = molini->select( ncac2 );
			num_atomsNCAC = molNCAC->get_num_atoms();

			// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
			// Warning, atoms are not copied, they're just pointers to the original atoms.
			// setmass = true --> adds masses to Occupancy and Bfactor, otherwise not left unchanged.
			// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
			//                   (unit mass is divided between NH or CO and their CAs)
			// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
			//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
			molf = cg_CA( molNCAC, true, nomass_switch ); // masses will be computed
			break;

		case 3: // N,CA,C-model
			printf( "imode> Output Coarse-Graining model: N,CA,C-model\n");
			// N,CA,C selection
			molf = molini->select( ncac2 );
			mass_NCAC( molf, nomass_switch );
			break;
		case 1: // 3BB2R model
			printf( "imode> Output Coarse-Graining model: 3BB2R\n");
			molf = new Macromolecule( molini );
			cg_3BBR2(molf, nomass_switch, nomass_switch);
			break;

		case 2: // Full-Atom
			printf( "imode> Output Coarse-Graining model: Full-Atom (no coarse-graining)\n");
			molf = molini;
			mass_FA(molf); // sets masses
			break;
		}
		num_atoms = molf->get_num_atoms();
		size_cart = num_atoms * 3;
		printf( "imode> Output model (pseudo)atoms: %d\n", num_atoms );

		// Saving "f" output model
		if(savemodel_switch)
		{
			sprintf(dummy_string,"%s_modelf.pdb",name);
			sprintf(text,"imode> Output model PDB:                    %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			molf->writePDB( dummy_string );
		}

		tri *propsf;
		switch(modelf)
		{
		case 0: // (first-N, CAs, last-C model)
		case 3: // (N,CA,C-model)
			properCA(molf,&propsf,&unat);
			break;
		case 1: // 3BB2R model
		case 2: // Full-Atom model
			properMFA(molf,&propsf,&unat,typef,modelf);
			break;
		}

		// Theoric number of DoFs= 3*T + 3*R -6 + Dihedral
		// (Before fixing, i.e. the fix-file format DoFs...)
		int sizef = 0;
		int old_size,seg_atoms=0;
		sizef = 0;
		old_size = sizef; // temp (Dihedral ICs)
		pdbIter *iterf = new pdbIter( molf ); // iter to screen fragments (residues)
		for( iterf->pos_fragment = 0; !iterf->gend_fragment(); iterf->next_fragment() ) // screen residues
			sizef += propsf[iterf->pos_fragment].nan; // Each residue may have different number of dihedral-angles!
		printf( "imode> Output model Dihedral angles: %d\n", sizef );
		old_size = sizef;
		for( iterf->pos_segment = 0; !iterf->gend_segment(); iterf->next_segment() ) // screen segments
		{
			seg_atoms = (iterf->get_segment())->num_atoms();
			if(seg_atoms > 1)
				sizef += 6;
			else
				sizef += 3;
		}
		iterf->clean_virtual();
		delete iterf;
		printf( "imode> Output model Rotational/Translational ICs (Non-Eckart): %d\n", sizef-old_size );
		sizef -= 6; // Eckart conditions
		printf( "imode> Output model predicted number of ICs (Eckart): %d\n", sizef );

		bool *fixedf=NULL;
		if(fixed != NULL)
		{
			fixedf = (bool *) malloc( sizeof( bool ) * sizef );

			// Change fix-arrays (internal coordinates) from a "i" model into a "f" model.
			// WARNING: "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
			//          (no memory allocation is performed)
			// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
			// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
			// "i"-model should be allways "finer" than "f"-model
			if(model > modelf) // i-model > f-model  --> normal
				changefixIC(mol, props, propsf, fixed, fixedf, model, type, modelf, typef, false);
			else if (model == modelf) // i-model == f-model  --> normal
			{
				if(type > typef) // i-type >= f-type  --> normal
					changefixIC(mol, props, propsf, fixed, fixedf, model, type, modelf, typef, false);
				else if(type < typef)  // i-model == f-model && i-type < f-type  --> inverse
					changefixIC(molf, propsf, props, fixedf, fixed, modelf, typef, model, type, true);
				else
				{
					free(fixedf);
					fixedf = fixed;
				}
			}
			else // i-model < f-model  --> inverse
				changefixIC(molf, propsf, props, fixedf, fixed, modelf, typef, model, type, true);

			old_size=0; // counter
			for(i=0; i<sizef; i++)
				if(!fixedf[i]) // if "f" DoF is fixed
					old_size++;
			printf( "imode> Output model fixed coordinates: %d \n", old_size );
			sizef -= old_size;
			printf( "imode> Output model mobile coordinates (sizef): %d \n", sizef );
		}

		if(save_fixfile)
		{
			sprintf(dummy_string,"%sf.fix",name);
			write_fixDH(dummy_string,molf,propsf,fixedf,typef,modelf); // Output fixation mask
			if(verb > 1)
				printf("imode> Output Fixation mask written: %s\n",dummy_string);
			sprintf(text,"imode> Output Fixation mask file:                 %29s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
		}

		// Creates two auxiliar arrays with segment properties (needed due to fixing):
		//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
		//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
		// (Allocates memory itself, if it's needed)
		bool *addrotf=NULL; // Should be added ROTATIONs? (with fixing)
		int *effsegf=NULL; // Effective segment indices
		sizef = seg_props(molf, propsf, fixedf, modelf, typef, &addrotf, &effsegf);
		printf( "imode> Number of ICs predicted by seg_props(): %d (final)\n", sizef);

		// "f" FLOATING precision eigenvectors
		floating *evecf;
		if(verb > 1)
			printf("imode> Eigenvector matrix memory-size (floating)= %.3f Mb\n",((float)sizeof(floating)*sizef*modes_saved)/1e6);
		if( !( evecf = (floating *) malloc(sizef*modes_saved * sizeof(floating)) ) )
		{
			printf("Sorry, unable to allocate Eigenvectors memory!!! Try a lower \"nevec\".\nForcing exit!!!\n\n");
			exit(1);
		}
		// "eigvec" initialization
		for(i=0; i<sizef*modes_saved; i++)
			evecf[i]=0.0;

		// Change normal modes (internal coordinates) from a "i" model into a "f" model.
		// WARNING: "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
		// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
		// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
		// "i"-model should be allways "finer" than "f"-model
		if(model > modelf) // i-model > f-model  --> normal
			changemodelIC(mol, props, propsf, eigvec, evecf, fixed, fixedf, nevec, size, sizef, model, type, modelf, typef, false);
		else if (model == modelf) // i-model == f-model  --> normal
		{
			if(type > typef) // i-type > f-type  --> normal
				changemodelIC(mol, props, propsf, eigvec, evecf, fixed, fixedf, nevec, size, sizef, model, type, modelf, typef, false);
			else if(type < typef) // i-type < f-type  --> inverse
				changemodelIC(molf, propsf, props, evecf, eigvec, fixedf, fixed, nevec, sizef, size, modelf, typef, model, type, true);
			else // same models, copying eigenvector reference
			{
				free(evecf);
				evecf = eigvec;
			}
		}
		else // i-model < f-model  --> inverse
			changemodelIC(molf, propsf, props, evecf, eigvec, fixedf, fixed, nevec, sizef, size, modelf, typef, model, type, true);

		sprintf(file_ptraj,"%s_icf.evec",name); // ptraj name
		if(verb > 1)
			printf("imode>\nimode> Saving the %d modes (Internal Coordinates): %s\n",modes_saved,file_ptraj);
		sprintf(text,"imode> Output model ICS eigenvector ptraj:  %35s\n",file_ptraj);
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		save_ptraj_modes(file_ptraj, sizef, 0, nevec, eigval, evecf, norm_modes); // Saves & Normalizes

		free(coord);
		molf->coordMatrix( &coord );
		if(model==0)
			free(coordNCAC);
		if(modelf==0)
			molNCAC->coordMatrix( &coordNCAC );
		iter2 = new pdbIter( molf ); // Iterator to screen atoms

		// Shift coordinates to their Center of Mass (Computing the CoM)
		// The structure should be centered on its CoM in order
		// to compute Kinetic-energy and Hessian Matrices.
		double mtot = 0.0;
		double mta = 0.0;
		double rd[3];
		rd[0] = rd[1] = rd[2] = 0.0;
		for ( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )  // screens all-atoms
		{
			mta = ( iter2->get_atom() )->getPdbocc(); // Load mass from occupancies...
			//		if( mta != 0.0 ) // some atoms are virtual, saving time...
			//		{
			mtot += mta;
			// Sum(mass*coord) before putting the CoM at 0
			rd[0] += mta * coord[iter2->pos_atom * 3];
			rd[1] += mta * coord[iter2->pos_atom * 3 + 1];
			rd[2] += mta * coord[iter2->pos_atom * 3 + 2];
			//		}
		}
		delete(iter2);
		rd[0] /= mtot;
		rd[1] /= mtot;
		rd[2] /= mtot;
		if(verb > 1)
			printf( "imode> Mass %8.8f Center %8.8f %8.8f %8.8f --> Shifting to origin\n", mtot, rd[0], rd[1], rd[2] );
		// Shifting CoM into origin
		for(int k = 0; k < num_atoms; k++)
		{
			coord[k * 3] -= rd[0];
			coord[k * 3 + 1] -= rd[1];
			coord[k * 3 + 2] -= rd[2];
		}

		if(modelf==0)
			for(int k = 0; k < num_atomsNCAC; k++)
			{
				coordNCAC[k * 3] -= rd[0];
				coordNCAC[k * 3 + 1] -= rd[1];
				coordNCAC[k * 3 + 2] -= rd[2];
			}

		if( conv_method == 0 && (save_cart || save_wcart ) ) // then K-matrix used
		{
			free(der);
			der = NULL;
			switch(modelf)
			{
			case 0:
				if(verb > 1)
					printf("imode> Computing output model CA-only K-matrix with dydqMCAx()\n");
				dydqMCAx(molf,coordNCAC, &der, sizef, fixedf);
				break;
			case 3:
				// (N,CA,C-model)
				if(verb > 1)
					printf("imode> Computing output model CA-only K-matrix with dydqMCA3x()\n");
				dydqMCA3x(molf,coord, &der, sizef, fixedf);
				break;
			case 1:
			case 2:
				if(verb > 1)
					printf("imode> Computing output model Full-Atom/3BB2R K-matrix with dydqMFAx()\n");
				dydqMFAx(molf,coord, propsf, &der, typef, modelf, sizef, fixedf, addrotf);
				break;
			}
		}

		V = NULL; // previously deleted
		W = NULL;
		body1 = NULL;
		if( conv_method == 1 && (save_cart || save_wcart ) && ( V == NULL || W == NULL ) ) // then V/W-arrays needed
			switch(modelf)
			{
			case 0: // (first-N, CAs, last-C model)
			case 3: // (N,CA,C-model)
				if(verb > -1)
					printf("imode> Computing output model CA-only K-matrix with vwMCA3x() \n");
				vwMCA3x(molNCAC,coordNCAC,&V,&W,&body1,sizef,fixedf,modelf); // <-- Check this!
				break;
			case 1: // 3BB2R model
			case 2: // Full-Atom model
				if(verb > -1)
					printf("imode> Computing output model Full-Atom/3BB2R K-matrix with vwMFAx() \n");
				vwMFAx(molf,coord,propsf,&V,&W,&body1,typef,modelf,sizef,fixedf,addrotf);
				break;
			}

		// "f" cartesian FLOATING precision eigenvectors
		floating *evecfc;
		printf("imode> Output model eigenvector matrix mem. (floating)= %.3f Mb\n",((float)sizeof(floating)*3*num_atoms*modes_saved)/1e6);
		if( !( evecfc = (floating *) malloc(3*num_atoms*modes_saved * sizeof(floating)) ) )
		{
			printf("Sorry, unable to allocate Eigenvectors memory!!! Try a lower \"nevec\".\nForcing exit!!!\n\n");
			exit(1);
		}
		// "eigvec" initialization
		for(int i=0; i<3*num_atoms*modes_saved; i++)
			evecfc[i]=0.0;

		// K-matrix is used
		if( conv_method == 0 && (save_cart || save_wcart) )
		{
			// Translating Internal Coords. Eigenvectors into Cartesian Coords.
			printf("imode> Converting ICS modes into CCS modes [di2cart()]\n");
			di2cart(evecf,der,evecfc,sizef,num_atoms,modes_saved);
		}

		// V/W-arrays are used
		if( conv_method == 1 && (save_cart || save_wcart) )
		{
			// Translating Internal Coords. Eigenvectors into Cartesian Coords.
			if(nthreads > 0) // Parallel
			{
				printf("imode> Converting ICS modes into CCS modes [di2cartVW_par()]\n");
				di2cartVW_par(evecf,coord,V,W,body1,evecfc,sizef,num_atoms,modes_saved,model,nthreads);
			}
			else // Not Parallel
			{
				printf("imode> Converting ICS modes into CCS modes [di2cartVW()]\n");
				di2cartVW(evecf,coord,V,W,body1,evecfc,sizef,num_atoms,modes_saved,model);
			}

//			// Translating Internal Coords. Eigenvectors into Cartesian Coords.
//			printf("imode> Converting ICS modes into CCS modes [di2cartVW()]\n");
//			di2cartVW(evecf,coord,V,W,body1,evecfc,sizef,num_atoms,modes_saved,model);
			// Lets Free V/W arrays (not used anymore)
			for(int i=0; i<sizef; i++)
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

		// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format (3BB2R-model)
		if(save_cart2)
		{
			sprintf(file_ptraj_cart,"%s_cartf.evec",name); // ptraj name
			if(verb > 1)
				printf("imode> Saving %d modes \"ptraj\" format (cartesian): %s\n",modes_saved,file_ptraj_cart);
			sprintf(text,"imode> Output model CCS eigenvector ptraj:  %35s\n",file_ptraj_cart);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			save_ptraj_modes(file_ptraj_cart, num_atoms*3, 0, modes_saved, eigval, evecfc, norm_modes);
		}

		free(mass); // WATCH OUT THIS
		if(save_wcart)
		{
			// Masses are needed for mass-weighting
			mass = ( double * ) malloc( num_atoms * sizeof( double ) );
			iter2 = new pdbIter( molf );
			Atom *atom;
			for( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
			{
				atom = ( Atom * ) iter2->get_atom();
				mass[iter2->pos_atom] = sqrt(atom->getPdbocc()); // get mass from occupancy
			}

			// Mass-weighting eigenvectors
			cart2wcart(evecfc,num_atoms,modes_saved,mass);

			if(save_wcart2)
			{
				// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format
				sprintf(file_ptraj_cart,"%s_wcartf.evec",name); // ptraj name
				if(verb > 1)
					printf("imode> Saving %d modes \"ptraj\" format (mass-weighted cartesian): %s\n",modes_saved,file_ptraj_cart);
				sprintf(text,"imode> Output model MW-CCS ptraj:           %35s\n",file_ptraj_cart);
				saved_files_len = strlen(saved_files);
				strcpy(saved_files+saved_files_len,text);
				save_ptraj_modes(file_ptraj_cart, num_atoms*3, 0, modes_saved, eigval, evecfc, norm_modes);
			}
		}

		// Computing & Saving the NMA predicted Covariance matrix (C)
		if(save_covarf)
		{
			// Allocating/Initializing Covariance matrix (C)
			printf("pcatool> Allocating the Covariance matrix (C) bytes=%d\n",sizeof(double) * size_cart*(size_cart+1)/2);
			if( !(covar = (double *) malloc(sizeof(double) * size_cart*(size_cart+1)/2)) )
			{
				printf("pcatool> Sorry, unable to allocate Covariance matrix (C): %d bytes\nExiting!\n",sizeof(double) * size_cart*(size_cart+1)/2);
				exit(1);
			}
			for(int i=0; i< size_cart*(size_cart+1)/2; i++)
				covar[i] = 0.0; // initialization

			// Computing Covariance matrix (C)
			timer.startTimer(); // timer
			printf("imode> Computing final model Cartesian Covariance matrix (%.1fK) ",Ta);
			fflush(stdout);
			for(int i=0; i<size_cart; i++) // iter rows. (coordinate)
				for(int j=i; j<size_cart; j++) // iter cols. (coordinates, including diagonal)
				{
					// Sum up each mode contribution
					for(int k=0; k<nevec; k++) // iter computed modes
						covar[i + j*(j+1)/2] += (evecfc[k*size_cart+i] * evecfc[k*size_cart+j])/eigval[k];
					// Botlzmann Temperature factor
					covar[i + j*(j+1)/2] *= KbT;
				}
			timer.stopTimer();
			printf("%lf s\n",timer.getElapsedTime());

			// Saving a triangular packed Covariance matrix
			timer.startTimer(); // timer
			if(save_matrices_text) // Text format
			{
				sprintf(dummy_string,"%s_covarf.txt",name);
				save_matrix(covar, size_cart, dummy_string, 2, true); // save indices and symmetric square covariance matrix
				printf("imode> Saving final covariance matrix (C=%d^2) (text): %s ",size_cart,dummy_string);
				sprintf(text,"imode> Final model Cov. matrix (text):      %35s\n",dummy_string);
			}
			else // Binary format
			{
				sprintf(dummy_string,"%s_covarf.bin",name); // ptraj name
				printf("imode> Saving final covariance matrix (C=%d^2) (binary): %s ",size_cart,dummy_string);
				save_matrixB(covar, size_cart, dummy_string);
				sprintf(text,"imode> Final model Cov. matrix (binary):    %35s\n",dummy_string);
			}
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			timer.stopTimer();
			printf("%lf s\n",timer.getElapsedTime());
			free(covar);
		}
	}

	printf("imode>\nimode> SAVED FILES:\n%s",saved_files);
	printf("imode>\nimode> Bye!\n");
	exit(0);
}
/*==============================================================================================*/

/*==============================================================================================*/
void parseOptions(int argc, char** argv)
{
        std::string temp;
        CmdLine cmd("imode","iMODE: Internal coordinates normal MODE analysis tool.", version );

   try {
        // Define required arguments not labeled
        UnlabeledValueArg<std::string> Pdb("pdb","PDB input file","default","pdb");
        cmd.add( Pdb );

        // Define labeled arguments

        // DEVELOPER INPUT
		ValueArg<double> DC("","dc", "Characteristic Distance of the Hardy's Quadric Interpolation used in Deformability computations (default=15). It should be > 0.",false,15,"double");
		cmd.add( DC );
		ValueArg<std::string> InEvec("","inevec", "Input IC Eigenvectors/values file (.evec). "
				"This disables Hessian and Kinetic energy matrices calculation and diagonalization.",false,"InEvec","string");
        cmd.add( InEvec );
		ValueArg<float> InterMolec("", "inter_molec","Sets the inter-molecular force constant factor (default=disabled)",false,1.0,"float");
		cmd.add( InterMolec );
        ValueArg<int> Eigensolver("","eigensolver", "Eigensolver (pay attention to this option specially for large systems):\n"
        										    "0= LAPACK/BLAS, fastest if more than 10% modes requested [DSPGVX],\n"
        		                                    "1= ARPACK, fastest if less than 5% modes requested [dsdrv1_AP_BP_W_mon] (default),\n"
        		                                    "2= ARPACK-square, [dsdrv1_A_B_W_mon] (experimental).",false,1,"int");
		cmd.add( Eigensolver );
        ValueArg<int> NThreads("","nthreads", "Number of threads for parallel processing (experimental)",false,1,"int");
		cmd.add( NThreads );
        ValueArg<std::string> FixIC("","fixIC", "Plain-text file defining the fixed Internal Coordinates. "
	    		"Each line will contain the index (0,1,...) of the ICs to be removed "
	    		"(DEVELOPER's).",false,"fixstring","string");
        cmd.add( FixIC );
        ValueArg<int> Debug("","debug", "Debug code <int> (default=disabled).",false,0,"int");
        cmd.add( Debug );
        ValueArg<int> Convert("","convert", "Conversion method from ICS to CCS: 0=K-matrix, 1=VW-arrays (DEVELOPER's) (default=1).",false,1,"int");
        cmd.add( Convert );
        ValueArg<int> Kinetic("","kinetic", "Kinetic energy matrix building method (DEVELOPER's) (default=2).",false,2,"int");
        cmd.add( Kinetic );
        ValueArg<int> Hessian("","hessian", "Hessian matrix building method (DEVELOPER's) (default=2).",false,2,"int");
        cmd.add( Hessian );
        ValueArg<float> SwapMatrixMem("","swapmatrix_mem", "Amount of RAM memory (in GB) to be used during Hessian matrix swapping (default=1GB) (DEVELOPER's).",false,1,"float");
		cmd.add( SwapMatrixMem );
		SwitchArg JustMatrices("","just_matrices", "Just computes and saves both matrices, then exit... (default=disabled) (DEVELOPER's).", true);
		cmd.add( JustMatrices );
		SwitchArg SaveMatrix("","save_matrices", "Saves both Hessian and Kinetic energy matrices in binary packed storage format (default=disabled) (DEVELOPER's).", true);
		cmd.add( SaveMatrix );
        ValueArg<float> FixRand2("R","fixRand2", "Randomly fixed ratio of Internal Coordinates (default=disabled). Example: 0.7 = 70% of the ICs will be randomly fixed (DEVELOPER's).",false,0.5,"float");
		cmd.add( FixRand2 );
        // END DEVELOPER INPUT

		// Less important parameters
		ValueArg<int> Verb("","verb", "Verbose level (0=low, 1=medium, 2=high) (default=0).",false,0,"int");
        cmd.add( Verb );

		SwitchArg DelHeteros("","delete_heteros", "Delete Hetero-atoms, including waters (default=disabled).", true);
		cmd.add( DelHeteros );
		SwitchArg KeepWaters("","keep_waters", "Disables Water molecules deletion (default=disabled).", true);
		cmd.add( KeepWaters );
		SwitchArg KeepHydrogens("","keep_hydrogens", "Disables Hydrogen atoms deletion (default=disabled).", true);
		cmd.add( KeepHydrogens );

        ValueArg<unsigned int> Seed("","seed", "Pre-define the random number generator SEED (Mersenne Twister) (default=random-seed from /dev/urandom)",false,386,"unsigned int");
		cmd.add( Seed );
		ValueArg<double> Temp("T","temperature", "Temperature [K] for covariance matrix computation (default=300).",false,300,"double");
		cmd.add( Temp );
		SwitchArg SaveCovarOut("","save_covar_out", "Computes and Saves the predicted covariance matrix for the output model at selected Temperature in binary packed storage format as <basename_covarf.bin>. "
				"If --save_wcart selected, then mass-weighted covariance matrix will be computed instead (default=disabled).", true);
		cmd.add( SaveCovarOut );
		SwitchArg ChiOut("","chi_out", "Considers first CHI dihedral angle in output modes (default=disabled).", true);
		cmd.add( ChiOut );
        ValueArg<int> ModelOut("","model_out", "Output Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom (default=disabled).",false,-1,"int");
        cmd.add( ModelOut );

		ValueArg<std::string> Funcfile("","func", "ASCII file defining the force constant functions to be applied "
				"according to Topology and/or Secondary Structure. "
				"The 5 cols. format is: <SS> <t> <k> <x0> <pow>\n"
				"Where <SS> is the two character pairwise interaction identifier, <t> is the topology, and "
				"<k>,<x0>,<pow> are the corresponding sigmoid function parameters. If --ss "
				"is not specified, the XX pairwise interaction identifier must be introduced. This way, only topologies will be considered. "
				"If <t> is \"-1\", any previously not-matched topology will be considered.",false,"","string");
		cmd.add( Funcfile );
		SwitchArg Norm("", "norm","Enables (norm=1) eigenvector normalization. "
				"Note this does not change vector direction (default=disabled).", true);
		cmd.add( Norm );
		SwitchArg NoTors("","notors", "Disables extra torsional potential (default=disabled).", true);
		cmd.add( NoTors );
		SwitchArg NoMass("","nomass", "Disables mass weighting (default=disabled).", true);
		cmd.add( NoMass );
		SwitchArg NoModel("", "nomodel","Disables PDB model building. "
				"Warning: introduced PDB model must match the CG selected with the -m option (default=disabled).", false);
		cmd.add( NoModel );

		ValueArg<float> k2_Cutoff("","k2_c","Non-bonding distance cutoff applied to --func option (default=10A).", false, 10,"float");
		cmd.add( k2_Cutoff );
		ValueArg<float> k1_Cte("", "k1_k","Tirion's method stiffness constant (default=1.0).",false,1.0,"float");
		cmd.add( k1_Cte );
		ValueArg<float> k1_Cutoff("","k1_c","Tirion's method distance cutoff (default=10A).", false, 10,"float");
		cmd.add( k1_Cutoff );
		ValueArg<float> k0_Power("","k0_p", "Sigmoid function power term (default=6).",false,6.0,"float");
		cmd.add( k0_Power);
		ValueArg<float> k0_X0("", "k0_x0","Sigmoid function inflexion point (default=3.8A).",false,3.8,"float");
		cmd.add( k0_X0 );
		ValueArg<float> k0_Cte("", "k0_k","Sigmoid function stiffness constant (default=1.0).",false,1.0,"float");
		cmd.add( k0_Cte );
		ValueArg<float> k0_Cutoff("","k0_c","Sigmoid function distance cutoff (default=10A).", false, 10,"float");
		cmd.add( k0_Cutoff );

		SwitchArg SaveCovarText("","save_covar_text", "Enables plain text output for the covariance matrices (default=disabled).", true);
		cmd.add( SaveCovarText );
		SwitchArg SaveDCovar("","save_dcovar", "Saves the CA-based distance-covariance matrix at selected Temperature ", true);
		cmd.add( SaveDCovar );
		SwitchArg SaveCovar("","save_covar", "Saves the predicted covariance matrix at selected Temperature "
				"in binary packed storage format as <basename_covar.bin>. If --save_wcart selected, "
				"then mass-weighted covariance matrix will be computed instead (default=disabled).", true);
		cmd.add( SaveCovar );
		SwitchArg SSOut("", "save_SSfile","Save secondary structure file as <basename.ss> "
				"(to be used with -S or -P=2 options) (default=disabled).", false);
		cmd.add( SSOut );
		SwitchArg KOut("", "save_Kfile","Save atom-pairwise force constants file as <basename_Kfile.dat> (to be used with -K option) (default=disabled)", false);
		cmd.add( KOut );
		SwitchArg SaveCA("", "save_ca","Save CA-based Cartesian modes as <basename_ca.evec> (default=disabled)", false);
		cmd.add( SaveCA );
		SwitchArg SaveWCart("", "save_wcart","Save Mass-weighted Cartesian modes as <basename_wcart.evec> (default=disabled)", false);
		cmd.add( SaveWCart );
		SwitchArg SaveCart("", "save_cart","Save Cartesian modes as <basename_cart.evec> (default=disabled)", false);
		cmd.add( SaveCart );
		SwitchArg SaveFixFile("", "save_fixfile","Save fixation file as <basename.fix> (to be used with -r or -S options; otherwise a fully mobile file will be generated) (default=disabled)", false);
		cmd.add( SaveFixFile );
		ValueArg<std::string> SSfile("","ss", "Secondary Structure ASCII file with 2 cols.: <n> <char>\n"
				"Where <n> is the corresponding residue index (0,1,...), and <char> is the "
				"single character SS identifier. "
				"By default SS will be computed internally (H=helix, E=strand, C=coil).",false,"","string");
		cmd.add( SSfile );

		// More important parameters (one letter)
		ValueArg<std::string> FixSS("S","fixSS", "All dihedral coordinates with a given secondary structure (SS) "
				"will be removed (see --ss). Ex: \"HE\" will fix the dihedrals corresponding to alpha-helices "
				"and beta-sheets.",false,"","string");
		cmd.add( FixSS );
		SwitchArg Chi("x","chi", "Considers first CHI dihedral angle (default=disabled).", true);
		cmd.add( Chi );
        ValueArg<float> Nevs("n","nevs", "Used modes range, either number [1,N] <integer>, or ratio [0,1) <float> (default=20).",false,20,"int/float");
        cmd.add( Nevs );
		ValueArg<std::string> KFile("K","Kfile", "Force constants ASCII file with 3 cols.: <i-atom> <j-atom> <K>\n"
				"Where <i/j-atom> are the corresponding atomic indices (1,2,...)\n"
			    "A demo file can be generated using --save_Kfile option.",false,"Kfile","string");
        cmd.add( KFile );
        ValueArg<int> Contact("P","potential", "Pairwise interaction potential: (default=0)\n"
        		"  0= Sigmoid function (= k/(1+(x/x0)^p), if x < c, else k=0).\n"
        		"  1= Tirion's cutoff (= k, if x < c, else k=0).\n"
        		"  2= Hinsen's function.\n"
        		"  3= Topology & Secondary Structure (--func is mandatory).\n"
        		"  4= edNMA formalism (CA-model only).\n"
        		"  By default an extra torsional potential will be added.",false,0,"int");
        cmd.add( Contact );
	    ValueArg<std::string> FixFile("f","fixFile", "ASCII file defining the ICs to be fixed with the format:\n"
	    		"Protein:     \"n phi chi psi\"\n"
	    		"NAcid:       \"n alpha beta gamma chi epsilon zeta\"\n"
	    		"Inter-chain: \"n 6D\"\n"
	    		"Where \"n\" is the residue index (0,1,..) and the coordinate name (phi, psi, etc...) "
	    		"can be set to 0(fixed) or 1(mobile). "
	    		"Each one of the 6 inter-chain variables should be specified on separate lines in the "
	    		"following order: x,y,z,Rx,Ry,Rz. "
	    		"Note \"n\" is just the sequential residue index (starting with 0) and NOT the PDB's residue index.\n"
	    		"A demo file can be generated using the --save_fixfile option.",false,"fixstring","string");
        cmd.add( FixFile );
        ValueArg<float> FixRand("r","fixRand", "Randomly fixed ratio of Dihedral Coordinates (default=disabled). Example: 0.7 = 70% of dihedrals will be randomly removed. Rotational/translational coords. always mobile.",false,0.5,"float");
		cmd.add( FixRand );
        SwitchArg Def("d","deform", "Turn on deformability calculations. CAUTION, only suitable for CA-model. (default=disabled).", true);
        cmd.add( Def );
        ValueArg<std::string> Name("o","name", "Output files basename (default=imode).",false,"imode","string");
        cmd.add( Name );
		ValueArg<int> Model("m","model", "Coarse-Grained model: 0=CA, 1=C5, 2=Heavy-Atom (default=2).",false,2,"int");
        cmd.add( Model );

		// Parse the command line.
        cmd.parse(argc,argv);

		// Getting the command line arguments.
        strcpy(file_pdb,((temp=Pdb.getValue()).c_str())); // Gets PDB file name
        strcpy(name,((temp=Name.getValue()).c_str())); // Gets Basename
        strcpy(Kfile,((temp=KFile.getValue()).c_str())); // Gets Force constants file-name
		Kout_switch = KOut.isSet();
		SSout_switch = SSOut.isSet();

		// Contacting method parameters
        power = k0_Power.getValue();
		cte_k0 = k0_Cte.getValue();
		x0 = k0_X0.getValue();
		cutoff_k0 = k0_Cutoff.getValue();
		cte_k1 = k1_Cte.getValue();
		cutoff_k1 = k1_Cutoff.getValue();
		cutoff_k2 = k2_Cutoff.getValue();

    	// Number of eigenvectors to be computed
        nevec_fact = Nevs.getValue();
        if(nevec_fact <= 0) // checking
        {
        	printf("imode> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
        	exit(1);
        }

        verb = Verb.isSet();
        if(Norm.isSet())
        {
        	norm_modes = true;
			if(parse_verb)
				printf("Parser input: Normal modes will be normalized into unit vectors.\n");
        }
        if(NoMass.isSet())
        {
    		nomass_switch = true;
			if(parse_verb)
				printf("Parser input: Masses will be set to 1.0.\n");
        }

		if(NoTors.isSet())
		{
			notors_switch = true;
			if(parse_verb)
				printf("Parser input: Aditional Torsional Springs DISABLED!\n");
		}

		// Setting model and chi
    	model = Model.getValue();
		if(Chi.isSet())
			type = 2; // phi,chi,psi
		else
			type = 0; // phi,psi

		contacts = Contact.getValue(); // Contact method
		if( contacts == 3 && !Funcfile.isSet() ) // checking
		{
        	printf("Parser error, you should include a Funcfile (see: --func)!\nForcing exit!\n");
        	exit(1);
		}
		if(SSfile.isSet())
		{
			ss_switch = true;
			strcpy(file_ss,((temp=SSfile.getValue()).c_str()));
			if(parse_verb)
				printf("Parser input: Secondary Structure File: --ss = %s\n",file_ss);
		}
		if(Funcfile.isSet())
		{
			contacts = 3; // override "Contact"
			func_switch = true;
			strcpy(file_func,((temp=Funcfile.getValue()).c_str()));
			if(parse_verb)
				printf("Parser input: Secondary Structure and Topology Functions File: --func = %s\n",file_func);
		}
		if( contacts == 4 && model != 0 ) // ED-NMA only valid for CA-model
		{
			printf("Parser> At this moment, the edNMA potential is only valid for CA-model!\nForcing exit!\n");
			exit(1);
		}
//		if(contacts == 5 && !KFile.isSet() ) // checking
//		{
//        	printf("imode> Error, you should include a KFile (see: --Kfile)!\nForcing exit!\n");
//        	exit(1);
//		}
		if( KFile.isSet() )
			contacts = 5;

		if (Def.isSet())
        	deform_switch = true;

		// Selecting Output-Modes Coarse-Graining model
    	if(ModelOut.isSet() && !ChiOut.isSet() ) // model_out is set and chi_out not-set
    	{
    		modelf = ModelOut.getValue();  // set model_out
       		typef = type; // set chi_out
    	}
    	else if(!ModelOut.isSet() && ChiOut.isSet() ) // model_out is not-set and chi_out set
    	{
    		typef = 2; // set chi_out
    		modelf = model; // set model_out
    	}
    	else if(ModelOut.isSet() && ChiOut.isSet() ) // both are set
    	{
    		modelf = ModelOut.getValue();  // set model_out
    		typef = 2; // set chi_out
    	}

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
        save_fixfile = SaveFixFile.isSet();
        if(FixRand.isSet())
        {
        	fixRand_prob = FixRand.getValue();
       		fixmodel = 5; // only dihedral coordinates would be fixed
        }
        if(FixRand2.isSet())
        {
        	fixRand_prob = FixRand2.getValue();
       		fixmodel = 4; // all internal coordinates would be fixed
        }
        if(FixSS.isSet())
        {
    		strcpy(fix_ss,((temp=FixSS.getValue()).c_str())); // Gets Fix-string
    		fixmodel = 6;
        }
        if(InEvec.isSet())
    	{
    		strcpy(file_inevec,((temp=InEvec.getValue()).c_str())); // Gets Eigenvector-file
    		inevec_switch = true;
    		save_dh = false; // Disabling IC modes saving... (it is senseles if InEvec enabled...)
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
			  if ((fp = fopen("/dev/urandom", "r")) == NULL)
			  {
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
		if(parse_verb)
			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);

		hessian_method = Hessian.getValue(); // Hessian matrix computation method (0=K-matrix, 1=VW-arrays, 2=Fast)
		kinetic_method = Kinetic.getValue(); // Kinetic energy matrix computation method (0=K-matrix, 1=VW-arrays, 2=Fast)
		conv_method = Convert.getValue(); // Conversion method from ICS to CCS (0=K-matrix or 1=VW-arrays)
		debug_code = Debug.getValue(); // Debugging code

		if(DelHeteros.isSet())
        	delHeteros_switch = true; // Delete heteroatoms
		if(KeepHydrogens.isSet())
			delHydrogens_switch = false; // Keep hydrogens
		if(KeepWaters.isSet())
        	delWaters_switch = false; // Keep waters

		if(SaveWCart.isSet())
		{
	    	save_wcart = true;
	    	save_wcart2 = true;
		}

	    if(SaveCart.isSet())
	    {
	    	save_cart = true;
	    	save_cart2 = true;
	    }

	    if(SaveCA.isSet())
	    {
	    	save_cart = true;
	    	save_CA = true;
	    	save_CA2 = true;
	    }

	    if(SaveMatrix.isSet())
	    	save_matrices = true;

	    if(JustMatrices.isSet())
	    {
	    	just_matrices = true;
	    	save_matrices = true;
	    }

	    swapmatrix_mem = SwapMatrixMem.getValue();

        Ta = Temp.getValue();
        dc = DC.getValue();
        if(dc <= 0.0)
        	dc = 0.001;
	    if(SaveCovar.isSet())
		{
	    	save_cart = true; // enables IC to CC eigenvector translation
	    	save_covar = true;
			if(parse_verb)
				printf("Parser input: Computing and Saving covariance matrix (T=%fK).\n",Ta);
		}
	    if(SaveDCovar.isSet())
		{
	    	save_cart = true; // enables IC to CC eigenvector translation
	    	save_CA = true; // enable CA-eigenvectors translation
	    	save_dcovar = true;
			if(parse_verb)
				printf("Parser input: Computing and Saving the CA-based distance-covariance matrix (T=%fK).\n",Ta);
		}
	    if(SaveCovarText.isSet())
		{
	    	save_cart = true; // enables IC to CC eigenvector translation
//	    	save_covar = true;
	    	save_matrices_text = true; // enable text format
			if(parse_verb)
				printf("Parser input: Saving covariance matrix in plain text format.\n");
		}
	    else
			if(parse_verb)
				printf("Parser input: Saving covariance matrix in binary packed storage format.\n");

	    if(SaveCovarOut.isSet())
		{
	    	save_cart = true; // enables IC to CC eigenvector translation
	    	save_covarf = true;
			if(parse_verb)
				printf("Parser input: Computing and Saving output model covariance matrix (T=%fK) in binary packed storage format.\n",Ta);
		}

	    if(NoModel.isSet())
	    	nomodel = true;

		if(NThreads.isSet())
		{
			nthreads = NThreads.getValue();
			printf("Parser> Using %d threads.\n",nthreads);
		}

		eigensolver = Eigensolver.getValue();
		printf("Parser> Using eigensolver: %d \n",eigensolver);

		if(eigensolver == 2)
		{
			squared_matrices = true;
			printf("Parser input: Using squared matrices\n");
		}

		if(InterMolec.isSet())
		{
			intermolec_factor = InterMolec.getValue();
			intermolec_switch = true;
		}

		if( Verb.isSet() )
		{
			verb = Verb.getValue();
			if(parse_verb)
				printf("Parser input: Verbose level: --verb = %d\n",verb);
		}

   } catch ( ArgException& e )
        { std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}

