/*************************************************************************
 *                        iMORPH/iMODFIT                                 *
 *************************************************************************
 * This program is part of iMOD: http://chaconlab.org/imod/index.html    *
 * (c) Jose Ramon Lopez-Blanco, Jose Ignacio Garzon and Pablo Chacon.    *
 * IQFR-CSIC's Structural Bioinformatics Group (2004-2011).              *
 *************************************************************************
 *                                                                       *
 *   NMA-based fitting/morphing program for Electron Microscopy maps     *
 *   and PDBs.                                                           *
 *   It applies linear combinations of selected Internal Coordinates     *
 *   Normal Modes to deform iteratively the input PDB, in order to       *
 *   minimize its cross-correlation/rmsd vs. the target map/structure.   *
 *                                                                       *
 *************************************************************************
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#  define CLOCKS_PER_SEC  1000000l // Check this again on new versions...

/*= INCLUDES ===================================================================================*/
#include <stdio.h> // needed by some compilers
#include <stdlib.h>
#include <cmdl/CmdLine.h> // Command-line parser
#include <libvolume/include/floatops.h> // ADP's floating point operations library (MUST BE ON TOP !!! ???)
#include <libpdb/include/Macromolecule.h> // ADP's floating point operations library (MUST BE ON TOP !!! ???)
#include <libnma/include/libnma_time.h> // Real-timer (Santi's)
#include <libnma/include/libnma_misc.h> // Mon's NMA related library
#include <libnma/include/libnma_io.h> // Mon's NMA Input-Output library
#include <libnma/include/libnma_cg.h> // Mon's NMA Coarse-Graining library
#include <libnma/include/libnma_move.h> // Mon's NMA IC-Motion library
#include <libnma/include/libnma_deriv.h> // Mon's NMA Derivatives library
#include <libnma/include/libnma_hessian.h> // Mon's NMA Derivatives library
#include <libnma/include/libnma_kinetic.h> // Mon's NMA Derivatives library
#include <libnma/include/libnma_diag.h> // Mon's NMA Derivatives library
#include <libpdb/include/SS_DISICL.h> // Pablo SS


#include <libtools/include/Io.h>
#include <libnmafit/include/nmafit.h>

#ifdef GAUSS
#include <libgausscorr/include/gausscorr.h>
#endif

//#include "libvolume/floatops.h" // ADP's floating point operations library (MUST BE ON TOP !!! ???)
//#include "tbb/tbb.h"

// Pre-processor definitions
// #define FITTING // "FITTING" CODE (#if not defined --> MORPHING)
// #define USE_TBB // Intel's "TBB" parallel routines enabled


/*= GLOBAL VARIABLES ===========================================================================*/
char *prog_name; // pointer to program name
char *prog_name2; // pointer to program name
char *prog_url; // pointer to program URL
char VERSION_NMAFIT[10]="v1.9";
char file_initial[FILE_NAME];
char file_ref_pdb[FILE_NAME];
char file_ref_vol[FILE_NAME];
char name[FILE_NAME];
char file_ptraj[FILE_NAME];
char file_out[FILE_NAME];
char file_score[FILE_NAME];
char text[FILE_NAME]; // temporal text array
char file_ss[FILE_NAME];
char file_func[FILE_NAME];
char fix_file[FILE_NAME];
char fix_ss[30]; // string which allocates single character SS identifiers
char file_aln[FILE_NAME]; // Sequence alignment file name
char dummy_string[FILE_NAME];

int contacts = 0; // Contacting method
float cte_k0 = 1.0;
float cte_k1 = 1.0;
float x0 = 3.8;
float power = 6.0;
float cutoff_k0 = 10;
float cutoff_k1 = 10;
float cutoff_k2 = 10;
float intermolec_factor = 1; // Inter-molecular force constant factor
double gauss_c = 20.0; // Distance parameter for gaussian weight (wrmsd)

int eigensolver = 1; // Eigensolver choice: 0= dspgvx, 1= dsdrv1_AP_BP_W_mon
int nthreads = 0; // Number of threads used for parallelization
int nevs = -1;
int addnevs = 0; // Number of added modes (in any case, this value will be added into "nevs", once computed)
int nex = 0.1; // Number of modes excited per iteration
int nex2 = 0.2; // Higher limit of modes excited per iteration
int excite;
int aln_level = 2; // Clustal alignment level  (1, 2 or 3 to match *, *: or *:. ,respectively)
float nevec_fact = -1; // % of eigenvectors to be computed (set by parser)
float nevs_ratio;
float nex_fact = -1; // % of eigenvectors computed to be excited (set by parser)
float nex_fact2 = -1; // % of eigenvectors computed to be excited (set by parser)

double dice;
double step = 8; // Constant step used when "rand_step_switch" and "step_end_switch" are FALSE
double step_begin = step; // Use with "rand_step_switch" and "step_end_switch"
double step_end = 2; // Use with "rand_step_switch" and "step_end_switch"
int type = 0; // Internal Coordinates Coarse-Graining type
int model = 1; // Model Coarse-Graining
int fixmodel = 0; // =0 --> no ICs fixation
float fixRand_prob=0.5; // random IC fixation probability
bool savemodel_switch = false; // = true --> save current coarse-grained-model initial PDB
bool savefitted_switch = false; // = true --> save fitted coarse-grained-model final PDB
bool saveref_switch = false; // = true --> save pdb_ref during PDB-MAP fitting
bool savemaps_switch = false; // = true --> save reference map (experimental), initial map (projected and filtered from input pdb) and fitted map.
bool notors_switch = false; // true --> Disables TORSional springs addition (hessianMDHx)
bool gauss_switch = false; // true --> Enable Gaussian overlap scoring
bool gauss_debug = true; // write/dump gaussian related stuff...
bool intermolec_switch = false; // Enables inter-molecular force constant factor
bool squared_matrices = false; // true --> Using squared matrices instead of the memory-efficient packed triangular storage.
bool eigensolver_switch = false; // true --> automatic eigensolver choice overridden.
float gauss_bound = 5.0; // Factor to scale Gaussian bounding boxes

// About Maps & Filtration
int filter=0; // sets the filtration method ( 0=test, 1=Fourier, 2=Kernel )
float resolution;
float gauss_sd;
float ref_thr;
float model_thr;

int max_iter;
double convergence = -1000; // very low default value for no convergence

int prob_method; // sets the normal mode probability selection method
double rediag = 0.1; // minimum score increment to trigger diagonalization

// Pose Refinement
float refine_tol = 0.125; // Thresold (tolerance) to avoid unsignificant movements (Rot.Tras-lational)
float shake_rot=1 * (M_PI/180); // radians // Rotational Shake amplitude (rad)
float shake_mov=2; // Traslational Shake amplitude (Amstrongs) (define depending on sampling, when using volumes)
float refine_prob = 0.05; // Rotational & Traslational refinement chance probability
int refine_each=200; // Each "refine_each" iters Rot.&Trans. will be refined! (Set to "-1" to enable "refine_prob")
//int refine_iter = 4; // Number of iterations per Rot./Trans.-lational refinement

float delta_save = 0.5; // score increment (normalized) to save frames
unsigned int seed; // random number generator seed

// remove W-corr
float w_corr=0.5;
float w_corr0=w_corr; // weight factor for correlation map

int conv_length=1000; // Segment to compute the averaged score
float trans_factor=0.1; // Defines, relatively, the number of transition iterations (correlation method)
float prob_cutoff=0; // Sets the Cutoff Threshold to take into account probabilities (for NM choice)
int corr_radius=2; // Sets the sphere "radius" to compute local cross-correlation (in voxels)
int score_method=1; // Sets the Cross-correlation method (score, PDB-VOL)
float min_prob=0.0; // Probability profile minimum value
float min_wrmsd=0.0; // Minimum RMSD Weighting
float norm_factor; // normalization factor (sigma) (It will be constant during all run!)
int verbose=0; // Sets the Verbose level( 0=minimum, 1=medium, 2=maximum)
int each_verb = 1; // Shows verbose each X iterations
int each_score = 20; // Saves scores each X iterations
float diag_choice_cutoff = 1.01; // If "nevs_ratio"(nevs) <= diag_choice_cutoff then nmad_dsygvx, else  nmad_dsygvd
float sigma_kernel = 3; // Kernel filter sigma (usually = 3)
float pad_fourier = 16.0; // 16.0; // Pablo's Padding Fourier Factor = 4.0
double pose_t0 = 0.0; // 6D pose refinement initial temperature
double pose_t; // 6D pose refinement temperature
double pose_ka = 1.0; // 6D pose refinement Ka cooling parameter
int pose_cool = 1; // 6D pose refinement cooling schedule
int sel_chain = -1; // if <0 the whole molecule will be 6D refined
//double p_x0 = 0.05;
//double p_s = 6;

bool unitmass_switch = false; // true --> Constant mass for every atom (=1.0)
bool unitnelec_switch = false; // true --> Constant number of electrons for every atom (=1.0)
bool rand_step_switch = false; // random step (within selected range), else "lineal step increase/decrease"
bool rand_excited_switch=false; // "nex2" must be also set!
bool pdb_switch=false; // Reference input method (PDB or VOLUME)
bool bench_switch=false; // To allow RMSDs calculation when reference PDB available! (for benchmark use)
bool nomodel = false; // Allows 3BB2R PDB input (initial model)
bool movie_switch = false; // Save the trajectory movie
bool pdblog_switch = false; // Enables PDB output each time it is requested
bool corr_nozero = false; // false = Frame / true = Target mask (!=0)
bool corrmap_switch = false; // Weighted Cross-Correlation Map Switch
bool wrmsd_switch = true; // Weighted RMSD switch
bool norm_switch = false; // Forces map normalization
bool out_files = false; // true = outputs many debug files...
bool time_switch = false; // true = output timers info to file (set in parser)
bool backward_switch = false; // false = only evaluates the "forward" movement
bool pad_switch = true; // true = pads volumes to fit motion
bool poseRand_switch = false; // true = random orientations pose refinement
int refine_delay = -1; // 6D pose refinement delay (number of iterations)
bool ss_switch = false; // the NMA force constants will be set according to SS rules
bool func_switch = false; // input file with function coefficients (Topology and SS)
bool debug = false;
bool fullatom_switch = false; // true = Full-atom output
bool deltasave_rmsd_switch = true; // true --> Save output when: delta_save < (delta-RMSD; respecting last saved)
// false --> Save output when: delta_save < (delta-Normalized_Score)
bool randweight = true; // true --> Weights randomly (in "select_mode_prob") each excited mode.
// false --> Excited mode weights will be either +1.0 or -1.0 (only random sign).
bool norm_evec_switch = false; // true --> norm=1 merged mode normalization
bool save_deviation = false; // true --> with PDB-VOL benchmark or PDB-PDB, it stores deviations in B-factors
bool more_rmsds = false; // true --> CA rmsd output
bool parse_verb = false;
bool fast_switch = false; // true --> Enable "fast" PDB into MAP projection method (not-trilinear)
bool dump_switch = false; // true --> Enables some verbose...
bool scv_weight = false; // true --> Enables SCV mode weighting method during mode merging.
bool aln_switch = false; // true --> Using only common atoms between initial and target atomic structures (morphing only)
bool server = false; // =true, to enable "server mode" (maximize the default automated options)
bool delHydrogens_switch = true; // Delete hydrogens
bool delHeteros_switch = false; // Delete heteroatoms
bool delWaters_switch = true; // Delete waters
bool chimera_switch = false; // Enable specific output for Chimera (e.g. dump current movie frame to <basename>_fitted.pdb

timer time_main; // timer for long timings
timerReal ht_stage,ht_substage; // ,ht_ali,ht_diag,ht_score,ht_mov,ht_wrmsd,ht_choice;
time_t t_stage=0;
float t_main=0,t_ali=0,t_nma,t_score=0,t_mov=0,t_select=0,t_merge=0,t_diag=0,t_kine=0,t_hess=0,
		t_mov_nm=0,t_mov_ali=0,t_mov_proj=0,t_next=0,t_sum=0,t_after=0,t_corr=0,t_springs=0,t_trial=0;

convoluteK_data *threads_data = NULL;
pthread_t *threads = NULL;

//extern CRandomMersenne * rg; // Mersenne Twister global object
// because it's already declared in "libnmadoc.cpp"

/*==============================================================================================*/
// input variables
using namespace TCLAP;
using namespace std;
void parseOptions(int argc, char** argv);
void parseOptionsMorph(int argc, char** argv);
void parseOptionsFit(int argc, char** argv);
char *printbool(bool var)
{
	if(var)
		return("true");
	else
		return("false");
}
/*==============================================================================================*/

//// Check maximum inter-atomic distance between two macromolecules, and returns its value.
//float checkmax(Macromolecule *mol1, Macromolecule *mol2)
//{
//	bool debug=true;
//	pdbIter *iter1 = new pdbIter( mol1 );
//	pdbIter *iter2 = new pdbIter( mol2 );
//	Tcoor pos1,pos2;
//	float current,max=0.0;
//	int i1,i2;
//
//	for ( iter1->pos_atom=0, iter2->pos_atom=0; !iter1->gend_atom(); iter1->next_atom(), iter2->next_atom() )  // screens all-atoms
//	{
//			( iter1->get_atom() )->getPosition(pos1);
//			( iter2->get_atom() )->getPosition(pos2);
//			current = sqrt( pow(pos1[0]-pos2[0],2) + pow(pos1[1]-pos2[1],2) + pow(pos1[2]-pos2[2],2) ); // distance 1-2
//			if(current > max)
//			{
//				max = current; // store maximum distance
//				i1=iter1->pos_atom; // store atomic indices
//				i2=iter2->pos_atom;
//			}
//	}
//
//	if(debug)
//		fprintf(stderr,"Msg(checkmax): d= %f  i1= %d  i2= %d\n",max,i1,i2);
//
//	delete(iter1);
//	delete(iter2);
//	return(max);
//}

// *****************************************************
// *                        MAIN                       *
// *****************************************************
int main(int argc, char **argv)
{
	// TESTING VARIABLES
	//	float current_dmax=0.0,max_dmax=0.0,avg_dmax=0.0;
	//	int iter_dmax=0;
	// END TESTING VARIABLES

	// Error messages for Chimera plugins
	char msg_missing_atoms[] = "Error, missing atoms detected in PDB! Please, you can either:\n\t"
			"1) select a less-restrictive atomic model, e.g. CA, or\n\t"
			"2) fix missing atoms, e.g. using Profix or Rosetta.\n";

	// Parsing Input (THIS MUST BE CALLED FIRST!)
#ifdef FITTING
	parseOptionsFit(argc,argv);
#else
	parseOptionsMorph(argc,argv);
#endif

	// SELECTING CA MODELS (just to compute RMSD)
	// Conditions
	Condition * calpha = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	calpha->add( " CA " );
	calpha->add( " P  " );
	Conditions * calpha2 = new Conditions();
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
		printf("Sorry, unable to open COMMAND-LOG FILE: %s\n",text);
		exit(1);
	}
	fprintf(f_com,"# %s> Welcome to %s %s\n",prog_name,prog_name2,VERSION_NMAFIT);
	bool already_seed=false;
	for(int i=0; i<argc; i++)
	{
		fprintf(f_com,"%s ",argv[i]);
		if(strcmp(argv[i],"--seed") == 0 )
			already_seed = true;
	}
	if(!already_seed)
		fprintf(f_com,"--seed %u\n",seed); // This allows user to carry out again the same run.
	else
		fprintf(f_com,"\n");

	if(eigensolver == 2)
	{
		squared_matrices = true;
		printf("Parser input: Using squared matrices\n");
	}

#ifdef USE_TBB // Enables Intel's Threading Building Blocks parallel routines
	tbb::task_scheduler_init init;  // Automatic number of threads
	if(nthreads > 0)
	{
		init.terminate();
		init.initialize(nthreads);  // Explicit number of threads
		printf("%s> Using TBB's %d threads from the %d maximum available.\n",prog_name,nthreads,tbb::task_scheduler_init::default_num_threads());
	}
	if(verbose > 0)
		if(init.is_active())
			printf("task_scheduler_init is active!\n");
		else
			printf("warning: task_scheduler_init is NOT active!\n");
#endif


	// Mersenne Twister Seed Initialization
	// It Outputs a random double (0 <= x < 1) with Random()
	rg = new CRandomMersenne( seed );

	int num_res,num_resR,num_atoms,num_atomsR,index,num_atomsNCAC;
	float *coord=NULL;
	float *coordNCAC=NULL;
	double *dist_matrix=NULL;
	trd *der=NULL;
	Atom *atom;
	Residue *res;
	Tcoor pos;
	int index_atom;
	float factor;
	int nev=0;
	double dummyd;
	// Flag...
	bool out_switch=false; // Outputs files...
	float matrix4[4][4];
	double last_diag=0; // RMSD from last diagonalized macromolecule
	vlVolume *vol_dummy,*volR,*volT,*volF,*volB,*vol_temp,*volD,*volD2;
	float rot[3],mov[3],trans[3]; // Needed for Docking & Shaking
	double *deviation=NULL; // array with each residue RMSD or Cross.Corr.
	float *relev = NULL;
	mode *selected = NULL; // stores selected modes
	twid *ipas;
	int nipa;
	float *kernel=NULL; // Filtration
	int dim_vox=0; // Filtration
	float score_lastsaved=0;
	float rmsd;
	bool accepted; // Simulated Annealing aceptation boolean
	double entropy=0,cost=0,delta_score; // Thermodynamic SA (TSA)
	vlVolume **next_vol;
	Macromolecule *next_mol,*next_mol2,*next_molT,*next_molT2;
	Macromolecule *dummy;
	double next_score;
	vlVolume *corrmap=NULL,*acorrmap=NULL;
	float score_avg1=0,score_avg2=0,delta_avgscore; // To check for convergence
	float *scores=NULL;
	int stage=0; // needed to perform different optimization stages
	bool do_refine = false; // Refinement decision
	bool trans_switch = false; // Transition between correlation methods switch
	int trans_iters;
	float trans_m=0;
	bool do_main_loop = true;
	float *wrmsd_prof=NULL;
	//	float first_corr=0;
	double first_corr=0.0;
	vlDim pad; // Defines "padding" during all program (speed-up)
	//	float avgR,sigR; // Store the Target average and sigma. (speed-up)
	double avgR,sigR; // Store the Target average and sigma. (speed-up) Double-precision necessary!
	float shake;
	int diag_info;// To handle diagonalization failures!
	int backup_iter=5; // number of re-diagonalizations before exiting due to diag. failures.
	//	double rand_backup = 0.3; // 70% prob. to save a backup model

	vlPoint3f Rmax,Tmax; // Corner positions (real space)
	vlPoint3f Rorig, Torig; // Target and Initial map origins.
	vlPoint3f pmin1,pmax1,pmin2,pmax2; // New Map corner positions
	vlDim Rdim,Tdim; // Map sizes
	vlUnit Rvox,Tvox; // Samplings (voxel size)
	timerReal ht_timer;
	timer t_timer; // timer
	float t_test1,t_test2;
	M4Rot *matrix4_op; // needed

	// COMMAND
	if( verbose > 0 )
	{
		printf("%s> INPUT COMMAND:\n%s> ",prog_name,prog_name);
		for(int i=0; i<argc; i++)
			printf("%s ",argv[i]);
		printf("\n%s>\n",prog_name);
	}

	// Initialize aminoacids and nucleotids
	init_aminoacids();

	//	// Put inside init_aminoacids();
	//	for(int i=0; i<N_AMINO; i++)
	//	{
	//		AA[i].nelec = 0.0;
	//		for(int k=0; k<AA[i].nheavyatoms; k++)
	//		{
	//			int atn = atom_types[AA[i].atom[k].fullatom_type-1].at;
	//			AA[i].nelec += Table_Elements::getElement(atn)->number;
	//		}
	//	}

	// Target Macromolecule pointer (needs to be outside!)
	Macromolecule *molr = new Macromolecule( "Target" );

	// Reading initial pdb
	if( verbose > -1 )
		printf( "%s> Model PDB file: %s\n",prog_name,file_initial );

	// Needed further to convert normal modes into the appropriate atomic model
	Macromolecule *molt,*molt2,*molr2;
	Macromolecule *moltini = new Macromolecule("Initial");
	Macromolecule *molrini = new Macromolecule("Reference");

	moltini->readPDB(file_initial); // Reading initial PDB
	if(delHydrogens_switch)
	{
		printf( "%s> Deleting Hydrogen atoms (if any) from Initial PDB...\n",prog_name );
		moltini->deleteHYDS();
	}
	if(delHeteros_switch || model == 0) // CA model can't deal with HETATM's... (TO DO)
	{
		printf( "%s> Deleting Hetero-atoms (if any) from Initial PDB...\n",prog_name );
		moltini->delete_heteros();
	}

	if(delWaters_switch)
	{
		printf( "%s> Deleting Water molecules (if any) from Initial PDB...\n",prog_name );
		moltini->delete_waters();
	}
	moltini->info(stdout); // checking initial model pdb

	// PRE-PROCESSING INITIAL MODEL (FORMAT)
	Macromolecule *mol0; // 0-chain, needed for full-atom
	Macromolecule *mol0_full; // full-chains, needed for full-atom in nma_type=2
	Macromolecule *mol0_model; // 0-chain-model, needed for full-atom
	tri *props0;
	int *unat0;
	int size0=0,size0_old;
	bool *fixed0=NULL;

	double *bf_mol0; // B-factors profile for Full-Atom model
	if(!nomodel) // Formatting (sorting) residue atoms (if it's needed)
	{
		// Formatting target pdb
		if( verbose > 1 )
			printf( "%s>  1) Format Initial Model residue names\n",prog_name );
		// moltini->format_residues(3); // non-hydrogen heavy atoms

		// Formatting target PDB and checking missing atoms according to selected CG-model
		if(moltini->format_residues(false, model) > 0)
		{
			printf( "%s> %s\n", prog_name, msg_missing_atoms);
			exit(2);
		}

		// Holds a copy of the formated HA input PDB
		if(fullatom_switch)
		{
			mol0 = new Macromolecule( moltini ); // Formatted copy of the initial HA input PDB
			properMFA(mol0,&props0,&unat0,type,2); // Computing its "properties"

			//			// Get B-factors profile for current CG-model. It should be applied later before saving each movie frame and initial or final model.
			//			bf_mol0 = (double *) malloc(mol0->get_num_atoms() * sizeof(double));
			//			check_pointer(bf_mol0,"B-factors profile for Full-Atom model"); // this checks the memory allocation
			//			mol0->get_Pdbfact(bf_mol0);

			// Computing Total number of degrees of freedom for the initial HA input PDB
			size0 = num_dofs(mol0,props0);
			printf( "%s> Number of Heavy Atom (HA) model Internal Coordinates: %d\n",prog_name,size0);
			size0_old = size0;

			// Allocating memory for the "fix" array of the initial HA input PDB
			fixed0 = (bool *)malloc(sizeof(bool)*size0);
			check_pointer(fixed0,"Fixation array (initial)"); // this checks the memory allocation
		}
	}

	// SETTING COARSE-GRAINED MODEL
	switch(model)
	{
	case 0: // CA-model model: CA + (NH and CO)-terminal model
	{
		printf( "%s> Coarse-Graining model: CA-model\n",prog_name);

		// N,CA,C selection
		molt2 = moltini->select( ncac2 ); // molt --> moving model
		num_atomsNCAC = molt2->get_num_atoms();

		// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
		// Warning, atoms are not copied, they're just pointers to the original atoms.
		molt = cg_CA( molt2, !nomodel, unitmass_switch, unitnelec_switch, !pdb_switch ); // masses and number of electrons will be computed

		// Saving "NCAC" model
		if(savemodel_switch)
		{
			sprintf(dummy_string,"%s_ncac.pdb",name);
			molt2->writePDB( dummy_string );
		}
		break;
	}
	case 3: // N,CA,C-model (C3-model)
	{
		printf( "%s> Coarse-Graining model: N,CA,C-model (Experimental)\n",prog_name);
		// N,CA,C selection
		molt = moltini->select( ncac2 );
		if(!nomodel) // Adds masses (if it's needed)
			mass_NCAC( molt, unitmass_switch, true, unitnelec_switch, !pdb_switch ); // masses and number of electrons will be computed
		break;
	}
	case 1: // C5 model
	{
		printf( "%s> Coarse-Graining model: C5\n",prog_name);
		molt = new Macromolecule(moltini);
		if(!nomodel) // Makes 3BB2R model (if it's needed)
		{
			// CREATES a 3BB2R reduced model
			//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
			//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
			//     There are a few exceptions: for Ala, Gly and Pro,
			//     and for the 1st and last residues.
			printf("%s> Creating C5 reduced model:\n",prog_name);
			cg_3BBR2(molt, unitmass_switch, unitnelec_switch, !pdb_switch);
		}
		else
			molt = moltini;
		break;
	}
	case 2: // HA
	{
		printf( "%s> Coarse-Graining model: HA \n",prog_name);
		molt = new Macromolecule(moltini);
		molt = moltini;
		mass_FA(molt, unitmass_switch, unitnelec_switch, !pdb_switch); // sets masses
		break;
	}
	}
	num_atoms = molt->get_num_atoms();
	num_res = molt->get_num_fragments();
	printf( "%s> Selected model number of residues: %d\n",prog_name,num_res );
	printf( "%s> Selected model number of (pseudo)atoms: %d\n",prog_name,num_atoms );

	// Holds a copy of the reduced (or not) model input PDB
	if(fullatom_switch)
	{
		if(model==0)
			mol0_model = new Macromolecule( molt2 ); // C3-model (NCAC)
		else
			mol0_model = new Macromolecule( molt );
	}

#ifdef GAUSS

	// GAUSSIAN MIXTURE MODEL RELATED STUFF... (Future work)
	GAUSS3D *gmolt,*gforward,*gbuff;
	// Make gaussian molt
	if(gauss_switch)
	{
		// See what Kovacs' DDFF paper sais about: Tama, Miyashita, Brooks III. (2004)
		// s = R/(2*sqrt(3)) --> 2*s^2 = (R^2)/6     (R= resolution)
		gauss_sd = pow(resolution,2)/6;

		gmolt = pdb2gauss(molt,gauss_sd);
		if(gauss_debug)
		{
			sprintf(dummy_string,"%s_molt.gauss",name);
			Write_Gaussian3D_File(dummy_string,num_atoms,gmolt,'I',"");
		}
		// Computing maximum size of gaussians to define their bounding boxes
		printf("Computing maximum size of gaussians to define their bounding boxes\n");
		gaussian_sizes(gmolt,num_atoms,gauss_bound);
	}
#endif

	// Saving initial (current CG-model) PDB
	if(savemodel_switch && !(chimera_switch && fullatom_switch) )
	{
		sprintf(dummy_string,"%s_model.pdb",name);
		molt->writePDB( dummy_string, false ); // warning, this renumbers the PDB
	}

	//	// Get B-factors profile for current CG-model. It should be applied later before saving each movie frame and initial or final model.
	//	double *bf_molt; // B-factors profile for current CG-model.
	//	if( !(bf_molt = (double *) malloc(num_atoms * sizeof(double))) )
	//	{
	//		printf("Error in memory allocation! Forcing exit!\n");
	//		exit(1);
	//	}
	//	molt->get_Pdbfact(bf_molt);

	// CURRENT MODEL - Selecting the first chain (zero)
	Macromolecule *molt0; // first chain (zero) selected by "select_chain"
	Macromolecule *molt02; // first chain (zero) NCAC selected by "select_chain"
	Macromolecule *p_molt, *p_molt0, *p_molt2, *p_molt02; // buffers

	char *ss_table; // SS-table
	TSfunc *funcs; // SS and Topology functions
	int nfunc=0; // number of functions
	int nss=0;

	// Reading Secondary structure file.
	if(ss_switch) // NMA with SS (optional)
	{
		// Reading ss-file (allocating table memory)
		read_ss(file_ss,&ss_table,&nss);
		// Some checking...
		if( nss != num_res)
		{
			printf("%s> Sorry... the SS-file should have the same number of entries (%d) as the number of residues (%d) in the input PDB\n",prog_name,nss,num_res);
			exit(1);
		}
	}
	else if(contacts == 3 || fixmodel == 6)
		ss_table = molt->secondary_structure(true); // Erney's SS method
	else
		ss_table = NULL; // disables SS-checking in "make_ipasTS" (Topology only)

	// Reading T/SS Functions.
	if(contacts == 3)
	{
		// Reads input TS functions (mandatory)
		if(func_switch)
			read_TSfunc(file_func, &funcs, &nfunc);
		else
		{
			printf("%s> Sorry... At this moment you should include always a connection functions file\n",prog_name);
			exit(1);
		}

		if(verbose > 0)
			for(int i=0; i<nfunc; i++)
				printf("%s> %c%c %d %f %f %f\n",prog_name, funcs[i].i, funcs[i].j, funcs[i].t, funcs[i].a, funcs[i].b, funcs[i].c );
	}

	// REFERENCE MODEL - Selecting the first chain (zero)
	Macromolecule *molr0,*molr02; // first chain (zero) selected by reference
	Macromolecule *p_molr, *p_molr0, *p_molr2, *p_molr02; // buffers


#ifdef GAUSS

	// GAUSSIAN MIXTURE MODEL RELATED STUFF... (Future work)
	GAUSS3D *gmolr;
	int num_gauss;
	// Reading reference
	if(gauss_switch)  // PDB-Gaussian-VOLUME Fitting
	{
		// Read gaussian mol from target file
		//		Read_Gaussian3D_File(file_ref_pdb, &num_atomsR, gmolr);
		readGAUSS3D("1oel_map100.gauss", &num_gauss, &gmolr);

		// Computing maximum size of gaussians to define their bounding boxes
		printf("Computing maximum size of gaussians to define their bounding boxes\n");
		gaussian_sizes(gmolr,num_gauss,gauss_bound);

	}
#endif

	//	else
	//	{
	// Processing the Reference PDB (reading, formatting, coarse graining)
	if(pdb_switch || bench_switch )  // PDB-PDB Morphing || Benchmark Fitting
	{
		// REFERENCE PDB
		// Reading reference PDB
		if( verbose > -1 )
			printf( "%s> Target PDB file: %s\n",prog_name,file_ref_pdb );
		molrini->readPDB(file_ref_pdb );
		if(delHydrogens_switch)
		{
			printf( "%s> Deleting Hydrogen atoms (if any) from Target PDB...\n",prog_name );
			molrini->deleteHYDS();
		}
		if(delHeteros_switch || model == 0) // CA model can't deal with HETATM's... (TO DO)
		{
			printf( "%s> Deleting Hetero-atoms (if any) from Initial PDB...\n",prog_name );
			moltini->delete_heteros();
		}

		if(delWaters_switch)
		{
			printf( "%s> Deleting Water molecules (if any) from Target PDB...\n",prog_name );
			molrini->delete_waters();
		}
		molrini->info(stdout); // checking initial model pdb

		// Formatting reference (target) PDB
		if( verbose > 1 )
			printf( "%s> Format Target Model residue names\n",prog_name );
		// molrini->format_residues(3); // non-hydrogen heavy atoms

		// Formatting target PDB and checking missing atoms according to selected CG-model
		if(molrini->format_residues(false, model) > 0)
		{
			printf( "%s> %s\n", prog_name, msg_missing_atoms);
			exit(2);
		}


		if(server)
		{
			// SETTING REFERENCE COARSE-GRAINED MODEL
			Macromolecule *target;
			target = new Macromolecule();

			switch(model)
			{
			case 0: // CA-IC model: CA + (NH and CO)-terminal model
			{
				// Just-CA selection
				target = molrini->select_cpy( calpha2 );
				break;
			}
			case 3: // N,CA,C-model (C3)
			{
				// N,CA,C selection
				target = molrini->select_cpy( ncac2 );
				break;
			}
			case 1: // C5 model
			{
				target = new Macromolecule(molrini);
				cg_3BBR2(target, unitmass_switch, unitmass_switch);
				break;
			}
			case 2: // HA model
			{
				target = new Macromolecule(molrini);
				mass_FA(target, unitmass_switch); // sets masses
				break;
			}
			}
			sprintf(dummy_string,"%s_target.pdb",name);
			target->writePDB( dummy_string, false ); // renumbers the PDB
			if( verbose > -1 )
				printf("%s> Target Model:          %35s\n",prog_name,dummy_string);
			delete target;
		}

		// SETTING REFERENCE COARSE-GRAINED MODEL
		switch(model)
		{
		case 0: // CA-IC model: CA + (NH and CO)-terminal model
		{
			if( verbose > 0 )
				printf( "%s> Coarse-Graining model: CA\n",prog_name);

			if(aln_switch)
			{
				// CA selection
				molr2 = molrini->select( calpha2 );
				molr = molrini->select( calpha2 );
			}
			else
			{
				// N,CA,C selection
				molr2 = molrini->select( ncac2 );

				// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
				// Warning, atoms are not copied, they're just pointers to the original atoms.
				molr = cg_CA( molr2, !nomodel, unitmass_switch ); // masses will be computed
			}
			break;
		}
		case 3: // N,CA,C-model (C3)
		{
			if( verbose > 0 )
				printf( "%s> Coarse-Graining model: N,CA,C-model (Experimental)\n",prog_name);

			// N,CA,C selection
			molr = molrini->select( ncac2 );

			if(!nomodel) // Adds masses (if it's needed)
				mass_NCAC( molr, unitmass_switch );
			break;
		}
		case 1: // C5 model
		{
			if( verbose > 0 )
				printf( "%s> Coarse-Graining model: C5\n",prog_name);
			if(aln_switch)
			{
				// CA selection
				molr = molrini->select( calpha2 );
			}
			else
			{
				molr = new Macromolecule(molrini);
				if(!nomodel) // Makes 3BB2R model (if it's needed)
				{
					// CREATES a 3BB2R reduced model
					//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
					//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
					//     There are a few exceptions: for Ala, Gly and Pro,
					//     and for the 1st and last residues.
					if( verbose > 0 )
						printf("%s> Creating C5 reduced model:\n",prog_name);
					cg_3BBR2(molr, unitmass_switch, unitmass_switch);
				}
				else
					molr = molrini;
			}
			break;
		}
		case 2: // HA model
		{
			if( verbose > 0 )
				printf( "%s> Coarse-Graining model: Full-Atom (no coarse-graining)\n",prog_name);
			molr = molrini;
			mass_FA(molr, unitmass_switch); // sets masses
			break;
		}
		}
		num_atomsR = molr->get_num_atoms();
		num_resR = molr->get_num_fragments();
		if( verbose > 0 )
		{
			printf( "%s> Selected model number of residues (target): %d\n",prog_name, num_resR );
			printf( "%s> Selected model number of (pseudo)atoms (target): %d\n",prog_name, num_atomsR );
		}

		// Some checking
		if(num_res != num_resR && !aln_switch)
		{
			printf("%s>\n%s> ERROR: Both PDB files should have the same number of residues!"
					"\n%s> Forcing exit!\n%s>\n",prog_name,prog_name,prog_name,prog_name);
			exit(1);
		}

		// Some checking
		if(num_atoms != num_atomsR && !aln_switch)
		{
			printf("%s>\n%s> ERROR: Both PDB files should have the same number of atoms!"
					"\n%s> Forcing exit!\n%s>\n",prog_name,prog_name,prog_name,prog_name);
			exit(1);
		}

		//			// Make gaussian molr
		//			if(gauss_switch)
		//				gmolr = pdb2gauss(molr,resolution);

		// Saving current model
		if(saveref_switch)
		{
			sprintf(dummy_string,"%s_ref.pdb",name);
			molr->writePDB( dummy_string, false ); // renumbers the PDB
		}
	}
	// PROCESSING THE TARGET MAP (reading, cropping, filtering)
	if( !pdb_switch ) // PDB-VOLUME Fitting
	{

		// Reading target (reference,final) volume
		// (Target volume should be filtered in the same way as during refinement!)
		// (We don't filter target volume at all!)
		if( verbose > -1 )
			printf( "%s> Target Map file: %s\n",prog_name,file_ref_vol);
		volR=FOPS::readFile(file_ref_vol); // read map

		FOPS::threshold( volR, ref_thr ); // apply threshold

		if( verbose > 1 )
			cout << prog_name << ">\t  resolution= " << resolution << "  cutoff= " << ref_thr
			<< "  size= " << volR->dim() << endl;
		volR=FOPS::crop( volR, ref_thr,true); // crop map

		if( verbose > 1 )
			cout << prog_name << "> Cropped Target Volume size " << volR->dim() << endl;

		vlDim pad1, pad2;
		int repeat_times=10; // Filtration is repeated "repeat_times" to improve time statistics.

		// Test to determine the optimal filtration method for current parameters:
		// resolution, size and voxel size.
		if(filter == 0 || filter == 1) // filtration timing test enabled
		{
			// **************
			//  FOURIER TEST
			// **************
			if( verbose > 1 && filter == 0)
				cout << prog_name << "> TESTING FOURIER FILTRATION METHOD" << endl;

			// Projecting initial PDB
			volT = molt->fillVolumeNE_3BB2R( volR->units().x() ); // works also for any CG-model
			if( verbose > 1 )
				cout << prog_name << "> Initial map Box-size " << volT->dim() << endl;

			// INITIAL VOLUME FOURIER PADDING
			pad1 = vlDim( (int)ceil( ( volT->dim().x()+resolution/volT->units().x()+2) / pad_fourier),
					(int)ceil( ( volT->dim().y()+resolution/volT->units().y()+2) / pad_fourier),
					(int)ceil( ( volT->dim().z()+resolution/volT->units().z()+2) / pad_fourier) );

			// INITIAL VOLUME PADDING
			if(pad1.x()!=0 || pad1.y()!=0 || pad1.z()!=0) // if we need padding !!!
			{
				volT = FOPS::padVolume( volT , pad1 ); // vol. must be padded to kernel filtration
				if( verbose > 1 )
					cout << prog_name << "> Padded Initial (and temporal) Volume: " << volT->dim() << endl;
			}

			// INITIAL VOLUME PROJECTION AND FILTRATION
			// (Note that the current volume will be built from the the Reduced 3BB2R model !!!)
			vol_dummy = new vlVolume( volT );
#ifdef USE_PTHREAD // Enables PThread parallel routines
			project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,1,resolution,fast_switch,nthreads,threads_data); // Fourier filtration
#else
			project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,1,resolution,fast_switch); // Fourier filtration
#endif
			delete(vol_dummy); // it will be allocated below (with the definitive size!)
			// SIZE PADDING - The whole movement must be contained !!! (due to sigmas...?)

			volT = FOPS::crop( volT, model_thr, true); // to see the proper size, it must be filtered!
			if( verbose > 1 )
				cout << prog_name << "> Cropped Initial (and temporal) Volume size " << volT->dim() << endl;

			// THE WHOLE MOTION MUST BE CONTAINED INSIDE VOLUMES
			// Setting new map sizes
			Rdim = volR->dim(); // Map sizes
			Tdim = volT->dim();
			Rvox = volR->units(); // Samplings (voxel size)
			Tvox = volT->units();
			volR->getPosition(&Rorig); // Origins
			volT->getPosition(&Torig);
			Rmax.x( (float)Rdim.x() * (float)Rvox.x() + Rorig.x() );
			Rmax.y( (float)Rdim.y() * (float)Rvox.y() + Rorig.y() );
			Rmax.z( (float)Rdim.z() * (float)Rvox.z() + Rorig.z() );
			Tmax.x( (float)Tdim.x() * (float)Tvox.x() + Torig.x() );
			Tmax.y( (float)Tdim.y() * (float)Tvox.y() + Torig.y() );
			Tmax.z( (float)Tdim.z() * (float)Tvox.z() + Torig.z() );
			// Origin (min corner)
			if( Torig.x() < Rorig.x() )
				pmin1.x( Torig.x() );
			else
				pmin1.x( Rorig.x() );
			if( Torig.y() < Rorig.y() )
				pmin1.y( Torig.y() );
			else
				pmin1.y( Rorig.y() );
			if( Torig.z() < Rorig.z() )
				pmin1.z( Torig.z() );
			else
				pmin1.z( Rorig.z() );
			// Max corner
			if( Tmax.x() > Rmax.x() )
				pmax1.x( Tmax.x() );
			else
				pmax1.x( Rmax.x() );
			if( Tmax.y() > Rmax.y() )
				pmax1.y( Tmax.y() );
			else
				pmax1.y( Rmax.y() );
			if( Tmax.z() > Rmax.z() )
				pmax1.z( Tmax.z() );
			else
				pmax1.z( Rmax.z() );
		}

		if(filter == 0) // filtration timing test enabled
		{
			// Initial vol. resize (now, the motion should be contained...)
			volT = FOPS::resize(volT,pmin1,pmax1,true);
			if( verbose > 1 )
				cout << prog_name << "> Re-sized Initial Map (volT) to fit motion " << volT->dim() << endl;

			// INITIAL VOLUME FOURIER PADDING
			pad1 = vlDim( (int)ceil( ( volT->dim().x()+resolution/volT->units().x()+2) / pad_fourier),
					(int)ceil( ( volT->dim().y()+resolution/volT->units().y()+2) / pad_fourier),
					(int)ceil( ( volT->dim().z()+resolution/volT->units().z()+2) / pad_fourier) );
			// INITIAL VOLUME PADDING
			if(pad1.x()!=0 || pad1.y()!=0 || pad1.z()!=0) // if we need padding !!!
			{
				volT = FOPS::padVolume( volT , pad1 ); // vol. must be padded to kernel filtration
				if( verbose > 1 )
					cout << prog_name << "> Padded Initial (and temporal) Volume: " << volT->dim() << endl;
			}

			// Filtration test (timing)
			vol_dummy = new vlVolume( volT );
			ht_timer.startTimer();
			for(int i=0; i < repeat_times; i++) // test is repeated "repeat_times" times
			{
#ifdef USE_PTHREAD // Enables PThread parallel routines
				project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,1,resolution,fast_switch,nthreads,threads_data); // Fourier filtration
#else
				project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,1,resolution,fast_switch); // Fourier filtration
#endif
				// fprintf(stderr,"volT avg= %f\n",FOPS::calc_average(volT));
			}
			delete(volT); // fillVolumeNE_3BB2R allocates memory...

			ht_timer.stopTimer();
			t_test1 = (float) ht_timer.getElapsedTime();
			delete(vol_dummy); // it will be allocated below (with the definitive size!)

			//				if( verbose > 0 )
			printf("%s> Fourier Filtration timing x%d: %fs  x1): %fs\n",prog_name,repeat_times,t_test1,t_test1/repeat_times);
			// End testing Fourier Filtration
		}

		if(filter == 0 || filter == 2) // filtration timing test enabled
		{
			// **************
			//  KERNEL TEST
			// **************
			if( verbose > 1 && filter == 0)
				cout << prog_name << "> TESTING KERNEL FILTRATION METHOD" << endl;

			// Projecting initial PDB
			volT = molt->fillVolumeNE_3BB2R( volR->units().x() );
			if( verbose > 1 )
				cout << prog_name << "> Initial map Box-size " << volT->dim() << endl;

			// INITIAL VOLUME KERNEL PADDING
			FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,volT->units().x(),resolution, sigma_kernel); // "3" sigma_factor
			pad2 = vlDim( (dim_vox-1)/2, (dim_vox-1)/2, (dim_vox-1)/2 );

			// INITIAL VOLUME PADDING
			if(pad2.x()!=0 || pad2.y()!=0 || pad2.z()!=0) // if we need padding !!!
			{
				volT = FOPS::padVolume( volT , pad2 ); // vol. must be padded to kernel filtration
				if( verbose > 1 )
					cout << prog_name << "> Padded Initial (and temporal) Volume: " << volT->dim() << endl;
			}

			// INITIAL VOLUME PROJECTION AND FILTRATION
			vol_dummy = new vlVolume( volT );
#ifdef USE_PTHREAD // Enables PThread parallel routines
			//				project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,2,resolution,fast_switch,nthreads,threads_data,threads); // Kernel filtration
			project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,2,resolution,fast_switch,nthreads,threads_data); // Kernel filtration
#else
			// FOPS::writeFile(volT,"volT_before.sit");
			project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,2,resolution,fast_switch); // Kernel filtration
#endif
			// FOPS::writeFile(volT,"volT_after.sit");
			// FOPS::writeFile(vol_dummy,"vol_dummy.sit");
			delete(vol_dummy); // it will be allocated below (with the definitive size!)
			free(kernel);

			// fprintf(stderr,"before crop\n");
			// SIZE PADDING - The whole movement must be contained !!! (due to sigmas...?)
			volT = FOPS::crop( volT, model_thr, true); // to see the proper size, it must be filtered!
			// fprintf(stderr,"after crop\n");
			if( verbose > 1 )
				cout << prog_name << "> Cropped Initial (and temporal) Volume size " << volT->dim() << endl;

			// THE WHOLE MOTION MUST BE CONTAINED INSIDE VOLUMES
			// Setting new map sizes
			Rdim = volR->dim(); // Map sizes
			Tdim = volT->dim();
			Rvox = volR->units(); // Samplings (voxel size)
			Tvox = volT->units();
			volR->getPosition(&Rorig); // Origins
			volT->getPosition(&Torig);
			Rmax.x( (float)Rdim.x() * (float)Rvox.x() + Rorig.x() );
			Rmax.y( (float)Rdim.y() * (float)Rvox.y() + Rorig.y() );
			Rmax.z( (float)Rdim.z() * (float)Rvox.z() + Rorig.z() );
			Tmax.x( (float)Tdim.x() * (float)Tvox.x() + Torig.x() );
			Tmax.y( (float)Tdim.y() * (float)Tvox.y() + Torig.y() );
			Tmax.z( (float)Tdim.z() * (float)Tvox.z() + Torig.z() );
			// Origin (min corner)
			if( Torig.x() < Rorig.x() )
				pmin2.x( Torig.x() );
			else
				pmin2.x( Rorig.x() );
			if( Torig.y() < Rorig.y() )
				pmin2.y( Torig.y() );
			else
				pmin2.y( Rorig.y() );
			if( Torig.z() < Rorig.z() )
				pmin2.z( Torig.z() );
			else
				pmin2.z( Rorig.z() );
			// Max corner
			if( Tmax.x() > Rmax.x() )
				pmax2.x( Tmax.x() );
			else
				pmax2.x( Rmax.x() );
			if( Tmax.y() > Rmax.y() )
				pmax2.y( Tmax.y() );
			else
				pmax2.y( Rmax.y() );
			if( Tmax.z() > Rmax.z() )
				pmax2.z( Tmax.z() );
			else
				pmax2.z( Rmax.z() );

			// Initial vol. resize
			volT = FOPS::resize(volT,pmin2,pmax2,true);
			if( verbose > 1 )
				cout << prog_name << "> Re-sized Initial Map (volT) to fit motion " << volT->dim() << endl;

			// INITIAL VOLUME KERNEL PADDING
			FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,volT->units().x(),resolution, sigma_kernel); // "3" sigma_factor
			pad2 = vlDim( (dim_vox-1)/2, (dim_vox-1)/2, (dim_vox-1)/2 );

			// INITIAL VOLUME PADDING
			if(pad2.x()!=0 || pad2.y()!=0 || pad2.z()!=0) // if we need padding !!!
			{
				volT = FOPS::padVolume( volT , pad2 ); // vol. must be padded to kernel filtration
				if( verbose > 1 )
					cout << prog_name << "> Padded Initial (and temporal) Volume: " << volT->dim() << endl;
			}

			// Initialization routine for "convoluteK_nopad_par"
#ifdef USE_PTHREAD // Enables PThread parallel routines
			if(nthreads > 0)
			{
				fprintf(stderr,"%s> Initializing convoluteK_nopad_par_init... nthreads= %d\n",prog_name,nthreads);
				FOPS::convoluteK_nopad_par_init(nthreads, &threads_data, &threads, kernel, dim_vox);
				//					fprintf(stderr,"Initializing project_RealNE_3BB2R_par_init... nthreads= %d\n",nthreads);

				//					fprintf(stderr,"before project_RealNE_3BB2R_par_init...\n");
				//					molt->project_RealNE_3BB2R_par_init(nthreads, &threads_data_proj, &threads_proj);
				//					fprintf(stderr,"after project_RealNE_3BB2R_par_init...\n");
			}
#endif

		}

		//fprintf(stderr,"KERNEL 1\n");
		//for(int x=0;x<pow(dim_vox,3); x++)
		//	fprintf(stderr,"%f ",kernel[x]);
		//fprintf(stderr,"\n");

		if(filter == 0) // filtration timing test enabled
		{
			// Mon made (27/5/2013)

			// Filtration test (timing)
			vol_dummy = new vlVolume( volT );
			ht_timer.startTimer();
			for(int i=0; i < repeat_times; i++) // test is repeated "repeat_times" times
				//				for(int i=0; i < 2000; i++) // test is repeated "repeat_times" times
			{
				//					fprintf(stderr,"Filtration iter: i= %d\n",i);
				//					fprintf(stderr,"before project_pdb...\n");
#ifdef USE_PTHREAD // Enables PThread parallel routines
				project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,2,resolution,fast_switch,nthreads,threads_data); // Fourier filtration
#else
				project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,2,resolution,fast_switch); // Fourier filtration
#endif
				//					fprintf(stderr,"after project_pdb...\n");
				//					sleep(1);
				// fprintf(stderr,"volT avg= %f\n",FOPS::calc_average(volT));
				//					char temp[50];
				//					sprintf(temp,"volT_%d.sit",i);
				//					FOPS::writeFile(volT,temp);
			}
			ht_timer.stopTimer();
			t_test2 = (float) ht_timer.getElapsedTime();
			delete(vol_dummy); // it will be allocated below (with the definitive size!)
			//				if( verbose > 0 )
			printf("%s> Kernel Filtration timing x%d: %fs  x1): %fs\n",prog_name,repeat_times,t_test2,t_test2/repeat_times);
			//FOPS::writeFile(volT,"volT_1.sit");

			//FOPS::writeFile(volT,"volT.sit");
			//exit(0);

			// Enable fastest filtration method
			if(t_test1 < t_test2)
				filter = 1; // Fourier
			else
				filter = 2; // Kernel

			if( verbose > 0 )
				cout << prog_name << "> Re-sized Target Map (volR) " << volR->dim() << endl;
			printf("%s> Best filtration method: %d FT(x%d)=%.3fs Kernel(x%d)=%.3fs\n",prog_name,filter,repeat_times,t_test1,repeat_times,t_test2);
			fprintf(f_com,"#%s> Best filtration method: %d FT(x%d)=%.3fs Kernel(x%d)=%.3fs\n",prog_name,filter,repeat_times,t_test1,repeat_times,t_test2);
		}
		// Filtering test end

		// FILTRATION RELATED PADDING (TARGET VOL.)
		// (once the motion is contained within the map)
		switch(filter)
		{
		case 1:
			volR = FOPS::resize(volR,pmin1,pmax1,true);
			// TARGET VOL. FOURIER PADDING
			pad = vlDim( (int)( ( volR->dim().x()+resolution/volR->units().x()+2) / pad_fourier),
					(int)( ( volR->dim().y()+resolution/volR->units().y()+2) / pad_fourier),
					(int)( ( volR->dim().z()+resolution/volR->units().z()+2) / pad_fourier) );
			break;
		case 2: // REMOVE THIS LATER...
			volR = FOPS::resize(volR,pmin2,pmax2,true);
			// Computing KERNEL (FILTERING)
			FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,volR->units().x(),resolution, sigma_kernel); // "3" sigma_factor
			if( verbose > 1 )
				printf("%s> Kernel size: dim_vox= %d\n",prog_name,dim_vox);
			// TARGET VOL. KERNELL PADDING
			pad = vlDim( (dim_vox-1)/2, (dim_vox-1)/2, (dim_vox-1)/2 );
			break;
		default:
			printf("Please, introduce a valid filtration method!\n");
			exit(1);
			break;
		}

		//fprintf(stderr,"KERNEL 2\n");
		//for(int x=0;x<pow(dim_vox,3); x++)
		//	fprintf(stderr,"%f ",kernel[x]);
		//fprintf(stderr,"\n");

		// PADDING REFERENCE, IF NECESSARY...
		if(pad.x()!=0 || pad.y()!=0 || pad.z()!=0) // if we need padding !!!
		{
			volR = FOPS::padVolume( volR, pad, true );
			if( verbose > 1 )
				cout << prog_name << "> Padded Target Volume (volR) size " << volR->dim() << endl;
		}

		// THE INITIAL AND MODEL VOLUMES WILL HAVE THE SAME SIZE AS TARGET ONE!
		// Volumes memory allocation (we wont need to allocate them any more!)
		volT = new vlVolume( volR );
		vol_dummy = new vlVolume( volR );
		volF = new vlVolume( volR );

		// INITIAL VOLUME PROJECTION AND FILTRATION
		// (Note that the current volume will be built from the the Reduced 3BB2R model !!!)
		if( verbose > 1 )
			printf( "%s> Projecting INITIAL PDB ( %s ) into Volume\n",prog_name,file_initial);
#ifdef USE_PTHREAD // Enables PThread parallel routines
		project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,filter,resolution,fast_switch,nthreads,threads_data);
#else
		project_pdb(molt,&volT,&vol_dummy,kernel,dim_vox,filter,resolution,fast_switch);
#endif
		if( verbose > 1 )
			printf("%s> Initial Volume Grid Size %f\n",prog_name,volT->units().x());
		//FOPS::writeFile(volT,"volT_2.sit");

		// If there are some negative values, cross-correlation may be negative... to avoid this:
		FOPS::threshold( volT, model_thr );

		// Pre-computing Target (reference) volume average and sigma
		// (to save time in cross-correlation computation)
		avgR = FOPS::avg_frameD(volR,pad);
		sigR = FOPS::sig_frameD(volR,pad);
		if( verbose > 0 )
			printf("%s> Precomputed Target Volume average= %f  and sigma = %f\n",prog_name,avgR,sigR);

		// INITIAL NORMALIZATION
		if(norm_switch)
		{
			norm_factor = FOPS::sigma(volT); // Computing normalization factor (sigma)
			if( verbose > 1 )
				printf( "%s> Computing Initial Model's Normalization-Factor = %f\n",prog_name,norm_factor);
			norm_vol(volR,1.0); // Normalization
			norm_vol(volT,1.0,norm_factor);
		}

		if(savemaps_switch)
		{
			// Writing the used target volume (it was not filtered)
			sprintf(text,"%s_target.sit",name);
			FOPS::writeFile(volR,text);
			if( verbose > 1 )
				printf("%s> Saved Target (reference) Volume: %s\n",prog_name,text);
			// Writing filtered initial volume
			sprintf(text,"%s_initial.sit",name);
			FOPS::writeFile(volT,text);
			if( verbose > 1 )
				printf("%s> Saved Initial (model) Volume: %s\n",prog_name,text);
		}
	}
	//	}
	fclose(f_com);

	// Initializing some residue variables (# atoms, # dihedrals, etc...)
	int size=0,i; // hessian rank
	tri *props;
	int *unat;

	// DEFINING PROPERTIES
	if( verbose > 1 )
		printf ("%s>  3) Initializing some CG-model related variables \n",prog_name);
	switch(model)
	{
	case 0:
	case 3:
		properCA(molt,&props,&unat);
		break;
	case 1:
	case 2:
		properMFA(molt,&props,&unat,type,model);
		break;
	}

	// Theoretic number of DoFs= 3*T + 3*R -6 + Dihedrals
	// (Before fixing, i.e. the fix-file format DoFs...)
	int n_seg,n_chain,old_size,seg_atoms=0;
	pdbIter *iter = new pdbIter( molt, true, true, true, true );
	size = 0;
	old_size = 0; // temp (Dihedral ICs)
	for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		size += props[iter->pos_fragment].nan; // Each residue may have different number of dihedral-angles!
	printf( "%s> Number of Dihedral angles: %d\n",prog_name,size );
	old_size = size;
	for( iter->pos_segment = 0; !iter->gend_segment(); iter->next_segment() ) // screen segments
	{
		seg_atoms = (iter->get_segment())->num_atoms();
		if(seg_atoms > 1)
			size += 6;
		else
			size += 3;
	}
	printf( "%s> Rotational/Translational ICs (Non-Eckart): %d\n",prog_name,size-old_size );
	size -= 6; // Eckart conditions
	printf( "%s> Predicted number of ICs (Eckart): %d\n",prog_name,size );
	n_seg = iter->num_segment();
	n_chain = iter->num_chain(); // WARNING! With SMOL's unknown behavior!!!
	delete iter;

	// FIXATION of VARIABLES
	// Mobile Internal Coordinates selection "i"-model
	// Fixation array is related to "i"-model, not to "current" one!
	// input: fixedi ---> output: fixedx
	bool *fixed=NULL;
	if(fixmodel != 0)
	{
		fixed = (bool *)malloc(sizeof(bool)*size);
		check_pointer(fixed,"Fixation array");
		old_size = size;
		switch(fixmodel)
		{
		case 2:
			size = read_fixIC(fix_file,molt,props,fixed);
			break;
		case 3:
			size = read_fixDH(fix_file,molt,props,fixed,type,model);
			break;
		case 4:
			size = fixRand(fixed,old_size,fixRand_prob);
			break;
		case 5:
			size = fixRandDH(molt,props,fixed,type,model,fixRand_prob);
			break;
		case 6:
			if(!ss_switch) // if not specified a SS file, it will be computed
				ss_table = molt->secondary_structure(true);
			size = fixSS(molt,props,fixed,ss_table,fix_ss,type,model);
			break;
		default:
			printf("Sorry, unknown fixation method!\nForcing exit\n");
			exit(1);
			break;
		}
		if(fullatom_switch) // UPDATE "fixed0" (full-atom mask) from "fixed" (current model mask)
			changefixIC(mol0, props0, props, fixed0, fixed, 2, type, model, type, true);
		printf("%s> Input CG-model Fixed Internal Coordinates: %d\n",prog_name, old_size-size );
		printf("%s> Input CG-model Mobile Internal Coordinates (size) = %d\n",prog_name,size);
	}
	else if(fullatom_switch) // Needed because CG-models-modes have different "size" than full-atom!
	{
		fixed = (bool *)malloc(sizeof(bool)*size);
		check_pointer(fixed,"Fixation array");
		for(int i=0; i<size; i++)
			fixed[i]=true; // all mobile
		changefixIC(mol0, props0, props, fixed0, fixed, 2, type, model, type, true);
	}

	// Creates two auxiliar arrays with segment properties (needed due to fixing):
	//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
	//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
	// (Allocates memory itself, if it's needed)
	bool *addrot=NULL; // Should be added ROTATIONs? (with fixing)
	int *effseg=NULL; // Effective segment indices
	size = seg_props(molt, props, fixed, model, type, &addrot, &effseg);
	printf( "%s> Number of ICs predicted by seg_props(): %d\n",prog_name,size );

	// Number of USED eigenvectors (defined here because we need to know "size" first)
	if(nevec_fact >= 1.0) // number of modes
		nevs = (int) nevec_fact;
	else
		nevs = (int) (nevec_fact * size);
	nevs += addnevs; // Increases the nevs value by addnevs to account for small structures with low DoFs number

	// Checking
	if(nevs > size)
	{
		printf("%s> Sorry, more --nevs requested (%d) than available (%d), forcing the maximum.\n",prog_name,nevs,size);
		nevs = size;
	}
	else if(nevs <= 0) // checking
	{
		printf("%s> Error, invalid --nevs requested %d (%f)!\nForcing exit!\n",prog_name,nevs,nevec_fact);
		exit(1);
	}
	nevs_ratio = (float)nevs/size; // to select diagonalization routine
	printf( "%s> Range of used modes: 1-%d (%.1f%%)\n",prog_name,nevs,nevs_ratio*100);

	if(eigensolver_switch) // Automatic eigensolver selection overridden
	{
		printf("%s> Using the user-selected eigensolver: %d\n",prog_name,eigensolver);
	}
	else  // Automatic eigensolver selection
	{
		if(size<100) // For less than 100 DoFs (small molecules) BLAS/LAPACK's eigensolver (option 0) will be selected regardless nevs_ratio
		{            // (because there must be some kind of bug in the ARPACK-based routine for small molecules)
			eigensolver = 0;
			printf("%s> Small numeber of DoFs, using 0 eigensolver!\n",prog_name);
		}
		else
		{
			if(nevs_ratio > 0.06) // tradeoff threshold is around 6% of the modes
				eigensolver = 0;
			else
				eigensolver = 1;
			printf("%s> Using the fastest eigensolver for %.1f%% eigenvectors: %d\n",prog_name,nevs_ratio*100,eigensolver);
		}
	}

	// Number of 2nd EXCITED eigenvectors (we need to know "nevs" first)
	if(nex_fact2 >= 1.0) // number of modes (integer), (<= 0 was checked before...)
		nex2 = (int) nex_fact2;
	else // ratio from available
		nex2 = (int) (nex_fact2 * nevs);

	// Number of EXCITED eigenvectors (we need to know "nevs" first)
	if(nex_fact >= 1.0) // number of modes (integer)
		nex = (int) nex_fact;
	else // ratio from available
		nex = (int) (nex_fact * nevs);

	if(nex==nex2)
		printf( "%s> Number of excited/selected modes: %d(nex)\n",prog_name, nex);
	else
		printf( "%s> Number of excited/selected modes: %d(nex)  final: %d(nex2)\n",prog_name, nex,nex2);

	// Checking
	if(nex > nevs)
	{
		printf("%s> Sorry, more excitation eigenvectors requested (%d) than available (%d), forcing maximum.\n",prog_name,nex,nevs);
		nex = nevs;
	}
	else if(nex <= 0) // checking
	{
		printf("%s> Warning, invalid number of excitation requested %d (%.2f), forcing nex=1\n",prog_name,nex,nex_fact);
		nex = 1;
	}

	//	// Number of 2nd EXCITED eigenvectors (we need to know "nevs" first)
	//    if(nex_fact2 >= 1.0) // number of modes (integer)
	//    {
	//    	nex2 = (int) nex_fact2;
	//		printf( "%s> Number of requested 2nd excited modes: %d\n",prog_name, nex2);
	//    }
	//    else // ratio from available
	//    {
	//    	nex2 = (int) (nex_fact2 * nevs);
	//    	printf( "%s> Number of requested 2nd excited modes: %d (%.1f%% from available)\n",prog_name, nex2, nex_fact2*100);
	//    }

	// Checking
	if(nex2 > nevs)
	{
		printf("%s> Warning, invalid --nex %d (from %d). Setting nex=%d\n",prog_name,nex2,nevs,nevs);
		nex2 = nevs;
	}
	else if(nex2 <= 0) // checking (could be 0 if nevs very low)
	{
		printf("%s> Warning, invalid --nex2 %d (%.2f). Setting nex2=%d\n",prog_name,nex2,nex_fact2,nex2);
		nex2 = 1;
	}

	int *chain_order; // stores chain order, during 6D pose refinement
	bool chain_exist=false; // bool to avoid chain repetition
	chain_order = (int *) malloc( sizeof(int) * n_chain);
	check_pointer(chain_order,"Chain order array");
	for(int i=0;i<n_chain;i++)
		chain_order[i] = -1; // initialization

	if( verbose > 0 )
		printf("%s> Using %d (%5.2f%%) lowest energy eigenvectors (%d)\n",prog_name,nevs,((float)nevs/size)*100,size);

	if( verbose > 1 )
		printf ("%s>  4) Allocating & Initializing memory( hess_matrix, cont_matrix, eigval)\n",prog_name);
	double *eigval=NULL, *eigval2=NULL, *hess_matrix=NULL, *hess_matrix2=NULL;

	double *evec;
	double *profile,*profile_soft; // "corr" method (from "prob_method")
	double *prof_acorr,*prof_acorr_soft;
	double *prof_corr_soft;
	double *prof_wrmsd=NULL;
	double *prof_corr=NULL;

	Macromolecule *forward; // if model==0, it's a selection from "forward2" (select_cg_CA)
	Macromolecule *forward2;
	Macromolecule *forwardT; // if model==0, it's a selection from "forwardT2" (select_cg_CA)
	Macromolecule *forwardT2;
	Macromolecule *p_fwd,*p_fwdT,*p_bcw,*p_bcwT,*molt_CA,*molt_CA0,*molr_CA,*molt_bb,*molr_bb,*molt_rediag;
	double score,score_CA=0,score_bb=0,score_f,score_old,score_aligned,first_score;
	int i_atom;

	// Opening Score file
	FILE *f_score;
	sprintf(file_score,"%s_score.txt",name);
	if( !(f_score=(FILE *)fopen(file_score,"w") ) )
	{
		printf("Sorry, unable to open output score file: %s\n",file_score);
		exit(1);
	}

	if(pdb_switch)
		fprintf(f_score,"# %s> Score file\n#%5s %11s %11s\n",prog_name,"Iter","Score","RMSD_CA");
	else
	{
		if(more_rmsds)
			fprintf(f_score,"# %s> Score file\n#%5s %11s %11s %11s %11s %11s\n",prog_name,"Iter","Score","RMSD_CA","cross.corr.","Delta_Score","D_Score_Avg");
		else
			fprintf(f_score,"# %s> Score file\n#%5s %11s %11s %11s %11s\n",prog_name,"Iter","Score","cross.corr.","Delta_Score","D_Score_Avg");
	}

	// Mon: Reading alignment mask file here?
	bool *maskres1,*maskres2,*maskmolt,*maskmolr,*maskmoltCA,*maskmolrCA;
	int nmatch,nali1,nali2;

	if(aln_switch)
	{
		printf( "%s> Reading Clustal's sequence alignment from: %s\n",prog_name,file_aln);
		read_aln(file_aln, &maskres1, &nali1, &maskres2, &nali2, &nmatch, aln_level);
		printf( "%s> %d residues matched!  nali1= %d  nali2= %d\n",prog_name,nmatch,nali1,nali2);

		maskmolt = molt->maskres2maskatom(maskres1);
		maskmolr = molr->maskres2maskatom(maskres2);


		// molt->writePDB("maskmolt.pdb");
		// molr->writePDB("maskmolr.pdb");

		// Counting number of non-hetero residues
		pdbIter *iter_seg;
		TMOL fragtype;
		Segment *seg;
		iter_seg = new pdbIter(molt, true, true, true, false );

		int nres1 = molt->get_num_fragments();
		for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
		{
			seg = ( Segment * ) iter_seg->get_segment();
			fragtype = seg->getMolType();
			if(fragtype == tmol_smol) // if hetero
				nres1--; // discount the number of not-hetero residues
		}
		iter_seg->clean_virtual();
		delete iter_seg;
		// Some checking...
		if(nali1 != nres1)
		{
			printf("%s>\n%s> ERROR: Number of residues mismatch between Initial PDB (%d) and the corresponding sequence (%d) in: %s"
					"\n%s> Forcing exit!\n%s>\n",prog_name,prog_name,nres1,nali1,file_aln,prog_name,prog_name);
			//exit(1);
		}

		// Counting number of non-hetero residues
		iter_seg = new pdbIter(molr, false, true, true, false );
		int nres2 = molr->get_num_fragments();
		for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
		{
			seg = ( Segment * ) iter_seg->get_segment();
			fragtype = seg->getMolType();
			if(fragtype == tmol_smol) // if hetero
				nres2--; // discount the number of not-hetero residues
		}
		iter_seg->clean_virtual();
		delete iter_seg;
		// Some checking...
		if(nali2 != nres2)
		{
			printf("%s>\n%s> ERROR: Number of residues mismatch between Target PDB (%d) and the corresponding sequence (%d) in: %s"
					"\n%s> Forcing exit!\n%s>\n",prog_name,prog_name,nres2,nali2,file_aln,prog_name,prog_name);
			// exit(1);
		}

		// THIS HAD TO BE REMOVED IN ORDER TO ALLOW LIGANDS (HETERO-ATOMS) MORPHING/FITTING
		//
		//		// Some checking...
		//		if(nali1 != num_res)
		//		{
		//			printf("%s>\n%s> ERROR: Number of residues mismatch between Initial PDB and the corresponding sequence in: %s"
		//					"\n%s> Forcing exit!\n%s>\n",prog_name,prog_name,file_aln,prog_name,prog_name);
		//			exit(1);
		//		}
		//		if(nali2 != num_resR)
		//		{
		//			printf("%s>\n%s> ERROR: Number of residues mismatch between Target PDB and the corresponding sequence in: %s"
		//					"\n%s> Forcing exit!\n%s>\n",prog_name,prog_name,file_aln,prog_name,prog_name);
		//			exit(1);
		//		}
		//
		// Testing masked RMSD
		//		printf("Masked RMSD: %f\n",molt->rmsd(molr,maskmolt,maskmolr));
		//		exit(0);
	}

	// COMPUTING INITIAL SCORES
	if( pdb_switch || bench_switch ) // PDB-PDB Morphing || Benchmark Fitting
	{
		molr_CA = molr->select(calpha2); // Needed for last rmsds
		//		molr_CA->writePDB("molr_CA.pdb");
		if(aln_switch)
			maskmolrCA = molr_CA->maskres2maskatom(maskres2); // translates residue level mask into atomic
		//		molr_CA->writePDB("maskmolrCA.pdb");
		if(more_rmsds)
		{
			//molt->writePDB("molt.pdb");
			molt_CA = molt->select_cpy(calpha2);
			//molt->writePDB("molt2.pdb");

			if(aln_switch)
			{
				// molt_CA->writePDB("molt_CA.pdb");
				// molr_CA->writePDB("molr_CA.pdb");

				sprintf(text,"%s_initial.rmsd",name);
				fprintf(stdout,"imorph> Saving %s_initial.rmsd\n", name);

				maskmoltCA = molt_CA->maskres2maskatom(maskres1); // translates residue level mask into atomic
				score_CA = molr_CA->rmsd_file(molt_CA,maskmolrCA,maskmoltCA,text); // first score CA's (RMSD)
				score_CA = molr_CA->rmsd(molt_CA,maskmolrCA,maskmoltCA); // first score CA's (RMSD)
				// output dihedrals && SS
				//
				float *di;
				int *ss, dangi;
				Fragment * res, *res2;
				int resn, Nres, simple;
				Segment * seg;
				int chino;
				float * chis;
				Chain *ch;
				pdbIter *iter1, *iter2;
				iter1 = new pdbIter( molrini,  false, false, true, false); // r maskres2
				iter2 = new pdbIter( moltini,  false, false, true, false ); // t maskres1
				// first SS assign

				molrini->all_dihedrals( &di);

				int ntemp;
				ntemp=iter1->num_fragment();
				// printf ("ntemp %d \n", ntemp);
				dihedrals2DISICL( di, &ss, ntemp );


				FILE *f_out=NULL;
				sprintf(text,"%s_reference.dih",name);
				fprintf(stdout,"imorph> Saving %s_reference.dih\n", name);
				if( !(f_out = fopen(text,"w")) )
				{
					fprintf(stderr,"Sorry, unable to write output file! Forcing exit!\n");
					exit(2);
				}
				fprintf(f_out, "# AA     resn    PHI      PSI     OMEGA     DISICL  simple    chi1     chi2     chi3      chi4\n");
				iter2->pos_fragment = 0;	dangi=0; Nres=0;
				for ( iter1->pos_fragment = 0; !iter1->gend_fragment(); iter1->next_fragment() )
				{
					res = ( Fragment * ) iter1->get_fragment();
					if(maskres2[iter1->pos_fragment]) {
						while(!iter2->gend_fragment() && !maskres1[iter2->pos_fragment]) // this places the index into the corresponding pas
							iter2->next_fragment();
						//						if(maskres1[iter2->pos_fragment]) {
						simple=DISICL_d2simple[ss[Nres]];  // change DISICL simple SS classes
						fprintf( f_out, "%3s %8d  %8.3f %8.3f %8.3f  %3d %3s %3d %3s  %8.3f %8.3f %8.3f %8.3f\n",
								res->getName(), res->getIdNumber(),
								di[dangi], di[dangi+1],  di[dangi+2], // (phi,psi,omega)
								ss[Nres],DISICL_d[ss[Nres]],	simple, DISICL_s[simple],
								di[dangi+3], di[dangi+4], di[dangi+5], di[dangi+6]  // (4x Chi)
						);
						iter2->next_fragment();
						//						}
					}
					dangi+=7;
					Nres++;
				}
				iter1->~pdbIter();
				iter2->~pdbIter();
				fclose(f_out);





				// Testing masked RMSD
				//				printf("Masked RMSD_CA: %f\n",molt_CA->rmsd(molr_CA,maskmoltCA,maskmolrCA));
			}
			else
				score_CA = molr_CA->rmsd(molt_CA); // first score CA's (RMSD)
			// pending rmsd....PABLO 2019
		}
	}
	//exit(0);

	// INITIAL SCORES
#ifdef GAUSS
	GAUSS3D *gmolr_CA,*gmolt_CA;
	double gauss_factor,gauss_factorR,gauss_factorT;
#endif

	int n_diag=0;
	float nmatime;
	printf("%s>\n",prog_name);
	if( pdb_switch ) // PDB-PDB Morphing
	{
		//		score = molr->rmsd(molt); // first score (RMSD)

		// GAUSSIAN MIXTURE MODEL RELATED STUFF... (Future work)
		if(gauss_switch)
		{
#ifdef GAUSS
			//			gauss_factorR = Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolr,num_atomsR,num_atomsR);
			gauss_factorR = Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolr,num_gauss,num_gauss);
			gauss_factorT = Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolt,gmolt,num_atoms,num_atoms);
			if(gauss_factorR > gauss_factorT)
				gauss_factor = gauss_factorR;
			else
				gauss_factor = gauss_factorT;
			//			score = Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_atomsR,num_atoms,gauss_factor);
			score = Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_gauss,num_atoms,gauss_factor);
#endif
		}
		else
			if(aln_switch)
				score = molr->rmsd(molt,maskmolr,maskmolt); // first score (RMSD)
			else
				score = molr->rmsd(molt); // first score (RMSD)

		first_corr = score;
		if( verbose > -1 )
		{
			if(more_rmsds)
				printf("%s> %5s %9s %9s %3s %10s\n",prog_name,"Iter","RMSD","RMSD(CA)","NMA","NMA_time");
			else
				printf("%s> %5s %9s %3s %10s\n",prog_name,"Iter","RMSD","NMA","NMA_time");
			if(!server)
			{
				if(more_rmsds)
					printf("%s> %5d %9.6f %9.6f",prog_name,0,score,score_CA);
				else
					printf("%s> %5d %9.6f",prog_name,0,score);
			}
		}
		if(!server)
		{
			if(more_rmsds)
				fprintf(f_score,"%6d %11.4e %11.4e\n",0,score,score_CA);
			else
				fprintf(f_score,"%6d %11.4e\n",0,score);
		}
	}
	else // PDB-Volume
	{
		// Initial score
		score = 1 - FOPS::correlation_frame(volT,volR,pad,corr_nozero,false,avgR,sigR);
		first_corr = 1-score;

		// First output lines
		if( verbose > -1 )
			if(bench_switch && more_rmsds)
			{
				printf("%s> %5s %9s %9s %9s %3s %10s\n",prog_name,"Iter","Score","RMSD(CA)","Corr.","NMA","NMA_time");
				printf("%s> %5d %9.6f %9.6f %9.6f",prog_name,0,score,score_CA,1-score);
			}
			else
			{
				printf("%s> %5s %9s %9s %3s %10s\n",prog_name,"Iter","Score","Corr.","NMA","NMA_time");
				printf("%s> %5d %9.6f %9.6f",prog_name,0,score,1-score);
			}
		if(bench_switch && more_rmsds)
			fprintf(f_score,"%6d %11.4e %11.4e %11.4e\n",0,score,score_CA,1-score);
		else
			fprintf(f_score,"%6d %11.4e %11.4e\n",0,score,1-score);
	}

	// Storing the initial CA-model
	if(deltasave_rmsd_switch)
	{
		//		if(more_rmsds)

		molt_CA0 = molt->select_cpy(calpha2); // selecting again because pointers were destroyed
		molt_CA = molt->select_cpy(calpha2); // selecting again because pointers were destroyed


		//		else
		//		{
		//			molt_CA0 = new Macromolecule(molt); // Initial reference structure for "deltasave_rmsd"
		//			molt_CA = new Macromolecule(molt); // Initial reference structure for "deltasave_rmsd"
		//		}
	}

	molt_rediag = new Macromolecule(molt); // Last diagonalization macromolecule (current CG-model)

	if(score > 0)
		score_old=1e12; // a very big initial score...(minimization)
	score_lastsaved = score_old;
	first_score = abs(score); // first score buffer

	// Initializing Normal Mode excitation probabilities: "PLAIN" METHOD
	relev = (float *) malloc( sizeof(float) * nevs ); // Normal Mode Relevances (or probabilities)
	check_pointer(relev,"Normal mode probabilities (relevances)");
	switch( prob_method )
	{
	case 1: // PLAIN
		if( verbose > 0 )
			printf("%s> Initializing Normal Mode excitation probabilities: \"plain\" METHOD\n",prog_name);
		for( int i=0; i<nevs; i++)
			relev[i] = 1;
		break;

	case 4: // LINE
		if( verbose > 0 )
			printf("%s> Initializing Normal Mode excitation probabilities: \"line\" METHOD\n",prog_name);
		for( int i=0; i<nevs; i++)
			relev[i] = 1.0 - ((float)1/nevs)*i;
		break;

		//		// UN-TESTED
		//	case 5: // Gaussian
		//		if( verbose > 0 )
		//			printf("%s> Initializing Normal Mode excitation probabilities: \"gaussian\" METHOD\n",prog_name);
		//		for( int i=0; i<nevs; i++)
		//			relev[i] = (double) exp( -pow((i-(nevs * p_x0))/((double) nevs * p_s * sqrt((double)2)),2)/2);
		//		break;
		//
		//	// UN-TESTED
		//	case 6: // Inverse Exponential
		//		if( verbose > 0 )
		//			printf("%s> Initializing Normal Mode excitation probabilities: \"invexp\" METHOD\n",prog_name);
		//		for( int i=0; i<nevs; i++)
		//			relev[i] = (double) 1/( 1 + pow(((double) nevs * p_x0)/i,p_s) );
		//		break;
	}

	// Showing probabilities...
	if( verbose > 1 && prob_method != 2 )
		for( int i=0; i<nevs; i++)
			printf("%s> %2d %8.6f\n",prog_name,i+1,relev[i]);

	// Internal Coordinates motion vector
	double *uu=NULL; // Stores current Internal Coords. displacement
	double *UU=NULL; // Stores the absolute (cummulative) Internal Coords. displacement
	uu = (double *) malloc( sizeof(double) * size );
	check_pointer(uu,"Current Internal Coords. displacement");
	UU = (double *) malloc( sizeof(double) * size );
	check_pointer(UU,"Absolute (cummulative) Internal Coords. displacement");
	for(int i=0;i<size;i++)
		UU[i] = 0.0; // it should be initialized

	///////////////////////////// BEGIN TESTING /////////////////////////////////
	//	{
	//		int repeats = 10000;
	//		float test_time;
	//
	//		sprintf(dummy_string,"test0.pdb");
	//		molt2->writePDB( dummy_string, false ); // renumbers the PDB
	//		fprintf(stderr,"Written initial pdb: %s\n",dummy_string);
	//
	//		// Creating residue-level array of iterators
	//		pdbIter *iter,**iter_array;
	//		Residue *res;
	//
	//		iter_array = (pdbIter **) malloc( sizeof(pdbIter *) * molt2->get_num_fragments() );
	//		check_pointer(iter_array,"Residue-level array of iterators");
	//
	//		iter = new pdbIter( molt2 );
	//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
	//		{
	//			res = ( Residue * ) iter->get_fragment();
	//			iter_array[ iter->pos_fragment ] = new pdbIter( res ); // iter residue atoms
	//		}
	//
	////		rand_vector(uu, size, -0.1, 0.1); // Fills a vector with random values from "low" to "high" [low:high)
	//
	//		ht_timer.startTimer();
	//		fprintf(stderr,"\nTESTING\n");
	//		for(int k=0; k<repeats; k++)
	//		{
	//			rand_vector(uu, size, -0.1, 0.1); // Fills a vector with random values from "low" to "high" [low:high)
	//			//		for(int n=0;n<size;n++)
	//			//			fprintf(stderr,"%f ",uu[n]);
	//			//		fprintf(stderr,"\n");
	//			move_dihedralMCAx(molt2, uu, props, 1.0, type, model, fixed, iter_array);
	//			//		fprintf(stderr,"Repeating %d times\r",k+1);
	//		}
	//		test_time = (float) ht_timer.clocks() / CLOCKS_PER_SEC;
	//
	//		fprintf(stderr,"\nTest_timer (x%d) = %fs (%fms x1)\n",repeats,test_time,1000*test_time/repeats);
	//
	//		sprintf(dummy_string,"test.pdb");
	//		molt2->writePDB( dummy_string, false ); // renumbers the PDB
	//		fprintf(stderr,"Written random pdb: %s\n\n",dummy_string);
	//
	//		exit(0);
	//	}
	///////////////////////////// END TESTING /////////////////////////////////

	Macromolecule *molx,*molx_model;
	int movie_index = 1; // Movie index to be fully compatible with Multi-PDB format (and Chimera's MD movie module)
	if(movie_switch) // Deleting old movie
	{
		FILE *f_temp;
		sprintf(text,"%s_movie.pdb",name);
		f_temp = fopen(text,"w");
		fclose(f_temp);
		if(fullatom_switch)
			mol0->writeMPDB(text,movie_index);
		else
			molt->writeMPDB(text,movie_index);
		if( verbose > 1 )
			printf("%s> Added First model to %s (iter=%d)\n",prog_name,text,0);
		if( verbose > 0 )
			printf("%s> Movie created! ( %s )\n",prog_name, text);
	}

	if(chimera_switch) // Specific output for UCSF-Chimera plugin
	{
		// Save "_fitted.pdb" at the beginning...
		sprintf(text,"%s_fitted.pdb",name); // Selected atom-model
		if(fullatom_switch)
		{
			mol0->writePDB(text); // Write "fitted"
			// Write a full-atom initial model (only for Chimera)
			sprintf(text,"%s_model.pdb",name); // Write Full-atom model into "model"
			mol0->writePDB(text);
		}
		else
			molt->writePDB(text);
	}

	double mtot;
	double mta;
	double rd[3];
	// Triangular packed matrix input arrays // Hessian generation routines will alocate the memory!
	double *mass_tr = NULL;
	double *hess_tr = NULL;
	// Squared matrices instead of triangular...
	double *mass_sq = NULL;
	double *hess_sq = NULL;
	pdbIter *iter_atom;

	time_main.restart();
	// ******************************************************************************************
	// * MAIN LOOP
	// ******************************************************************************************
	int n_iter=1;
	bool apply_pose;
	bool rediag_next;
	int rediag_trials=0;
	float corr;
	pose_t = pose_t0; // setting initial pose refinement temperature

	// For aligning every chain independently (with fullatom_switch=true)
	Macromolecule *selm; // selected chain "scoring" model
	Macromolecule *selx; // selected chain "full-atom" model
	Macromolecule *selxm; // selected chain "moving" model

	// Allocating trial macromolecule/s only once!
	Macromolecule *buff; // swapping buffer
	if(model == 0) // CA-model
	{
		// NCAC-model
		forward2 = new Macromolecule(molt2); // forward movement PDB
		// CA-model (references to NCAC-model's atoms)
		forward = select_cg_CA(forward2); // Copying CA-model pointer tree from the NCAC-model.
	}
	else // C5- or HA- models
		forward = new Macromolecule(molt); // forward movement PDB

#ifdef GAUSS
	// GAUSSIAN MIXTURE MODEL RELATED STUFF... (Future work)
	// Make gaussian trial molecule
	if(gauss_switch)
		gforward = pdb2gauss(forward,gauss_sd);
#endif

	// Eigenvectors/values backup (in case diagonalization fails)
	if(eigval2!=NULL) // already allocated
		free(eigval2);
	eigval2 = (double *) malloc(size * sizeof(double));
	check_pointer(eigval2,"Eigenvalue backup memory");
	if(hess_matrix2!=NULL)
		free(hess_matrix2);
	hess_matrix2 = (double *) malloc(size*nevs * sizeof(double));
	check_pointer(hess_matrix2,"Eigenvectors backup memory");

	bool *maskmolt2,*maskmolr2;
	if(aln_switch && model == 0)
	{
		maskmolt2 = molt2->maskres2maskatom(maskres1);
		maskmolr2 = molr2->maskres2maskatom(maskres2);
	}


	// The MAIN LOOP !!!
	while( do_main_loop && n_iter < max_iter )
	{
		if( verbose > 1 )
			printf("# Iter: %3d   ###############################################################\n",n_iter);

		// PREVIOUS COMPUTATIONS
		if(time_switch)
			ht_stage.startTimer();
		do_refine = false;

		if( refine_each > 0 ) // Equi-distant refinements
		{
			if( n_iter%refine_each == 0 )
			{
				if( verbose > 0 )
					printf("%s> POSITION & ORIENTATION REFINEMENT SELECTED (iter=%d)\n",prog_name,n_iter);
				do_refine = true;
			}
		}
		else // Randomly separated refinements
		{
			dice = (double) rg->Random(); // [0:1) Playing dice with Mersenne!
			if( verbose > 1 )
				printf("%s> POSITION & ORIENTATION REFINEMENT CHANCE (iter=%d)\n",prog_name,n_iter);
			if( dice < refine_prob )
			{
				if( verbose > 1 )
					printf("%s> ACEPTED! (dice=%5.3f < %5.3f)\n",prog_name,dice,refine_prob);
				do_refine = true;
			}
			else
				if( verbose > 1 )
					printf("%s> REJECTED! (dice=%5.3f >= %5.3f)\n",prog_name,dice,refine_prob);
		}

		if(n_iter < refine_delay)
			do_refine = false; // cancels 6D pose refinement

		if(n_iter == refine_delay)
				do_refine = true; // cancels 6D pose refinement

		// Do 6D pose refinement
		if( do_refine )
		{

			if(pdb_switch) // Equivalent, but for PDB-PDB Morphing
			{

				// Align PDBs
				if( verbose > 1 )
					printf("%s> minRmsd - Aligning all!\n",prog_name);

				if(wrmsd_switch)
				{
					if(model==0) // If CA model
					{
						float score_0;

						if(aln_switch)
							score_0 = molt2->rmsd(molr2,maskmolt2,maskmolr); // first score (RMSD)
						else
							score_0 = molt2->rmsd(molr2); // score computed with selected model


						if(aln_switch)
							molr2->minWRmsd(molt2,matrix4,prof_wrmsd,maskmolr2,maskmolt2);
						else
							molr2->minWRmsd(molt2,matrix4,prof_wrmsd);

						molx = new Macromolecule(molt2); // molecule copy from the HA readed one
						matrix4_op = new M4Rot(matrix4);
						molx->applyAtoms(matrix4_op);

						float score_R;
						if(aln_switch) {
						score_R = molx->rmsd(molr2,maskmolt2,maskmolr2);
						} else {
						score_R = molx->rmsd(molr2); // score computed with selected model
						}
						delete molx;

						//fprintf(stderr,"   scoret0= %f %f %e\n",score_R, score_0, score);

						//if (score_R<score) {
							if (score_0-score_R>0.000001) //
							{
							molt2->applyAtoms(matrix4_op);

							if(gauss_switch)
							{
#ifdef GAUSS

								pdb2gauss_update(molt,gmolt);
								//								score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_atomsR,num_atoms,gauss_factor);
								score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_gauss,num_atoms,gauss_factor);
#endif
							}
							else {
								score_R =score;
								if(aln_switch)
									score = molt->rmsd(molr,maskmolt,maskmolr); // first score (RMSD)
								else
									score = molt->rmsd(molr); // score computed with selected model
							}
							// fprintf(stderr,"   scoret= %f %f %e\n",score_R, score, score-score_R);


						}
					}
					else // If C5 or HA models
					{
						if(aln_switch)
							molr->minRmsd(molt,matrix4,maskmolr,maskmolt);
						else
							molr->minRmsd(molt,matrix4);

						molx = new Macromolecule(molt); // molecule copy
						matrix4_op = new M4Rot(matrix4); // rotate
						molx->applyAtoms(matrix4_op);

						float score_R;
						if(aln_switch) {
						score_R = molx->rmsd(molr,maskmolt,maskmolr);
						} else {
						score_R = molx->rmsd(molr); // score computed with selected model
						}
						delete molx;
						/*

						float matrix0[4][4];



						for(int i=0; i < 4; i++)
							for(int j=0; j < 4; j++)
							    matrix0[i][j]=matrix4[i][j];


					if ((score-score_R) > delta_save/2) {

							for(int rn=1; rn <= 10; rn++) {

							for(int i=0; i < 3; i++)
								matrix4[i][3]=matrix0[i][3]*rn/10.0;

							//fprintf(stderr,"  %d scoret0= %f %f %f %f\n", rn, score_R, score, score-score_R, 0.1/(score-score_R));

							molx = new Macromolecule(molt); // molecule copy
							matrix4_op = new M4Rot(matrix4); // rotate
							molx->applyAtoms(matrix4_op);

							score_R = molx->rmsd(molr,maskmolt,maskmolr);
							delete molx;

							if ((score-score_R) > delta_save/2) break;

							// fprintf(stderr,"  %d scoret1= %f %f %f\n", rn, score_R, score, score-score_R);

							}

						}
*/

						/*
						for(int i=0; i < 4; i++)
						{
							for(int j=0; j < 4; j++)
								printf("%6.3f ", matrix4[i][j]);
							printf("\n");
						}
                        */


						if (score-score_R>0.00000001) {
							//fprintf(stderr,"   scoret2= %f %f %e\n",score_R, score, score-score_R);

							molt->applyAtoms(matrix4_op);
							//score_R = molt->rmsd(molr,maskmolt,maskmolr);
							//fprintf(stderr,"   scoret2= %f %f %e\n",score_R, score, score-score_R);

							score=score_R;
							if(gauss_switch)
							{
#ifdef GAUSS

								pdb2gauss_update(molt,gmolt);
								score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_gauss,num_atoms,gauss_factor);
#endif
							}
						}
					}



					//					if(aln_switch)
					//						molr->minWRmsd(molt,matrix4,prof_wrmsd,maskmolr,maskmolt);
					//					else
					//						 molr->minWRmsd(molt,matrix4,prof_wrmsd);
					//					//						molt->writePDB("molt1.pdb");
					//					matrix4_op = new M4Rot(matrix4);
					//
					//					score = molt->rmsd(molr,maskmolt,maskmolr); // first score (RMSD)
					//					fprintf(stderr,"   scoret0= %f %f\n",score);
					//
					//
					//
					//					if(model==0) {  // CA-Model
					//						molt2->applyAtoms(matrix4_op);
					//						fprintf(stderr,"  rotate\n",);
					//
					//					}
					//					else
					//						molt->applyAtoms(matrix4_op);
					//					//						molt->writePDB("molt2.pdb");
					//
					//					if(gauss_switch)
					//					{
					//
					//						pdb2gauss_update(molt,gmolt);
					//						//							score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_atomsR,num_atoms,gauss_factor);
					//						score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_gauss,num_atoms,gauss_factor);
					//
					//					}
					//					else
					//					{
					//						//							molt->writePDB("molt.pdb");
					//						//							molr->writePDB("molr.pdb");
					//						if(aln_switch)
					//							score = molt->rmsd(molr,maskmolt,maskmolr); // first score (RMSD)
					//						else
					//							score = molt->rmsd(molr); // score computed with selected model
					//					}
					//					//						exit(0);
				}
				else
				{
					if(model==0) // If CA model
					{
						if(aln_switch)
							molr2->minRmsd(molt2,matrix4,maskmolr2,maskmolt2);
						else
							molr2->minRmsd(molt2,matrix4);
						matrix4_op = new M4Rot(matrix4);
						molt2->applyAtoms(matrix4_op);
						if(gauss_switch)
						{
#ifdef GAUSS

							pdb2gauss_update(molt,gmolt);
							//								score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_atomsR,num_atoms,gauss_factor);
							score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_gauss,num_atoms,gauss_factor);
#endif
						}
						else
							if(aln_switch)
								score = molt->rmsd(molr,maskmolt,maskmolr); // first score (RMSD)
							else
								score = molt->rmsd(molr); // score computed with selected model
					}
					else // If C5 or HA models
					{
						if(!gauss_switch)
							if(aln_switch)
								score = molr->minRmsd(molt,matrix4,maskmolr,maskmolt);
							else
								score = molr->minRmsd(molt,matrix4);
						matrix4_op = new M4Rot(matrix4); // <-- Why this?
						molt->applyAtoms(matrix4_op); // <-- Why this?
						if(gauss_switch)
						{
#ifdef GAUSS

							pdb2gauss_update(molt,gmolt);
							score =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gmolt,num_gauss,num_atoms,gauss_factor);
#endif
						}
					}
				}


				// Pablo
				// free(matrix4_op);
				delete(matrix4_op);
			}
			else // POSE REFINEMENT (PDB-VOL)
			{
				// Initial correlation
				corr = 1 - score;

				// Chain 6D pose refinement ordering choice
				if( sel_chain >= 0 )
				{
					//						chain_order[0] = (int) ( (float) rg->Random() * n_chain ); // [0:1) Playing dice with Mersenne!
					chain_order[0] = (int) ( rg->Random() * n_chain ); // [0:1) Playing dice with Mersenne!
					for(int i=1;i<n_chain;i++)
					{
						if( chain_order[i] < 0 || chain_exist ) // if not-selected yet
						{
							chain_exist = false;
							//								chain_order[i] = (int) ( (float) rg->Random() * n_chain ); // [0:1) Playing dice with Mersenne!
							chain_order[i] = (int) ( rg->Random() * n_chain ); // [0:1) Playing dice with Mersenne!
							for(int j=0; j<i; j++) // checking previous existence
								if( chain_order[i] == chain_order[j] )
									chain_exist = true; // already selected
							if(chain_exist)
								i--;
						}
					}
				}
				else
					chain_order[0] = -1; // Full molecule 6D pose refinement

				// Showing some info
				if( verbose > 1 )
				{
					printf("%s> chain_order[%d]= ",prog_name,n_chain);
					for(int index_chain=0; index_chain<n_chain; index_chain++)
						printf("%d ",chain_order[index_chain]);
					printf("\n");
				}

				// REVISAR ESTO CUANDOSEA...
				// Lo de "select_cg_CA" se puede evitar!!! Pasar: Mol & MolCA(reference)
				//
				// Pose refinement of a Macromolecule into a Map using: INVERSE PARABOLIC INTERPOLATION
				// Based on: Numerical Recipes, pag 402. 10.2 Parabolic Interpolation (1D x 6)
				// Returns: correlation --> if pose improvement is significative
				//          -1.0 --> if pose improvement is neglected (due to "refine_tol")
				if(model==0) // CA-model
				{
					corr = pose_ref(
							&molt2, // p_mol --> Scoring and moving Macromolecule (the same for both)
							&volF, // Map successfully moved and updated (if return > 0)
							&volR, // Reference map (untouched)
							&vol_dummy, // Map buffer
							corr, // Initial correlation (fx - point, it saves some time)
							chain_order, // array with non-repeated chain indices for multi-chain pose refinement
							n_chain, // number of chains (saves some time...)
							kernel, dim_vox, filter, model_thr, pad, resolution,
							corr_nozero, avgR, sigR, // filtration related
							shake_mov, shake_rot, // step-size for translation and rotation
							refine_tol, // significative-ness threshold
							select_cg_CA, // Macromolecule *(*fptr)(Macromolecule *)=NULL,
							&molt // p_mol2 --> Initial Scoring Macromolecule
							// (Note that "mol2" will be a reference to an atoms sub-set in "mol")
#ifdef USE_PTHREAD // Enables PThread parallel routines
							,nthreads // Number of threads used in parallelization
#endif
					);
				}
				else
				{
					corr = pose_ref(
							&molt, // p_mol --> Scoring and moving Macromolecule (the same for both)
							&volF, // Map successfully moved and updated (if return > 0)
							&volR, // Reference map (untouched)
							&vol_dummy, // Map buffer
							corr, // Initial correlation (fx - point, it saves some time)
							chain_order, // array with non-repeated chain indices for multi-chain pose refinement
							n_chain, // number of chains (saves some time...)
							kernel, dim_vox, filter, model_thr, pad, resolution,
							corr_nozero, avgR, sigR, // filtration related
							shake_mov, shake_rot, // step-size for translation and rotation
							refine_tol, // significative-ness threshold
							NULL, NULL
#ifdef USE_PTHREAD // Enables PThread parallel routines
							,nthreads // Number of threads used in parallelization
#endif
					);
				}

				if( corr > 0 ) // Pose refinement improved correlation
				{
					vol_temp = volT;
					volT = volF; // Volume sucessfully moved and updated
					volF = vol_temp;
					score = 1 - corr; // update score
				}
				// else --> nothing to do
			}
		}
		else
			if( verbose > 1 )
				printf("%s> Rot.Tras-lational Refinement rejected!\n",prog_name);

		if(time_switch)
		{
			ht_stage.stopTimer();
			t_ali += (float) ht_stage.getElapsedTime();
		}

		// ************************
		// * NORMAL MODE ANALYSIS *
		// ************************
		last_diag = molt->rmsd(molt_rediag); // Current RMSD from last diagonalized macromolecule (molt_rediag)

		if( verbose > 1 )
			//			printf("%s> Diag. normalized score: %f > %f (rediag)? \n",prog_name,(last_diag-score)/first_score,rediag);
			printf("%s> RMSD increment form last diagonalization: %f > %f (rediag)? \n",prog_name,last_diag,rediag);
		//		if( (last_diag-score)/first_score > rediag || rediag_next || n_iter == 1 ) // minimum score increment to trigger diagonalization
		if( last_diag > rediag || rediag_next || n_iter == 1 ) // minimum score increment to trigger diagonalization
		{
			// NMA timer
			ht_stage.startTimer();

			molt_rediag->copy_coordinates(molt); // Update last diagonalized macromolecule

			// If random IC fix
			if( (fixmodel==4 || fixmodel==5) && n_iter != 1)
			{
				// Under full-atom output demand, the initial full-atom model should be updated
				if(fullatom_switch)
				{
					// HA movie...
					molx = new Macromolecule(mol0); // molecule copy from the HA readed one
					molx_model = new Macromolecule(mol0_model); // molecule copy from the readed (C3 if CA or C5 if C5) (it comes from "molt" ol "molt2" if CA)
					// Moving macromolecule from the initial conformation
					move_dihedralMFAx(molx, UU, props0, 1.0, type, 2, fixed0);

					switch(model)
					{
					case 0:
					case 3:
						move_dihedralMCAx(molx_model, UU, props, 1.0, type, model, fixed);
						break;
					default:
						move_dihedralMFAx(molx_model, UU, props, 1.0, type, model, fixed);
						break;
					}

					// Computing transformation into current model
					if(model==0)
						if(aln_switch)
							molt2->minRmsd(molx_model,matrix4); // molt2 is C3 (if CA-model)
						else
							molt2->minRmsd(molx_model,matrix4); // molt2 is C3 (if CA-model)
					else
						molt->minRmsd(molx_model,matrix4); // molt is C5
					matrix4_op = new M4Rot(matrix4);
					// Applying transformation to current C3 or C5 model
					molx_model->applyAtoms(matrix4_op);
					// Applying transformation to the HA model
					molx->applyAtoms(matrix4_op);
					delete(matrix4_op);

					// Update "initial" models
					delete mol0;
					delete mol0_model;
					mol0 = molx;
					mol0_model = molx_model;
				}

				// Re-defining fixation mask
				if(fixmodel == 4)
					size = fixRand(fixed,old_size,fixRand_prob);
				else
					size = fixRandDH(molt,props,fixed,type,model,fixRand_prob);

				if(verbose > 0)
					printf("%s> Randomizing fixed dofs: old_size= %d  fix_prob= %f  --> size= %d\n",prog_name,old_size,fixRand_prob,size);
				if( verbose > 1 )
				{
					printf("%s> Input CG-model Fixed Internal Coordinates: %d\n",prog_name, old_size-size );
					printf("%s> Input CG-model Mobile Internal Coordinates (size) = %d\n",prog_name,size);
				}

				// Free Hessian and Kinetic matrices (if random ICs choice!)
				// MON: These free and malloc are only needed if different "size", consider removal...
				if(squared_matrices)
				{
					free(hess_sq);
					free(mass_sq);
					hess_sq=NULL;
					mass_sq=NULL;
				}
				else
				{
					free(hess_tr);
					free(mass_tr);
					hess_tr=NULL;
					mass_tr=NULL;
				}

				if(uu!=NULL)
					free(uu);
				uu = (double *) malloc( sizeof(double) * size ); // MON: Where is this initialized? --> In: merge_modes()
				check_pointer(uu,"Current displacement vector");

				if(fullatom_switch)
				{
					// MON: This free and malloc are only needed if different "size", consider removal...
					free(UU);
					UU = (double *) malloc( sizeof(double) * size );
					check_pointer(UU,"Cumulative displacement vector (UU)");
					for(int i=0;i<size;i++)
						UU[i] = 0.0; // initialization

					// UPDATE "fixed0" (full-atom mask) from "fixed" (current model mask)
					changefixIC(mol0, props0, props, fixed0, fixed, 2, type, model, type, true);
					// UPDATE "size0"
					size0 = size; // they share every mobile variable!
				}
			}

			// Getting coordinates single row (pseudo-atom model)
			if( verbose > 1 )
				printf ("%s> Getting coordinates single row\n",prog_name);
			free(coord);
			molt->coordMatrix( &coord );
			if(model==0)
			{
				free(coordNCAC);
				molt2->coordMatrix( &coordNCAC );
			}

			// Computing the PDB's Center of Mass (CoM)
			iter_atom = new pdbIter( molt );
			mtot = 0.0;
			mta = 0.0;
			rd[0] = rd[1] = rd[2] = 0.0;
			for ( iter_atom->pos_atom = 0; !iter_atom->gend_atom(); iter_atom->next_atom() )  // screens all-atoms
			{
				mta = ( iter_atom->get_atom() )->getPdbocc(); // Load mass...
				mtot += mta;
				/* Sum(mass*coord) before putting the CoM at 0 */
				rd[0] += mta * coord[iter_atom->pos_atom * 3];
				rd[1] += mta * coord[iter_atom->pos_atom * 3 + 1];
				rd[2] += mta * coord[iter_atom->pos_atom * 3 + 2];
			}
			delete(iter_atom);
			rd[0] /= mtot;
			rd[1] /= mtot;
			rd[2] /= mtot;
			if( verbose > 1 )
				printf( "Msg(nmafit): Mass %8.8f Center %8.8f %8.8f %8.8f --> Shift to 0,0,0\n", mtot, rd[0], rd[1], rd[2] );

			// shift CM of pdb_model to 0,0,0
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

			// PUT THIS OUTSIDE WHEN RAND= CONSTANT SIZE !!!
			// Eigenvalues memory allocation
			if(eigval!=NULL) // already allocated
				free(eigval);
			eigval  = (double *) malloc(size * sizeof(double));
			check_pointer(eigval,"Eigenvalue memory");

			// Eigenvectors memory allocation and initialization
			if(hess_matrix!=NULL)
				free(hess_matrix);
			hess_matrix = (double *) malloc(size*nevs * sizeof(double));
			check_pointer(hess_matrix,"Eigenvectors memory");

			// NMA begins...
			if(verbose > 0)
				printf("\n");

			// Creating Interacting Pair of Atoms list (IPAs)
			if( verbose > 1 )
				printf("%s> Creating Interacting Pair of Atoms list (IPAs)\n",prog_name);
			ipas = ( twid * ) malloc( 1 * sizeof( twid ) );
			check_pointer(ipas,"Interacting Pair of Atoms list (IPAs) memory");

			// ******************************************
			// * CONTACTING
			// ******************************************
			ht_timer.startTimer();
			t_timer.restart(); // timer
			if( verbose > 0 )
			{
				printf("%s> Building elastic network: ",prog_name);
				fflush(stdout);
			}
			switch(contacts)
			{
			case 0: // INVERSE EXPONENTIAL (power of distance for contact matrix)
				// Making Interacting Pair of (non-virtual) Atoms (ipas)
				make_ipas_new0(molt, &ipas, &nipa, (float) cutoff_k0 );
				if( verbose > 1 )
					printf("\n%s> Inverse Exponential method (%d nipas) cutoff= %f\n",prog_name,nipa,cutoff_k0);
				for(int i=0; i<nipa; i++)
					ipas[i].C = inv_exp(cte_k0,ipas[i].d,x0,power); // setting Force Constants
				break;

			case 1: // DISTANCE CUTOFF METHOD
				// Making Interacting Pair of (non-virtual) Atoms (ipas)
				make_ipas_new0( molt, &ipas, &nipa, (float) cutoff_k1 );
				if( verbose > 1 )
					printf("\n%s> Cutoff Distance method (%d nipas) cutoff= %f\n",prog_name,nipa,cutoff_k1);
				for(int i=0; i<nipa; i++)
					ipas[i].C = cte_k1; // setting Force Constants
				break;

			case 2:	// HINSEN DISTANCE METHOD

				if( verbose > 1 )
					printf("\n%s> Hinsen's distance criterion\n",prog_name);
				for(int i=0; i<nipa; i++)
				{
					if(ipas[i].d <= 4.0)
						ipas[i].C = 86000*ipas[i].d-23900; // setting Force Constants
					else
						ipas[i].C = 128/ (pow( ipas[i].d/10, 6 ) );
				}
				break;

			case 3: // Setting Force Constants according to Secondary Structure and Topology
				if( verbose > 1 )
					printf("\n%s> Topology and/or Secondary Structure method (make_ipasTS)\n",prog_name);
				make_ipasTS(molt, &ipas, &nipa, cutoff_k2, funcs, nfunc, ss_table);
				//				// Write Force constants file (Kfile)
				//				sprintf(text,"%s_kfile.txt",name);
				//				write_Kfile(ipas,nipa,text);
				break;

			case 4: // Laura's "Mixed" model
				// Making Interacting Pair of (non-virtual) Atoms (ipas)
				make_ipas_mix0(molt, &ipas, &nipa);
				if( verbose > 1 )
					printf("imode> \"Mixed\" (%d nipas)\n",nipa);
				break;

			default:
				printf("%s> Please, introduce a valid Connection method to continue!!!\n\nForcing exit!\n\n",prog_name);
				exit(1);
				break;
			}

			// Still working on this... Please check and review!!!
			if(intermolec_switch)
			{
				printf("\n%s> Multiply each inter-molecule force constant by factor= %f\n",prog_name,intermolec_factor);
				// Multiplies by "factor" every ipa force constant ("C") if both atoms belong to different molecules.
				modify_intermolec_ipas(molt, ipas, nipa, intermolec_factor);
			}

			if( verbose > 0 )
			{
				ht_timer.stopTimer();
				if( verbose > 1)
					printf("%s> Elastic Network building time: %lf s\n",prog_name,ht_timer.getElapsedTime() );
				else
					printf(" %lf s\n",ht_timer.getElapsedTime());
				fflush(stdout);
			}

			if(time_switch)
			{
				ht_timer.stopTimer();
				t_springs += (float) ht_timer.getElapsedTime();
			}

			ht_timer.startTimer();
			t_timer.restart(); // timer

			if( verbose > 0 )
			{
				printf("%s> Fast Hessian (rank= %d) ",prog_name,size);
				fflush(stdout);
			}

			// ******************************************
			// * HESSIAN BUILDING
			// ******************************************
			switch(model)
			{
			case 0:
				if( verbose > 1)
					printf("%s> Fast CA-only Hessian Matrix Building O(n^2) [hessianMCAx()] size= %d\n",prog_name,size);
				hessianMCAx(molt,ipas,nipa,size,coordNCAC,coord,&hess_tr,props,unat,fixed);
				break;
			case 1:
			case 2:
				if( verbose > 1)
					printf("%s> Fast Hessian Matrix Building O(n^2) [hessianMFAx()] size= %d\n",prog_name,size);
				hessianMFAx(molt,ipas,nipa,size,coord,&hess_tr,props,unat,type,model,fixed,addrot);
				break;
			case 3:
				if( verbose > 1)
					printf("%s> N,CA,C-model K-matrix with hessianMCA3x() NOT IMPLEMENTED YET\n"
							"Try \"naive methods: K-matrix or V/W-arrays\" instead.\n",prog_name);
				exit(1);
				break;
			}
			free(ipas); // avoids memory leakage!!!

			// Adds torsional springs to a previously built hessian matrix (triangular packing)
			// (ec.5) from:   Lu, Poon and Ma. J.Chem. Theory Comput. 2006, 2, 464-471.
			if(!notors_switch)
				hessianMDHx(molt,props,hess_tr,1,size,fixed,addrot);

			if( verbose > 0 )
			{
				ht_timer.stopTimer();
				if( verbose > 1)
					printf("%s> Hessian matrix building time: %lf s\n",prog_name,ht_timer.getElapsedTime() );
				else
					printf(" %lf s\n",ht_timer.getElapsedTime());
				fflush(stdout);
			}

			if(time_switch)
			{
				ht_timer.stopTimer();
				t_hess += (float) ht_timer.getElapsedTime();
			}

			if(squared_matrices) // If symmetric squared matrices instead of packed triangular storage...
			{
				// Allocating single precision memory
				if( !(hess_sq = (double *) malloc( sizeof( double ) * size*size ) ) )
				{
					printf("Msg(diag): I'm sorry, unable to allocate %d bytes\nForcing exit\n",sizeof( double )*size*size);
					exit(1);
				}
				for ( int i = 0; i < size; i++ )
					for ( int j = i; j < size; j++ ) // diagonal included
					{
						hess_sq[i+size*j] = hess_tr[i + j*(j+1)/2]; // upper triangle + diagonal
						if(i != j) // avoids writing diagonal twice...
							hess_sq[j+size*i] = hess_tr[i + j*(j+1)/2]; // lower triangle
					}
				free(hess_tr); // not needed any more
				hess_tr = NULL; // this way it will be allocated later...
			}

			// ******************************************
			// * KINETIC ENERGY MATRIX BUILDING
			// ******************************************
			ht_timer.startTimer();
			t_timer.restart(); // timer

			if( verbose > 0 )
			{
				printf("%s> Fast Kinetic energy (rank= %d) ",prog_name,size);
				fflush(stdout);
			}
			switch(model)
			{
			case 0:
				if( verbose > 1)
					printf("%s> Fast CA-model Kinetic Matrix Building O(n^2) [kineticMCAx()] size= %d\n",prog_name,size);
				// kineticMCAx( molt, coordNCAC, props, size, &mass_tr, model, fixed ); // H-matrix (triangular)
				kineticMCAxHD( molt, coordNCAC, props, size, &mass_tr, model,fixed); // H-matrix (triangular)

				break;
			case 1:
			case 2:
				if( verbose > 1)
					printf ("%s> Fast Kinetic-Energy matrix Building O(n^2) [kineticMFAx()] size= %d\n",prog_name,size);
				kineticMFAx( molt, coord, props, size, &mass_tr, type, model, fixed, addrot ); // H-matrix (triangular)
				break;
			case 3:
				if( verbose > 1)
					printf("%s> N,CA,C-model K-matrix with kineticMCA3x() NOT IMPLEMENTED YET\n"
							"Try \"naive methods: K-matrix or V/W-arrays\" instead.\n",prog_name);
				exit(1);
				break;
			}

			if( verbose > 0 )
			{
				ht_timer.stopTimer();
				if( verbose > 1)
					printf("%s> Kinetic Energy matrix building time: %lf s\n",prog_name,ht_timer.getElapsedTime() );
				else
					printf(" %lf s\n",ht_timer.getElapsedTime());
				fflush(stdout);
			}

			if(time_switch)
			{
				ht_timer.stopTimer();
				t_kine += (float) ht_timer.getElapsedTime();
			}

			if(squared_matrices) // If symmetric matrices instead of packed storage...
			{
				// Allocating single precision memory
				if( !(mass_sq = (double *) malloc( sizeof( double ) * size*size ) ) )
				{
					printf("Msg(diag): I'm sorry, unable to allocate %d bytes\nForcing exit\n",sizeof( double )*size*size);
					exit(1);
				}
				for ( int i = 0; i < size; i++ )
					for ( int j = i; j < size; j++ ) // diagonal included
					{
						mass_sq[i+size*j] = mass_tr[i + j*(j+1)/2]; // upper triangle + diagonal
						if(i != j) // avoids writing diagonal twice...
							mass_sq[j+size*i] = mass_tr[i + j*(j+1)/2]; // lower triangle
					}
				free(mass_tr); // not needed any more
				mass_tr = NULL; // this way it will be allocated later...
			}

			ht_timer.startTimer();
			t_timer.restart();

			//			// this should be outside!!!
			//			if(nevs_ratio <= diag_choice_cutoff)
			//			{
			//				if( verbose > 0)
			//					printf("%s> Diagonalization DSPGVX  nevs= %d size= %d  ",prog_name,nevs,size);
			//				fflush(stdout);
			//				diag_info = diag_dspgvx(hess_tr, mass_tr, eigval, hess_matrix, size, nevs); // selected eigenvectors
			//			}
			//			else
			//			{
			//				if( verbose > 0)
			//					printf("%s> Diagonalization DSPGVD  ",prog_name);
			//				fflush(stdout);
			//				diag_dspgvd(hess_tr,mass_tr,eigval,hess_matrix,size); // triangular packing storage (eigenvectors -> hess_matrix)
			//			}

			switch(eigensolver)
			{
			case 0:
				if( verbose > 0)
				{
					printf("%s> Diagonalization with LAPACK/BLAS-based XSPGVX()... ",prog_name);
					fflush(stdout);
				}
				diag_xspgvx(hess_tr, mass_tr, eigval, hess_matrix, size, nevs); // detects with sizeof() the floating point precision
				break;
			case 1:
				if( verbose > 0)
				{
					printf("%s> Diagonalization with ARPACK-based dsdrv1_AP_BP_W_mon()... ",prog_name);
					fflush(stdout);
				}
				dsdrv1_AP_BP_W_mon(hess_tr, mass_tr, eigval, hess_matrix, size, nevs);
				break;
			case 2:
				if( verbose > 0)
				{
					printf("%s> Diagonalization with ARPACK-based dsdrv1_A_B_W_mon()... ",prog_name);
					fflush(stdout);
				}
				dsdrv1_A_B_W_mon(hess_sq, mass_sq, eigval, hess_matrix, size, nevs);
				break;
			}

			if(squared_matrices) // If symmetric matrices instead of packed storage...
			{
				free(hess_sq);
				free(mass_sq);
			}

			// Diagonalization "2nd chances"
			// MON: If the "bug" was removed, the following lines may be deleted... check it!
			if(diag_info != 0) // If diagonalization fails...
			{
				if(rediag_trials == 0 || rediag_trials == backup_iter)
				{
					if(dump_switch) // Enables verbose...
					{
						sprintf(dummy_string,"%s_i%d_molt_%d.pdb",name,n_iter,rediag_trials);
						molt->writePDB( dummy_string, false ); // renumbers the PDB
						fprintf(stderr,"\nSorry, diag_info= %d. Critic PDB saved: %s\n",diag_info,dummy_string);
						sprintf(dummy_string,"%s_i%d_molt2_%d.pdb",name,n_iter,rediag_trials);
						molt2->writePDB( dummy_string, false ); // renumbers the PDB
						fprintf(stderr,"Sorry, diag_info= %d. Critic PDB saved: %s\n",diag_info,dummy_string);
					}
				}

				if(rediag_trials >= backup_iter) // If more than "backup_iter" times or so... then exit...
				{
					fprintf(stderr,"Sorry, I've tried %d times to diagonalize... (critic conformation).\nCowardly exiting!\n",rediag_trials);
					exit(3);
				}

				// Restoring previous valid eigenvectors/values (current were destroyed by diag-routine)
				rediag_next = true;
				for(int i=0; i<size; i++)
					eigval[i] = eigval2[i];
				for(int i=0; i<size*nevs; i++)
					hess_matrix[i] = hess_matrix2[i];

				rediag_trials++;
			}
			else // If diagonalization Succeeded!
			{
				rediag_next = false;
				rediag_trials = 0;
				// Backup eigenvectors/values
				for(int i=0; i<size; i++)
					eigval2[i] = eigval[i];
				for(int i=0; i<size*nevs; i++)
					hess_matrix2[i] = hess_matrix[i];
			}

			// Showing output...
			if( verbose > 1 )
			{
				printf("%s> Showing the first 10 eigenvalues:\n",prog_name);
				printf("%s>\n%s> %5s %12s\n","MODE","EIGENVALUE",prog_name,prog_name);
				for(int i=0; i<10; i++)
					printf("%s> %5d %12.5e\n",prog_name,i+1, eigval[i]);
			}

			// MON: Consider removing these lines, I think they are never used...
			// Normalization is not "Normal"-Mode Analysis like !!!! (14/9/2009)
			if(norm_evec_switch)
				norm_evec(hess_matrix,nevs,size); // Standard (norm=1) normalization (strictly needed!)
			// We should optimize "step" without this! But it seems to be neeededdd!!!
			// NOTE: NM-amplitude should be set according to its energy! Check further, please!

			// MON SERVER
			if(server)
			{
				if(more_rmsds)
					printf("%s> %5d %9.6f %9.6f",prog_name,n_iter,score,score_CA);
				else
					printf("%s> %5d %9.6f",prog_name,n_iter,score);
				fflush(stdout);
			}

			ht_timer.stopTimer();
			ht_stage.stopTimer();
			nmatime = ht_stage.getElapsedTime();
			printf(" %3d %6.2f s\n",n_diag,nmatime);
			fflush(stdout);

			if(time_switch)
			{
				t_diag += (float) ht_timer.getElapsedTime();
				t_nma += (float) ht_stage.getElapsedTime();
			}
			ht_stage.startTimer();

			// Set NM excitation probabilities: "VARiance" METHOD and "gaussian" and "invexp"
			if( prob_method == 2 ) // if "var"
			{
				if( verbose > 2 )
					printf("%s> Normal Mode excitation probabilities: \"var\" METHOD\n",prog_name);
				// Normalizing selected eigenvalues to 1 probability
				for( int i=0; i<nevs; i++)
					relev[i] = eigval[0]/eigval[i]; // MSD (variance, normalized respect first eigenvalue)
				// Modesto asked about this! --> MSD or RMSD ???
				//					relev[i] = 1/sqrt(eigval[i]); // RMSD

				if( verbose > 2 )
					for( int i=0; i<nevs; i++)
						printf("%s> %2d %8.6f\n",prog_name,i+1,relev[i]);
			}

			// Normalizes the probability profile (Allows to specify a minimum value)
			if( prob_method != 1 ) // Not needed when "prob = --plain"
				norm_profile(relev,nevs,min_prob,1);

			if(out_files)
			{
				FILE *f_relev;
				sprintf(text,"%s_prob_%d.txt",name,n_diag); // ptraj name
				f_relev=fopen(text,"w");
				for( int i=0; i<nevs; i++ )
					// printf("%5d %12.6f\n",i+1,relev[i]);
					fprintf(f_relev,"%5d %12.6f\n",i+1,relev[i]);
				if( verbose > 2 )
					printf("%s> Relevance file %s written!\n",prog_name,text);
				fclose(f_relev);
			}

			// WEIGHTED-RMSD computation
			if(wrmsd_switch)
				if(pdb_switch) // PDB-PDB WRMSD
				{
					if( verbose > 1 )
					{
						printf("%s>  Computing the \"weighting RMSD\" profile for PDB-PDB--> ",prog_name);
						fflush(stdout);
					}
					if(aln_switch)
						molt->gaussian_weight(molr,&prof_wrmsd,maskmolt,maskmolr,gauss_c);
					else
						molt->gaussian_weight(molr,&prof_wrmsd,gauss_c);
					//					printf("\nJUST CALCULATED: ");
					//					for(int kk=0;kk<num_atoms;kk++)
					//						printf("%f ",prof_wrmsd[kk]);
					//					printf("\n");

					//					for(int x=0;x<molt->get_num_atoms();x++)
					//						printf("%f ",prof_wrmsd[x]);
					//					printf("\n");

					norm_profile(prof_wrmsd, molt->get_num_atoms(), min_wrmsd, 1.0);
				}
				else // PDB-VOL WRMSD
				{
					if( verbose > 1 )
					{
						printf("%s>  Computing the \"weighting RMSD\" profile for PDB-VOL--> ",prog_name);
						fflush(stdout);
					}

					project_map2pdb(volR, molt,&prof_corr,5,dim_vox); // It's more less the same to do this with "Target" instead of "corrmap"
					prof_corr_soft = smooth_profile(prof_corr, num_res, 2, 0);
					norm_profile(prof_corr_soft, num_res, min_wrmsd, 1.0);

					// Atomic level profile for weighted minWRMSD
					array_res2atom(molt,prof_corr_soft,&prof_wrmsd); // for Weighted minRMSD

					if(out_files) // CHECK THIS...
					{
						dummy = new Macromolecule(molt);
						dummy->exchange_Pdbfact(prof_corr_soft);
						sprintf(text,"%s_corr_%d.pdb",name,n_diag);
						dummy->writePDB(text);
						delete dummy;
					}
					free(prof_corr_soft);

					if( verbose > 1 )
					{
						printf("DONE!\n");
						fflush(stdout);
					}
				}
			n_diag++;

			if(time_switch)
			{
				ht_stage.stopTimer();
				t_after += (float) ht_stage.getElapsedTime();
			}
		}

		// LETs MOVE
		if(time_switch)
		{
			ht_stage.startTimer();
			ht_substage.startTimer(); // Select timer
		}

		// Copy trial macromolecule coordinates (faster than allocation...)
		if(model == 0)
		{   // CA-model
			// should have N and C atoms (just to define dihedral directions)
			forward2->copy_coordinates(molt2); // copy current coords.
			// "forward" references "forward2", so nothing to be done here!
		}
		else // 3BB2R- or Full-atom- models
			forward->copy_coordinates(molt); // copy current coords. (faster)

		if(time_switch)
		{
			ht_substage.stopTimer();
			t_trial += (float) ht_substage.getElapsedTime();
			ht_substage.startTimer(); // Merging timer
		}

		// SELECTING MODES
		// Choosing step
		if(rand_step_switch)
			if(step_end >= step_begin) // when random step choice, it's the same...
				step = (double) (step_end-step_begin) * rg->Random() + step_begin; // [step_begin,step_end) (Mersenne)
			else
				step = (double) (step_begin-step_end) * rg->Random() + step_end; // [step_begin,step_end) (Mersenne)
		else
			step = ( ( step_end - step_begin )/max_iter ) * (n_iter+1) + step_begin; // Linearly decreasing step size (arbitrary units)

		if( verbose > 1 )
			printf("%s> Current Step = %f\n",prog_name,step);

		// Choosing number of excited modes
		if(rand_excited_switch)
			if(nex2>=nex)
				excite = (int) ( (nex2-nex+1) * (float)rg->Random() ) + nex; // [nex,nex2]
			else
				excite = (int) ( (nex-nex2+1) * (float)rg->Random() ) + nex2; // [nex,nex2]
		else
			excite = ( (float)( nex2 - nex )/(float)max_iter ) * (float) (n_iter+1) + nex; // Linearly decreasing number of excited modes

		if( verbose > 1 )
			printf("%s> Current excited modes = %d (nex= %d  nex2= %d)\n",prog_name,excite,nex,nex2);

		// It selects and weights "num" modes according to a given probability profile.
		if(scv_weight)
			select_mode_prob(relev, eigval, nevs, excite, &selected, prob_cutoff, randweight); // Uses "eigval" to apply SCV scheme in addition to random weights
		else
			select_mode_prob(relev, NULL, nevs, excite, &selected, prob_cutoff, randweight); // Just using random weights

		if(time_switch)
		{
			ht_substage.stopTimer();
			t_select += (float) ht_substage.getElapsedTime();
			ht_substage.startTimer(); // Merging timer
		}

		// MERGING MODES
		// It Merges "nev" modes, selected in "nm_props"
		// (warning: it does not allocate memory! (uu))
		merge_modes(uu, hess_matrix, size, selected, excite);
		if(time_switch)
		{
			ht_substage.stopTimer();
			t_merge += (float) ht_substage.getElapsedTime();
		}

		if(time_switch)
			ht_substage.startTimer(); // Moving timer

		// MOVING MACROMOLECULE (with dihedrals)
		switch(model)
		{
		case 0:
			move_dihedralMCAx(forward2, uu, props, step, type, model, fixed); // "forward2" is C3 model
			break;
		case 3:
			move_dihedralMCAx(forward, uu, props, step, type, model, fixed);
			break;
		default:
			move_dihedralMFAx(forward,uu,props,step,type,model,fixed,addrot);
			break;
		}

		if(time_switch)
		{
			ht_substage.stopTimer();
			t_mov_nm += (float) ht_substage.getElapsedTime();
		}

		// ALIGNMENT (align_method)
		if(time_switch)
			ht_substage.startTimer(); // Alignment timer

		// Alignment respect to the current PDB model !!! IT WORKS!!!
		if( verbose > 1 )
			printf("%s> Aligning Next model to the current one!\n",prog_name);

		//		forward->writePDB("forward0.pdb");

		if(wrmsd_switch)
			molt->minWRmsd(forward,matrix4,prof_wrmsd); // Not-cheating
		else
			molt->minRmsd(forward,matrix4); // Not-cheating
		matrix4_op = new M4Rot(matrix4);
		if(model==0)
			forward2->applyAtoms(matrix4_op); // "forward" will be also moved because it references to "forward2" atoms
		else
			forward->applyAtoms(matrix4_op);
		delete(matrix4_op);

		//		forward->writePDB("forward1.pdb");

		if(time_switch)
		{
			ht_substage.stopTimer();
			t_mov_ali += (float) ht_substage.getElapsedTime();
			ht_substage.startTimer();
			ht_stage.stopTimer();
			t_mov += (float) ht_stage.getElapsedTime();
			ht_stage.startTimer();
		}

		// CHECKING MAXIMUM DISPLACEMENT - BEGIN
		//		current_dmax = checkmax(molt,forward);
		//		avg_dmax += current_dmax;
		//		if(current_dmax > max_dmax)
		//		{
		//			max_dmax = current_dmax;
		//			iter_dmax = n_iter;
		//		}
		// CHECKING MAXIMUM DISPLACEMENT - END

		// COMPUTING NEW SCORES
		if(pdb_switch) // PDB-PDB Morphing
		{
			if(gauss_switch)
			{
#ifdef GAUSS
				pdb2gauss_update(forward,gforward);
				//				score_f =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gforward,num_atomsR,num_atoms,gauss_factor);
				score_f =	Overlap_Integral_Bwn_GAUSS_MOLECULEs(gmolr,gforward,num_gauss,num_atoms,gauss_factor);
#endif
			}
			else
				if(aln_switch)
					score_f = forward->rmsd(molr,maskmolt,maskmolr);
				else
					score_f = forward->rmsd(molr);
			//			forward->writePDB("forward.pdb");
			//			molr->writePDB("molr.pdb");
			//			molt->writePDB("molt.pdb");
			//			fprintf(stderr,"score_f= %f\n",score_f);
			//			exit(0);
		}
		else // PDB-VOL
		{
			if(time_switch)
				ht_substage.startTimer(); // Projection timer

			if( verbose > 1 )
				printf( "%s> Projecting Model into Volume (filter=%d)\n",prog_name,filter);

			// PROJECTING MOVED MACROMOLECULE INTO A MAP
#ifdef USE_PTHREAD // Enables PThread parallel routines
			project_pdb(forward,&volF,&vol_dummy,kernel,dim_vox,filter,resolution,fast_switch,nthreads,threads_data);
#else
			project_pdb(forward,&volF,&vol_dummy,kernel,dim_vox,filter,resolution,fast_switch);
#endif


			//if(n_iter==100)
			//	FOPS::writeFile(volT,"volT_3.sit");
			//if(n_iter==201)
			//{
			//	FOPS::writeFile(volT,"volT_4.sit");
			//	exit(0);
			//}
			// If there are some negative values, cross-corr may be negative... to avoid this:
			if(filter==1)
				FOPS::threshold( volF, model_thr );

			if(time_switch)
			{
				ht_substage.stopTimer();
				t_mov_proj += (float) ht_substage.getElapsedTime();
				ht_substage.startTimer(); // Projection timer
			}

			if(norm_switch)
				norm_vol(volF,1.0,norm_factor); // Normalization

			score_f = 1 - FOPS::correlation_frame(volF,volR,pad,corr_nozero,false,avgR,sigR);
			if(time_switch)
			{
				ht_substage.stopTimer();
				t_corr += (float) ht_substage.getElapsedTime(); // correlation timer
				ht_substage.startTimer();
			}
		}

		// SHOWING THE BEST SCORES
		if( verbose > 1 )
			printf("%s> SCORES:\tOld= %9.6f\tCurrent= %9.6f\t\n",prog_name,score,score_f);
		if(time_switch)
		{
			ht_stage.stopTimer();
			t_score += (float) ht_stage.getElapsedTime();
			ht_stage.startTimer();
		}

		// Selecting the NEXT SOLUTION
		next_score = score_f;

		if(!pdb_switch)
			next_vol = &volF;

		// Fixes low-resolution BUG (15/1/2009)
		// During easy fittings and with low-resolutions, the scores are so low
		// (arround 1e-6/-7, i.e. the float precision range) that can be
		// lower than zero (score<0), thus producing senseless score values.
		if(next_score<0)
		{
			next_score = 0;
			do_main_loop = false; // Forcing exit
			n_iter = max_iter-1; // just to be compatible with some scripts...
		}

		// Acceptance test
		if( next_score-score < 0 || rediag_next) // better next score || (rediag.needed) --> ACCEPTED!
		{
			if( verbose > 1 )
				printf("%s> Trial motion Acepted (score= %f  next_score= %f)\n",prog_name,score,next_score);

			// Updating absolute conformation in internal coords. (UU)
			if(fullatom_switch)
				for(int i=0; i<size; i++)
					UU[i] += uu[i] * step; // Merged-Mode * step = absolute conf.

			// Swapping models (faster)
			buff=molt;
			molt=forward;
			forward=buff;
			if(model==0)
			{
				buff=molt2;
				molt2=forward2;
				forward2=buff;
			}
			if(gauss_switch)
			{
#ifdef GAUSS

				gbuff = gmolt;
				gmolt = gforward;
				gforward = gbuff;
#endif
			}

			if(!pdb_switch) // PDB-Map fitting
			{
				vol_temp = volT; // Swapping T-->Next & Next-->T (to avoid memory leaks!)
				volT = *next_vol;
				*next_vol = vol_temp;
			}
			score_old = score;
			score = next_score;
		}
		else // if it's not accepted or T-ambient determination... (next move rejected!)
		{
			if( verbose > 1 )
				printf("%s> Trial motion Rejected (score= %f  next_score= %f)\n",prog_name,score,next_score);
		}

		// RMSDs computation for CA-model
		//		if(more_rmsds && (pdb_switch || bench_switch) && ( verbose > 1 || n_iter % each_verb == 0))
		//		{
		//			delete molt_CA;
		//			if(model==0)
		//				molt_CA = molt2->select_cpy(calpha2); // selecting again because pointers were destroyed
		//			else
		//				molt_CA = molt->select_cpy(calpha2); // selecting again because pointers were destroyed
		//			if(aln_switch)
		//				score_CA = molr_CA->rmsd(molt_CA,maskmolrCA,maskmoltCA); // CA score (RMSD)
		//			else
		//				score_CA = molr_CA->rmsd(molt_CA); // CA score (RMSD)
		////			fprintf(stderr,"score_CA= %f\n",score_CA);
		//		}

		// Mon (2/2/2015): Some bug related to the trigger of movie-frames-saving solved...
		if(movie_switch || more_rmsds || pdb_switch || bench_switch)
		{
			delete molt_CA;
			if(model==0)
				molt_CA = molt2->select_cpy(calpha2); // selecting again because pointers were destroyed
			else
				molt_CA = molt->select_cpy(calpha2); // selecting again because pointers were destroyed
		}
		if(more_rmsds && (pdb_switch || bench_switch) && ( verbose > 1 || n_iter % each_verb == 0))
		{
			if(aln_switch)
				score_CA = molr_CA->rmsd(molt_CA,maskmolrCA,maskmoltCA); // CA score (RMSD)
			else
				score_CA = molr_CA->rmsd(molt_CA); // CA score (RMSD)
			//			fprintf(stderr,"score_CA= %f\n",score_CA);
		}

		// Score Visualization
		if(pdb_switch) // PDB-PDB (morphing)
		{
			if(!server)
			{
				if( verbose > 1 || n_iter % each_verb == 0)
				{
					if(more_rmsds)
						printf("\r%s> %5d %9.6f %9.6f",prog_name,n_iter,score,score_CA);
					else
						printf("\r%s> %5d %9.6f",prog_name,n_iter,score);
					fflush(stdout);
				}
			}
			if(n_iter % each_score == 0)
				if(more_rmsds)
					fprintf(f_score,"%6d %11.4e %11.4e\n",n_iter,score,score_CA); // saving score
				else
					fprintf(f_score,"%6d %11.4e\n",n_iter,score); // saving score
		}
		else // PDB-VOLUME (fitting)
		{
			if( verbose > 1 || n_iter % each_verb == 0)
			{
				if(bench_switch && more_rmsds)
					printf("\r%s> %5d %9.6f %9.6f %9.6f",prog_name,n_iter,score,score_CA,1-score);
				else
					printf("\r%s> %5d %9.6f %9.6f",prog_name,n_iter,score,1-score);
				fflush(stdout);
			}
			// Score Saving (PDB-VOL)
			if(n_iter % each_score == 0)
				if(bench_switch && more_rmsds)
					fprintf(f_score,"%6d %11.4e %11.4e %11.4e %11.4e %11.4e\n",n_iter,score,score_CA,1-score,next_score-score,delta_avgscore);
				else
					fprintf(f_score,"%6d %11.4e %11.4e %11.4e %11.4e\n",n_iter,score,1-score,next_score-score,delta_avgscore);
		}

		// Convergence stuff...
		if(!rediag_next) // only if rediagonalization is not needed
		{

			// Convergence criteria
			scores = ( float * ) realloc( scores, n_iter * sizeof( float ) ); // resizes scores list
			scores[n_iter - 1] = score; // Current score

			if( n_iter <= conv_length + stage )
				score_avg1 += score;
			if( n_iter == conv_length + stage )
				score_avg1 /= conv_length;
			if( n_iter > conv_length + stage && n_iter <= 2*conv_length + stage )
				score_avg2 += score;
			if( n_iter == 2*conv_length + stage)
				score_avg2 /= conv_length;
			if( n_iter > 2*conv_length + stage)
			{
				// Updating moving averages
				score_avg1 = score_avg1 + ( (scores[n_iter - conv_length -1] - scores[n_iter - 2*conv_length -1] ) / conv_length );
				score_avg2 = score_avg2 + ( (scores[n_iter -1] - scores[n_iter - conv_length -1] ) / conv_length );
				delta_avgscore = score_avg1 - score_avg2;

				if( verbose > 1 )
					//					printf("%s> delta_avgscore/first_score= %f: < %f (convergence) ?\n",prog_name,delta_avgscore/first_score,convergence);
					printf("%s> delta_avgscore= %f: < %f (convergence) ?\n",prog_name,delta_avgscore,convergence);

				//				if(delta_avgscore/first_score < convergence)
				if(delta_avgscore < convergence)
				{
					do_main_loop = false; // Forces to exit from the main loop
					//					printf("\n%s> Convergence reached at iter %d:  %.2e < %.2e",prog_name,n_iter,delta_avgscore/first_score,convergence);
					printf("\n%s> Convergence reached at iter %d (%5e<%.2e)",prog_name,n_iter,delta_avgscore,convergence);
				}
			}
		}
		if( verbose > 1 )
			printf("%s>  score_avg1= %f  score_avg2= %f  delta_avgscore= %f  score=%f\n",prog_name,score_avg1,score_avg2,delta_avgscore,scores[n_iter - 1]);

		// SAVING Output when significative improvement is detected
		if(!rediag_next) // only if rediagonalization is not needed
			if(delta_save>0)
			{
				if(deltasave_rmsd_switch)
				{
					// computing current delta-rmsd (rmsd)
					//				delete molt_CA;
					//				molt_CA = molt->select_cpy(calpha2); // selecting again because pointers were destroyed
					//					molt_CA->copy_coordinates(molt);
					rmsd = molt_CA0->rmsd(molt_CA); // CA score (RMSD)
					//					molt_CA->writePDB("molt_CA.pdb");
					//					molt_CA0->writePDB("molt_CA0.pdb");
					//exit(0);
					if( verbose > 1 )
						printf("%s> Saving deltasave_rmsd_switch? rmsd=%f delta=%f delta_save=%f\n",prog_name,rmsd,rmsd,delta_save);
					// Mon (2/2/2015): Some bug related to the trigger of movie-frames-saving solved...
					if( rmsd > delta_save )
						out_switch = true;
				}
				else
				{
					if( verbose > 1 )
						printf("%s> Saving norm. delta_save= %f > %f ?\n",prog_name,(score_lastsaved-score)/first_score,delta_save);
					if( (score_lastsaved-score)/first_score > delta_save ) // saving score after "save_output" iterations
						out_switch = true;
				}
			}
			else // Output files each "-delta_save" iterations
			{
				if( verbose > 1 )
					printf("%s> Saving output each %d (-delta_save) iters\n",prog_name,(int)-delta_save);
				if( n_iter % ((int)-delta_save) == 0 ) // saving score after "-delta_save" iterations
					out_switch = true;
			}

		if(out_switch && !rediag_next) // saving output files // only if rediagonalization is not needed
		{
			if( verbose > 1 )
				printf("%s> Saving Output! (iter= %d)\n",prog_name,n_iter);
			fflush(stdout);
			if(deltasave_rmsd_switch)
			{
				// Updates last-saved rmsd
				//score_lastsaved = rmsd;
				// Updates saved frame CA-model (from which the following RMSD's will be computed)
				//				delete molt_CA0;
				//				molt_CA0 = new Macromolecule(molt_CA);
				molt_CA0->copy_coordinates(molt_CA);
			}
			else
				score_lastsaved = score;  // updates last-saved score

			fclose(f_score); // trick... (saving output files, the scores file is also updated and flushed)
			if(	!( f_score = (FILE *)fopen(file_score,"a+") ) )
			{
				printf("\nAn ERROR ocurred during \"file_score\" append...\nForcing Exit!\n\n");
				exit(1);
			}

			if(movie_switch)
			{
				sprintf(text,"%s_movie.pdb",name);
				if(fullatom_switch)
				{
					// Full-Atom movie...
					molx = new Macromolecule(mol0); // molecule copy from "readed" one
					molx_model = new Macromolecule(mol0_model); // molecule copy from "readed-3BB2R" one
					// Moving macromolecule from the initial conformation
					move_dihedralMFAx(molx, UU, props0, 1.0, type, 2, fixed0);
					switch(model)
					{
					case 0:
					case 3:
						move_dihedralMCAx(molx_model, UU, props, 1.0, type, model, fixed);
						break;
					default:
						// Watch out full-atom "addrot" here...
						move_dihedralMFAx(molx_model, UU, props, 1.0, type, model, fixed);
						break;
					}

					// Internal DoFs are always applied over the same directions and positions
					// but Inter-segment DoFs are applied over the center of mass of current
					// model (when it happened), thus, given we don't store center of mass trajectory,
					// it's necessary to align each segment separately!
					if(n_seg > 1)
					{
						for(int i=0; i<n_seg; i++)
						{
							// Computing transformation into current model ("model")
							// Watch out SMOL's here...
							selxm = select_segment(molx_model,i);
							selx = select_segment(molx,i);
							if(model==0)
								selm = select_segment(molt2,i);
							else
								selm = select_segment(molt,i);
							selm->minRmsd(selxm,matrix4);

							matrix4_op = new M4Rot(matrix4);
							// Applying transformation to current model
							selxm->applyAtoms(matrix4_op);
							// Applying transformation to the Full-Atom model
							selx->applyAtoms(matrix4_op);
							delete(matrix4_op);
							selxm->erase_level(pdb_segment);
							delete selxm;
							selx->erase_level(pdb_segment);
							delete selx;
							selm->erase_level(pdb_segment);
							delete selm;
						}
					}
					else
					{
						if(model==0)
							molt2->minRmsd(molx_model,matrix4); // "molx_model" comes from "molt2" if CA
						else
							molt->minRmsd(molx_model,matrix4); // "molx_model" comes from "molt" if C5 or HA
						matrix4_op = new M4Rot(matrix4);
						// Applying transformation to the Full-Atom model
						molx->applyAtoms(matrix4_op);
						delete(matrix4_op);
					}
					// Saving aligned Full-Atom model into Multi-PDB
					movie_index++; // movie frame index to conform Multi-PDB format an Chimera's stuff
					molx->writeMPDB(text,movie_index);

					if(chimera_switch) // Specific output for UCSF-Chimera plugin
					{
						// Save "_fitted.pdb" every time a movie frame is generated (Full-atom)
						sprintf(text,"%s_fitted.pdb",name); // Selected atom-model
						molx->writePDB(text);
					}

					delete molx;
					delete molx_model;
				}
				else
				{
					movie_index++; // movie frame index to conform Multi-PDB format an Chimera's stuff
					molt->writeMPDB(text,movie_index);

					if(chimera_switch) // Specific output for UCSF-Chimera plugin
					{
						// Save "_fitted.pdb" every time a movie frame is generated (at CG-level of the simulation)
						sprintf(text,"%s_fitted.pdb",name); // Selected atom-model
						molt->writePDB(text);
					}
				}

				// "New frame" tag for Chimera auto-refresh...
				if(chimera_switch) // Specific output for UCSF-Chimera plugin
#if defined  _WIN32 || defined _WIN64 // Windows version
					fprintf(stdout,"new_frame"); // In Windows, python's pipes behavior seems different wrt stdout and stderr...
#else // Linux version
				fprintf(stderr,"new_frame");
#endif

				if( verbose > 1 )
					printf("%s> Added model %5d to %s\n",prog_name,n_iter,text);
			}
			if(pdblog_switch)
			{
				// Save iteration PDB
				sprintf(text,"%s_%d.pdb",name,n_iter);
				molt->writePDB(text);
				if( verbose > 1 )
					printf("%s> Written model %5d (%s)\n",prog_name,n_iter,text);
			}
			out_switch = false;
		}

		if(time_switch)
		{
			ht_stage.stopTimer();
			t_next += (float) ht_stage.getElapsedTime();
			ht_stage.startTimer();
		}

		if(do_main_loop && !rediag_next) // and rediag. is not needed
			n_iter++;
	}
	t_main = time_main.elapsed();
	// ******************************************************************************************
	// * END MAIN LOOP
	// ******************************************************************************************
	delete molt_rediag; // delete last diagonalized macromolecule

	// Last RMSDs computation
	if(pdb_switch || bench_switch) // PDB-PDB Morphing || Benchmark Docking
	{
		delete molt_CA;
		if(model==0)
			molt_CA = molt2->select_cpy(calpha2); // selecting again because pointers were destroyed
		else
			molt_CA = molt->select_cpy(calpha2); // selecting again because pointers were destroyed


		if(aln_switch) {
			// bug Pablo
			sprintf(text,"%s_final.rmsd",name);
			score_CA = molr_CA->rmsd_file(molt_CA,maskmolrCA,maskmoltCA, text); // CA score (RMSD)
			score_CA = molr_CA->rmsd(molt_CA,maskmolrCA,maskmoltCA); // CA score (RMSD)

		}
		else
			score_CA = molr_CA->rmsd(molt_CA); // CA score (RMSD)
	}


	// Adding last line to score-file
	if(pdb_switch) // PDB-PDB
		fprintf(f_score,"%6d %11.4e %11.4e\n",n_iter,score,score_CA); // saving score
	else // PDB-VOLUME
		fprintf(f_score,"%6d %11.4e %11.4e %11.4e %11.4e %11.4e\n",n_iter,score,score_CA,1-score,next_score-score,delta_avgscore);

	// Showing last score line...
	if( verbose > -1 )
		if(pdb_switch) // PDB-PDB
		{
			printf("\n%s>\n%s> %5s %9s %9s\n",prog_name,prog_name,"Iter","RMSD","RMSD_CA");
			printf("%s> %5d %9.6f %9.6f\n",prog_name,n_iter,score,score_CA);
		}
		else // PDB-VOL
		{
			if(bench_switch)
			{
				printf("\n%s>\n%s> %5s %9s %9s %9s\n",prog_name,prog_name,"Iter","score","RMSD_CA","Corr.");
				printf("%s> %5d %9.6f %9.6f %9.6f\n",prog_name,n_iter,score,score_CA,1-score);
			}
			else
			{
				printf("\n%s>\n%s> %5s %9s %9s\n",prog_name,prog_name,"Iter","score","Corr.");
				printf("\r%s> %5d %9.6f %9.6f\n",prog_name,n_iter,score,1-score);
			}
		}
	printf("%s>\n",prog_name);

	// Final Full-atom model generation
	if( fullatom_switch && (savefitted_switch || movie_switch) )
	{
		// Full-Atom
		molx = new Macromolecule(mol0); // molecule copy from "readed" one
		molx_model = new Macromolecule(mol0_model); // molecule copy from "readed-model" one

		// Moving macromolecule from the initial conformation
		move_dihedralMFAx(molx, UU, props0, 1.0, type, 2, fixed0);
		switch(model)
		{
		case 0:
		case 3:
			move_dihedralMCAx(molx_model, UU, props, 1.0, type, model, fixed);
			break;
		default:
			move_dihedralMFAx(molx_model, UU, props, 1.0, type, model, fixed, addrot);
			break;
		}

		// Internal DoFs are always applied over the same directions and positions
		// but Inter-segment DoFs are applied over the center of mass of the current
		// model when it occurred, thus given we don't store center of mass trajectory,
		// it's necessary to align each segment separately!
		if(n_seg > 1)
		{
			for(int i=0; i<n_seg; i++)
			{
				// Computing transformation into current model ("model")
				// Watch out SMOL's here...
				selxm = select_segment(molx_model,i);
				selx = select_segment(molx,i);
				if(model==0)
					selm = select_segment(molt2,i);
				else
					selm = select_segment(molt,i);
				selm->minRmsd(selxm,matrix4);

				matrix4_op = new M4Rot(matrix4);
				// Applying transformation to current model
				selxm->applyAtoms(matrix4_op);
				// Applying transformation to the Full-Atom model
				selx->applyAtoms(matrix4_op);
				delete(matrix4_op);
				selxm->erase_level(pdb_segment);
				delete selxm;
				selx->erase_level(pdb_segment);
				delete selx;
				selm->erase_level(pdb_segment);
				delete selm;
			}
		}
		else
		{
			// Computing transformation into current model
			if(model==0)
				molt2->minRmsd(molx_model,matrix4);
			else
				molt->minRmsd(molx_model,matrix4);
			matrix4_op = new M4Rot(matrix4);
			// Applying transformation to the Full-Atom model
			molx->applyAtoms(matrix4_op);
			delete(matrix4_op);
		}

		if(bench_switch)
		{
			float temp = molx->dist_profile(molrini,&deviation);
			if( verbose > 0 )
			{
				printf("%s> Residue Averaged Deviations applied to B-factors!\n",prog_name);
				printf("%s> dist_profile's RMSD = %f\n",prog_name,temp);
			}
			molx->exchange_Pdbfact(deviation,true); // Now "molx" has the deviations in B-factors
		}

		if(savefitted_switch)
		{
			// Saving fitted Full-Atom final model
			//			sprintf(text,"%s_fitFULL.pdb",name); // Warning: Full-Atom model
			sprintf(text,"%s_fitted.pdb",name); // Warning: Full-Atom model
			molx->writePDB(text);
			if( verbose > -1 )
				printf("%s> Final Full-Atom model: %35s\n",prog_name,text);


			if(aln_switch) {

			//
			// output dihedrals && SS
			//
			float *di2;
			int *ss, dangi;
			Fragment * res;
			int resn, Nres, simple;
			Segment * seg;
			int chino;
			float * chis;
			Chain *ch;
			pdbIter *iter1, *iter2;
			iter1 = new pdbIter( molrini,  false, false, true, false); //  r maskres2
			iter2 = new pdbIter( molx,  false, false, true, false);    //  t maskres1
			// first SS assign
			molx->all_dihedrals( &di2);
			dihedrals2DISICL( di2, &ss, (iter2->num_fragment()) );


			// first SS assign
			FILE *f_out=NULL;
			sprintf(text,"%s_fitted.dih",name);
			if( !(f_out = fopen(text,"w")) )
			{
				fprintf(stderr,"Sorry, unable to write output file! Forcing exit!\n");
				exit(2);
			}

			fprintf(f_out, "# AA     resn    PHI      PSI     OMEGA     DISICL  simple    chi1     chi2     chi3      chi4\n");
			iter2->pos_fragment = 0;	dangi=0; Nres=0;
			for ( iter1->pos_fragment = 0; !iter1->gend_fragment(); iter1->next_fragment() )
			{
				if(maskres2[iter1->pos_fragment]) {
					while(!iter2->gend_fragment() && !maskres1[iter2->pos_fragment]) // this places the index into the corresponding pas
					{
						iter2->next_fragment();
						dangi+=7;
						Nres++;
					}
					res = ( Fragment * ) iter2->get_fragment();
					//					if(maskres1[iter2->pos_fragment]) {
					simple=DISICL_d2simple[ss[Nres]];  // change DISICL simple SS classes
					fprintf( f_out, "%3s %8d  %8.3f %8.3f %8.3f  %3d %3s %3d %3s  %8.3f %8.3f %8.3f %8.3f\n",
							res->getName(), res->getIdNumber(),
							di2[dangi], di2[dangi+1],  di2[dangi+2], // (phi,psi,omega)
							ss[Nres],DISICL_d[ss[Nres]],	simple, DISICL_s[simple],
							di2[dangi+3], di2[dangi+4], di2[dangi+5], di2[dangi+6]  // (4x Chi)
					);
					iter2->next_fragment();
					dangi+=7;
					Nres++;
					//					}
				}

			}
			iter1->~pdbIter();
			iter2->~pdbIter();


			fclose(f_out);

		 }
		}

		if(movie_switch)
		{
			sprintf(text,"%s_movie.pdb",name);
			// Saving aligned Full-Atom model into Multi-PDB
			movie_index++; // movie frame index to conform Multi-PDB format an Chimera's stuff
			molx->writeMPDB(text,movie_index);
			if( verbose > -1 )
				printf("%s> Movie made:            %35s\n",prog_name,text);
		}

		delete molx;
		delete molx_model;
	}

	// Saving last frame into the movie
	if(movie_switch && !fullatom_switch)
	{
		sprintf(text,"%s_movie.pdb",name);
		movie_index++; // movie frame index to conform Multi-PDB format an Chimera's stuff
		molt->writeMPDB(text,movie_index); // Saving aligned 3BB2R model into Multi-PDB
		if( verbose > -1 )
			printf("%s> Movie file:            %35s\n",prog_name,text);
	}

	if(savemodel_switch)
	{
		sprintf(text,"%s_model.pdb",name);
		printf("%s> Initial Model:         %35s\n",prog_name,text);
		if(model == 0) // in CA-model case, a NCAC model is also written
		{
			sprintf(text,"%s_ncac.pdb",name);
			printf("%s> N,CA,C Model:          %35s\n",prog_name,text);
		}
	}

	if(savefitted_switch)
	{
		if(fullatom_switch)
			sprintf(text,"%s_fitCG.pdb",name);
		else
			sprintf(text,"%s_fitted.pdb",name);

		if( save_deviation && (pdb_switch || bench_switch) ) // (PDB or bench) and save deviation
		{
			printf("%s> Residue Averaged Deviations applied to B-factors!\n",prog_name);
			printf("%s> dist_profile's RMSD = %f\n",prog_name,molt->dist_profile(molr,&deviation));
			molt->exchange_Pdbfact(deviation,true); // now "molt" has distance deviations
		}

		// Writting final model
		molt->writePDB(text); // Warining: it has charges and masses!

		if( verbose > -1 )
			printf("%s> Final Model:           %35s\n",prog_name,text);
	}

	// "New frame" tag for Chimera auto-refresh...
	if(chimera_switch) // Specific output for UCSF-Chimera plugin
#if defined  _WIN32 || defined _WIN64 // Windows version
		fprintf(stdout,"new_frame"); // In Windows, python's pipes behavior seems different wrt stdout and stderr...
#else // Linux version
	fprintf(stderr,"new_frame");
#endif

	if(!pdb_switch) // PDB-VOL
	{
		if(savemaps_switch)
		{
			sprintf(text,"%s_fitted.sit",name);
			FOPS::writeFile(volT,text);
			if( verbose > 1 )
				printf("%s> Fitted & Filtered map: %35s\n",prog_name,text);
		}

		fprintf(f_score,"# Std.corr:  FIRST= %12.8f   LAST= %12.8f  TIME= %6.0f s\n",first_corr,FOPS::correlation_frame(volT,volR,pad,corr_nozero,false,avgR,sigR),(float)t_main);
	} // PDB-PDB
	else
		if(aln_switch)
			fprintf(f_score,"# RMSD:  FIRST= %12.8f   LAST= %12.8f  TIME= %6.0f s\n",first_corr,molt->rmsd(molr,maskmolt,maskmolr),(float)t_main);
		else
			fprintf(f_score,"# RMSD:  FIRST= %12.8f   LAST= %12.8f  TIME= %6.0f s\n",first_corr,molt->rmsd(molr),(float)t_main);
	fclose(f_score);
	sprintf(text,"%s_score.txt",name);
	printf("%s> Score file:            %35s\n",prog_name,text);

	t_sum = t_ali+t_nma+t_after+t_mov+t_score+t_next;
	if(time_switch)
		if( verbose > -1 )
		{
			printf("\nFINAL TIMINGS: (iters: %6d)\n",n_iter);
			printf("%s %8s %7s\n",    "Stage        ","time(s)","time(%)");
			printf("%s %8.2f %7.2f\n","Aligning     ",(float)t_ali, 100*(float)t_ali/t_main);
			printf("%s %8.2f %7.2f\n","NMA:         ",(float)t_nma, 100*(float)t_nma/t_main);
			printf("%s %8.2f %7.2f\n","   Springs   ",(float)t_springs, 100*(float)t_springs/t_main);
			printf("%s %8.2f %7.2f\n","   Hessian   ",(float)t_hess, 100*(float)t_hess/t_main);
			printf("%s %8.2f %7.2f\n","   Kinetic   ",(float)t_kine, 100*(float)t_kine/t_main);
			printf("%s %8.2f %7.2f\n","   Diag.     ",(float)t_diag, 100*(float)t_diag/t_main);
			printf("%s %8.2f %7.2f\n","After-Diag.  ",(float)t_after, 100*(float)t_after/t_main);
			printf("%s %8.2f %7.2f\n","Moving:      ",(float)t_mov, 100*(float)t_mov/t_main);
			printf("%s %8.2f %7.2f\n","   Trial_MDL ",(float)t_trial, 100*(float)t_trial/t_main);
			printf("%s %8.2f %7.2f\n","   Selecting ",(float)t_select, 100*(float)t_select/t_main);
			printf("%s %8.2f %7.2f\n","   Merging   ",(float)t_merge, 100*(float)t_merge/t_main);
			printf("%s %8.2f %7.2f\n","   Mov.NM.   ",(float)t_mov_nm, 100*(float)t_mov_nm/t_main);
			printf("%s %8.2f %7.2f\n","   Mov.Align ",(float)t_mov_ali, 100*(float)t_mov_ali/t_main);
			printf("%s %8.2f %7.2f\n","Scoring:     ",(float)t_score, 100*(float)t_score/t_main);
			printf("%s %8.2f %7.2f\n","   Project   ",(float)t_mov_proj, 100*(float)t_mov_proj/t_main);
			printf("%s %8.2f %7.2f\n","   Corr.     ",(float)t_corr, 100*(float)t_corr/t_main);
			printf("%s %8.2f %7.2f\n","Next-stuff   ",(float)t_next, 100*(float)t_next/t_main);
			printf("%s %8.2f %7.2f\n","SUM          ",(float)t_sum, 100*(float)t_sum/t_main);
			printf("%s %8.2f %7.2f (%6.3f h)\n","Main Loop    ",(float)t_main, 100*(float)t_main/t_main, (float)t_main/3600);

			// Opening Timings file
			FILE *f_time;
			if( time_switch )
			{
				sprintf(text,"%s_time.txt",name);
				if( !(f_time=(FILE *)fopen(text,"w") ) )
				{
					printf("Sorry, unable to open output time file: %s\n",text);
					exit(1);
				}
			}

			fprintf(f_time,"\n# FINAL TIMINGS: (iters: %6d)\n",n_iter);
			fprintf(f_time,"# %s %8s %7s\n",    "Stage        ","time(s)","time(%)");
			fprintf(f_time,"# %s %8.2f %7.2f\n","Aligning     ",(float)t_ali, 100*(float)t_ali/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","NMA:         ",(float)t_nma, 100*(float)t_nma/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Springs   ",(float)t_springs, 100*(float)t_springs/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Hessian   ",(float)t_hess, 100*(float)t_hess/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Kinetic   ",(float)t_kine, 100*(float)t_kine/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Diag.     ",(float)t_diag, 100*(float)t_diag/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","After-Diag.  ",(float)t_after, 100*(float)t_after/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","Moving:      ",(float)t_mov, 100*(float)t_mov/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Trial_MDL ",(float)t_trial, 100*(float)t_trial/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Selecting ",(float)t_select, 100*(float)t_select/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Merging   ",(float)t_merge, 100*(float)t_merge/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Mov.NM.   ",(float)t_mov_nm, 100*(float)t_mov_nm/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Mov.Align ",(float)t_mov_ali, 100*(float)t_mov_ali/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","Scoring      ",(float)t_score, 100*(float)t_score/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Project   ",(float)t_mov_proj, 100*(float)t_mov_proj/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","   Corr.     ",(float)t_corr, 100*(float)t_corr/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","Next-stuff   ",(float)t_next, 100*(float)t_next/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f\n","SUM          ",(float)t_sum, 100*(float)t_sum/t_main);
			fprintf(f_time,"# %s %8.2f %7.2f (%6.3f h)\n","Main Loop    ",(float)t_main, 100*(float)t_main/t_main,(float)t_main/3600);

			fclose(f_time);
			printf("%s> Timings file:          %35s\n",prog_name,text);
		}

	sprintf(text,"%s.log",name);
	printf("%s> Log file:              %35s\n",prog_name,text);

	printf("%s>\n%s> Success! Time elapsed %s\n",prog_name,prog_name,time_main.print_time());
	printf("%s> Bye!\n",prog_name);

	//	avg_dmax /=n_iter;
	//	printf("TESTING> dmax= %f avg_dmax= %f  iter_dmax= %d n_iter= %d\n",max_dmax,avg_dmax,iter_dmax, n_iter);
}

void parseOptions(int argc, char** argv)
{
	prog_name = "nmafit";
	prog_name2 = "NMAFIT";
	prog_url = "http://sbg.cib.csic.es/Software/NMAFIT";
	printf("%s>\n%s> Welcome to %s %s\n%s>\n",prog_name,prog_name,prog_name2,VERSION_NMAFIT,prog_name);

	std::string temp;
	CmdLine cmd(prog_name,prog_url, VERSION_NMAFIT );

	try {
		// Define required arguments (not labeled)
		UnlabeledValueArg<std::string> Pdb("pdb","PDB initial model","default","pdb");
		cmd.add( Pdb );
		UnlabeledValueArg<std::string> Ref("map","Target EM map file","default","map");
		cmd.add( Ref );
		UnlabeledValueArg<float> Res("resolution", "Resolution",10,"resolution"); // Sets the resolution (EMAN-like). It should be similar to target map resolution.
		cmd.add( Res );
		UnlabeledValueArg<float> RefThr("cutoff", "EM density map threshold",0.001,"cutoff");
		cmd.add( RefThr );

		SwitchArg Dump("", "dump","Enable dump verbose (default=disabled)", true);
		cmd.add( Dump );
		SwitchArg Timer("", "time","Enable clocks (default=disabled)", true);
		cmd.add( Timer );
		SwitchArg FastProj("", "fast_proj","Enable fast pdb to map projection method (not-trilinear) (default=disabled)", true);
		cmd.add( FastProj );
		SwitchArg NoModel("", "nomodel","Disables PDB model building. Warning: introduced PDB model should match the selected CG (-m option). (default=disabled)", false);
		cmd.add( NoModel );
		//		ValueArg<int> MinNevs("","min_nevs", "Increases nevs value by min_nevs (default=10)",false,10,"int");
		//		cmd.add( MinNevs );
		SwitchArg NoRandWeight("", "norand_weight","Disables random excited modes weighting (default=disabled)", true);
		cmd.add( NoRandWeight );
		SwitchArg RandExcited("", "rand_excited","Random number of excited modes between nex and nex2 (default=disabled)", false);
		cmd.add( RandExcited );
		ValueArg<float> Excited2("","nex2", "Final number of excited eigenvectors, either number [1,nevs] <integer>, or ratio [0,1) <float>. Number of excited modes will change lineally from \"nevs\" to \"nevs2\" (default=disabled)",false,0.2,"int");
		cmd.add( Excited2 );
		SwitchArg RandStep("", "rand_step","Random amplitude selection between step and step2 (default=disabled)", false);
		cmd.add( RandStep );
		ValueArg<double> StepEnd("","step2", "Final amplitude applied to excited modes (amplitude will decrease lineally from \"step\" to \"step2\") (default=step/10)",false,1,"double");
		cmd.add( StepEnd );
		ValueArg<double> Step("","step", "Initial amplitude applied to excited modes (default=0.2)",false,0.2,"double");
		cmd.add( Step );

		ValueArg<int> Verbose("","verb", "Verbosity display level (0=low, 1=medium, 2=high) (default=0)",false,0,"int");
		cmd.add( Verbose );
		ValueArg<unsigned int> Seed("","seed", "Pre-define the random number generator SEED (Mersenne Twister) (default=random-seed from /dev/urandom)",false,386,"unsigned int");
		cmd.add( Seed );

		ValueArg<int> Conv_Length("","conv_win", "Window length (number of iterations) to estimate convergence. (default=5000)",false,5000,"int");
		cmd.add( Conv_Length );
		ValueArg<double> Conver("","conv", "Convergence threshold (relative to first score) (default=1e-8)",false,1e-8,"double");
		cmd.add( Conver );
		ValueArg<float> WRmsd("", "wrmsd","Sets the RMSD weighting [0:1] during model alignment (wrmsd=0, means maximum weighting) (default=0)",false,0.0,"float");
		cmd.add( WRmsd );
		ValueArg<float> DeltaSave("","delta_save", "RMSD increment from previous saved frame to trigger frame saving. (default=0.5)",false,0.5,"float");
		cmd.add( DeltaSave );
		ValueArg<double> Rediag("","rediag", "RMSD ratio to trigger diagonalization. (last_diag_RMSD-current_RMSD)/first_RMSD > rediag (default=0.05)",false,0.05,"double");
		cmd.add( Rediag );

		// Selection probability related
		//		ValueArg<float> FuncX0("","p_x0", "Prob's p_x0 factor (default=0.1)",false,0.1,"float");
		//		cmd.add( FuncX0 );
		//		ValueArg<float> FuncSigma("","p_s", "Prob's p_s factor (default=0.05)",false,0.05,"float");
		//		cmd.add( FuncSigma );
		//		ValueArg<float> MinProb("","min_prob", "Lineal normalization of the selection probability profile between min_prob and 1.0 (default=disabled)",false,0.0,"float");
		//		cmd.add( MinProb );
		ValueArg<std::string> NM_prob("","prob", "Normal Mode Selection Probabitity. "
				"--nex modes will be selected and merged from --nevs subset according to "
				"the following probabilities (p): (default=var)\n"
				"\tplain: equiprobability (--min_prob will be ignored)\n"
				"\tvar: proportional to the mode variance, p(i)= 1/eigenvalue(i)\n"
				"\tline: lineally decreasing probability, p(i)= 1-i/nevs",false,"var","string");
		//				"\tgaussian: p(i)= exp(-(((i-nevs*p_x0)/(nevs*sqrt(2)*p_s))^2)/2)\n"
		//				"\tinvexp: p(i)= 1/(1+nevs*p_x0/i)^p_s )\n"
		//				"Note --min_prob option will normalize the probability.",false,"var","string"); // (plain,var,out,line)
		cmd.add( NM_prob );

		// Pose refinement related
		ValueArg<float> ShakeRot("","6Dref_rot", "Maximum rotational increment (degrees) for pose refinement (default=1)",false,1,"float");
		cmd.add( ShakeRot );
		ValueArg<float> ShakeMov("","6Dref_trans", "Maximum translational increment (Angstroms) for pose refinement. (default=2)",false,2,"float");
		cmd.add( ShakeMov );
		SwitchArg PoseChain("", "6Dref_chain","Independent local 6D pose refinement for each single chain.", false);
		cmd.add( PoseChain );
		ValueArg<int> PoseDelay("","6Dref_delay", "Number of iterations before first local 6D pose refinement (default=disabled)",false,-1,"int");
		cmd.add( PoseDelay );
		ValueArg<int> RefEach("","6Dref", "Number of iterations between local 6D pose refinement (default=200)",false,200,"int"); // Rotational & Traslational refinement each \"refine_each\" iters
		cmd.add( RefEach );

		// Less important parameters (string)
		ValueArg<int> Filter("","filter", "Select filtration method: 1-Fourier, 2-Kernel (By default the fastest will be selected)",false,0,"int");
		cmd.add( Filter );
		ValueArg<float> ModelThr("","cutoff2", "Density cutoff of simulated map (default = 0.001)",false,0.001,"float");
		cmd.add( ModelThr );
		//		SwitchArg NoNormEvec("", "nonormevec","Disables eigenvector \"norm=1\" normalization (default=disabled)", false);
		//		cmd.add( NoNormEvec );
		SwitchArg NoWRmsd("", "nowrmsd","Disables Gaussian Weighted RMSD (default=disabled)", false);
		cmd.add( NoWRmsd );
		SwitchArg NoTors("","notors", "Disables extra torsional potential (default=disabled)", true);
		cmd.add( NoTors );
		SwitchArg UnitMass("","unitmass", "Sets unit masses for every atom (default=disabled)", true);
		cmd.add( UnitMass );
		SwitchArg UnitNElec("","unitnelec", "Sets unit electron density for every atom (default=disabled)", true);
		cmd.add( UnitNElec );
		SwitchArg Rmsds("","morermsds", "Enables C-alpha RMSD computation (default=disabled)", true);
		cmd.add( Rmsds );
		SwitchArg PdbsOut("","morepdbs", "Saves initial (basename_model.pdb) and final (basename_fitted.pdb) models. With -F option, also the fitted CG-model (basename_fitCG.pdb). (default=disabled)", true);
		cmd.add( PdbsOut );
		//		ValueArg<std::string> Funcfile("","ss_func", "Force constants function file based on Topology and Secondary Structure.",false,"","string");
		//		cmd.add( Funcfile );
		//		ValueArg<std::string> SSfile("","ss", "Secondary Structure file with the same format as dssp2ss.pl perl script.",false,"","string");
		//		cmd.add( SSfile );
		ValueArg<std::string> PdbRef("","pdb_ref", "Reference PDB file (test only)",false,"","string");
		cmd.add( PdbRef ); // Check this if PDB-PDB morphing program...

		ValueArg<float> k2_Cutoff("","k2_c","Non-bonding distance cutoff applied to --func option (default=10A).", false, 10,"float");
		cmd.add( k2_Cutoff );
		ValueArg<float> k1_Cte("", "k1_k","Distance cutoff method stiffness constant (default=1.0)",false,1.0,"float");
		cmd.add( k1_Cte );
		ValueArg<float> k1_Cutoff("","k1_c","Distance cutoff method distance cutoff (default=10A)", false, 10,"float");
		cmd.add( k1_Cutoff );
		ValueArg<float> k0_Power("","k0_p", "Inverse Exponential's power term (default=6)",false,6.0,"float");
		cmd.add( k0_Power);
		ValueArg<float> k0_X0("", "k0_x0","Inverse Exponential's inflexion point (default=3.8A)",false,3.8,"float");
		cmd.add( k0_X0 );
		ValueArg<float> k0_Cte("", "k0_k","Inverse Exponential's stiffness constant (default=1.0)",false,1.0,"float");
		cmd.add( k0_Cte );
		ValueArg<float> k0_Cutoff("","k0_c","Inverse Exponential's distance cutoff (default=10A)", false, 10,"float");
		cmd.add( k0_Cutoff );

		ValueArg<std::string> Funcfile("","func", "ASCII file defining the force constant functions to be applied "
				"according to Topology and/or Secondary Structure. "
				"The 5 cols. format is: <SS> <j-i> <k> <x0> <pow>\n"
				"Where <SS> is the two character pairwise interaction identifier, <j-i> is the topology, and "
				"<k>,<x0>,<pow> are the corresponding inverse exponential function parameters. If --ss "
				"is not specified, only the XX pairwise interaction identifier will be considered. "
				"If <j-i> is 0, any previously not-matched topology will be matched.",false,"","string");
		cmd.add( Funcfile );
		ValueArg<std::string> SSfile("","ss", "Secondary Structure ASCII file with 2 cols.: <index> <char>\n"
				"Where <index> is the corresponding residue index (1,2,...), and <char> is the "
				"single character SS identifier. By default SS will be computed internally (H=helix, E=strand, C=coil).",false,"","string");
		cmd.add( SSfile );

		// More important parameters (one letter)
		ValueArg<std::string> FixSS("S","fixSS", "All dihedral coordinates belonging to residues matching the indicated SS identifieres will be fixed. "
				"Ex: \"HE\" will fix the dihedrals corresponding to alpha-helices and beta-sheets.",false,"","string");
		cmd.add( FixSS );
		ValueArg<float> FixRand2("R","fixRand2", "Randomly fixed ratio of Internal Coordinates. Example: 0.7 = 70% of IC will be randomly fixed.",false,0.5,"float");
		cmd.add( FixRand2 );
		ValueArg<float> FixRand("r","fixRand", "Randomly fixed ratio of Dihedral Coordinates. Example: 0.7 = 70% of dihedrals will be randomly fixed. (Rotational/Translational coords. always mobile)",false,0.5,"float");
		cmd.add( FixRand );
		ValueArg<std::string> FixFile("f","fixFile", "ASCII file defining the ICs to be fixed with the format:\n"
				"Protein:     \"index phi chi psi\"\n"
				"NAcid:       \"index alpha beta gamma chi epsilon zeta\"\n"
				"Inter-chain: \"index IC\"\n"
				"Where IC-name= 0(fixed) or 1(mobile)\n"
				"A demo file can be generated using iMODE with the --save_fixfile option.",false,"fixstring","string");
		cmd.add( FixFile );
		ValueArg<int> Contact("P","potential", "Pairwise interaction potential: (default=0)\n"
				"  0= Sigmoid function (= k/(1+(x/x0)^p), if x < c, else k=0)\n"
				"  1= Tirion's cutoff (= k, if x < c, else k=0)\n"
				"  2= Hinsen's function\n"
				"  3= Topology & Secondary Structure (--func is mandatory)\n"
				"  By default an extra torsional potential will be added.",false,0,"int");
		cmd.add( Contact );
		ValueArg<std::string> Name("o","name", "Output files basename. (default=nmafit)",false,"nmafit","string");
		cmd.add( Name );
		SwitchArg Full("F", "full","Enables full-atom output for movie and fitted model.", false);
		cmd.add( Full );
		ValueArg<float> Excited("e","nex", "Excited modes range, either number [1,nevs] <integer>, or ratio [0,1) <float>. (default=0.2)",false,0.2,"int/float");
		cmd.add( Excited );
		ValueArg<float> Nevs("n","nevs", "Used modes range, either number [1,N] <integer>, or ratio [0,1) <float>. (default=0.1)",false,0.1,"int/float");
		cmd.add( Nevs );
		SwitchArg Movie("t", "otraj","Outputs a Multi-PDB trajectory movie (basename_movie.pdb)", false);
		cmd.add( Movie );
		SwitchArg Chi("x","chi", "Considers first CHI dihedral angle. (default=disabled)", true);
		cmd.add( Chi );
		ValueArg<int> Model("m","model", "Coarse-Grained model: 0=CA, 1=3BB2R, 2=Full-Atom, 3=NCAC(experimental) (default=2)",false,2,"int");
		cmd.add( Model );
		SwitchArg Morph("p", "morph","PDB-PDB morphing (introduce a PDB instead of a target <map>)", true); // , and any value at <resolution> and <cutoff>
		cmd.add( Morph ); // Remove it if PDB-PDB morphing program...
		ValueArg<int> Iter("i","iter", "Maximum number of iterations (default=10000)",false,10000,"int");
		cmd.add( Iter );

		// Parse the command line.
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		// =============================================================================
		strcpy(file_initial,((temp=Pdb.getValue()).c_str())); // Gets Initial PDB file name
		if(parse_verb)
			printf("Parser input: Initial PDB file: %s\n",file_initial);

		// Setting Coarse-Graining models
		// Setting model and chi
		model = Model.getValue();
		if(Chi.isSet())
			type = 2; // phi,chi,psi
		else
			type = 0; // phi,psi

		if(FixFile.isSet())
		{
			strcpy(fix_file,((temp=FixFile.getValue()).c_str())); // Gets Fix-file
			fixmodel = 3;
		}
		if(FixRand.isSet())
		{
			fixRand_prob = FixRand.getValue();
			fixmodel = 5; // only dihedral coordinates would be fixed
		}
		if(FixRand2.isSet())
		{
			fixRand_prob = FixRand.getValue();
			fixmodel = 4; // all internal coordinates would be fixed
		}

		if(FixSS.isSet())
		{
			strcpy(fix_ss,((temp=FixSS.getValue()).c_str())); // Gets Fix-string
			fixmodel = 6;
		}

		if(UnitMass.isSet())
		{
			unitmass_switch = true;
			printf("Parser input: Masses will be set to 1.0.\n");
		}

		if(UnitNElec.isSet())
		{
			unitnelec_switch = true;
			printf("Parser input: Electron density per atom will be set to 1.0.\n");
		}

		if( NoRandWeight.isSet() )
			randweight = false;
		printf("Parser input: randweight= %s\n",printbool(randweight));

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
			printf("Parser input: Secondary Structure File: --ss = %s\n",file_ss);
		}
		if(Funcfile.isSet())
		{
			contacts = 3; // override "Contact"
			func_switch = true;
			strcpy(file_func,((temp=Funcfile.getValue()).c_str()));
			printf("Parser input: Secondary Structure and Topology Functions File: --func = %s\n",file_func);
		}

		// Contacting method parameters
		power = k0_Power.getValue();
		cte_k0 = k0_Cte.getValue();
		x0 = k0_X0.getValue();
		cutoff_k0 = k0_Cutoff.getValue();
		cte_k1 = k1_Cte.getValue();
		cutoff_k1 = k1_Cutoff.getValue();
		cutoff_k2 = k2_Cutoff.getValue();

		if(NoTors.isSet())
		{
			notors_switch = true;
			printf("Parser input: Aditional Torsional Springs DISABLED!\n");
		}

		// Selecting Target input (PDB-PDB or PDB-VOL)
		if( Morph.isSet() ) // PDB-PDB
		{
			pdb_switch = true;
			strcpy(file_ref_pdb,((temp=Ref.getValue()).c_str())); // Gets Target PDB file name
			if(parse_verb)
				printf("Parser input: PDB Input Reference file: %s (PDB mode)\n",file_ref_pdb);
		}
		else // PDB-VOL
		{
			pdb_switch = false;
			strcpy(file_ref_vol,((temp=Ref.getValue()).c_str())); // Gets Target Map file name
			savefitted_switch = true; // fitted pdb should be always saved
			if(parse_verb)
				printf("Parser input: VOLUME Input Reference file: %s (VOLUME mode)\n",file_ref_vol);
		}

		filter = Filter.getValue();
		if(parse_verb)
			printf("Parser input: Filtration method: %d (Map mode)\n",filter);

		ref_thr = RefThr.getValue();
		if(parse_verb)
			printf("Parser input: Target Map threshold: %f (Map mode)\n",ref_thr);

		model_thr = ModelThr.getValue();
		if(parse_verb)
			printf("Parser input: Model Map threshold: %f (Map mode)\n",model_thr);

		resolution = Res.getValue();
		if(parse_verb)
			printf("Parser input: Resolution set to: %f\n",resolution);

		strcpy(name,((temp=Name.getValue()).c_str())); // Gets Basename
		if(parse_verb)
			printf("Parser input: Output files Base-Name: %s\n",name);

		convergence = Conver.getValue();
		if(parse_verb)
			printf("Parser input: Convergence cut-off (relative to 1st score): %e\n",convergence);

		// Set NM selection method
		strcpy(text,((temp=NM_prob.getValue()).c_str()));
		if( strcmp(text,"plain") == 0 )
		{
			if(parse_verb)
				printf("Parser input: PLAIN selection probability.\n");
			prob_method = 1;
		}
		else if( strcmp(text,"var") == 0 )
		{
			if(parse_verb)
				printf("Parser input: VAR selection probability.\n");
			prob_method = 2;
		}
		else if( strcmp(text,"line") == 0 )
		{
			prob_method = 4;
			if(parse_verb)
				printf("Parser input: LINE selection probability.\n");
		}
		else if( strcmp(text,"gaussian") == 0 )
		{
			prob_method = 5;
			if(parse_verb)
				printf("Parser input: GAUSSIAN selection probability.\n");
		}
		else if( strcmp(text,"invexp") == 0 )
		{
			prob_method = 6;
			if(parse_verb)
				printf("Parser input: INVEXP selection probability.\n");
		}
		else // Not selected
		{
			printf("\nSorry, you must specify a Normal Mode probability selection method!\nForcing Exit!\n\n");
			exit(1);
		}

		// Number of eigenvectors to be computed 1
		nevec_fact = Nevs.getValue();
		if(nevec_fact <= 0) // checking
		{
			printf("nmafit> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
			exit(1);
		}
		if(parse_verb)
			printf("Parser input: Ratio of Eigenvectors: --nevs = %8.6f\n",nevec_fact);

		nex_fact = Excited.getValue();
		if(nex_fact <= 0) // checking
		{
			printf("nmafit> Error, invalid number excited 1 eigenvectors requested (%f)!\nForcing exit!\n",nex_fact);
			exit(1);
		}
		if(parse_verb)
			printf("Parser input: Ratio of Excited Modes: --nex = %5f\n",nex_fact);

		if(Excited2.isSet())
		{
			nex_fact2 = Excited2.getValue();
			if(nex_fact2 <= 0) // checking
			{
				printf("nmafit> Error, invalid number excited 2 eigenvectors requested (%f)!\nForcing exit!\n",nex_fact2);
				exit(1);
			}
		}
		else
			nex_fact2 = Excited.getValue();
		if(parse_verb)
			printf("Parser input: Ratio of Excited Modes: --nex2 = %5f\n",nex_fact2);


		max_iter = Iter.getValue(); // Gets number of iters
		if(parse_verb)
			printf("Parser input: Maximum Number of Iterations: --iter = %d\n",max_iter);

		refine_each = RefEach.getValue(); // for both, PDB-PDB and PDB-VOL
		if(parse_verb)
			if(pdb_switch) // PDB-PDB
				printf("Parser input: PDB-PDB alignment (minRmsd) each X iters: --6Dref = %d\n",refine_each);
			else // PDB-VOL
				printf("Parser input: Rot.Tras-lational refinement each X iters: --6Dref = %d\n",refine_each);

		//		if(Seed.isSet()) // Fixed seed
		//			seed = (unsigned int) Seed.getValue(); // Gets seed
		//		else // Random seed (time initialization)
		//			seed = (unsigned int) time(0); // int32 seed = (int32)time(0); // random seed (Mersenne.h)
		//		if(parse_verb)
		//			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);
		if(Seed.isSet()) // Fixed seed
			seed = (unsigned int) Seed.getValue(); // Gets seed
		else // Random seed (time initialization)
		{
			// seed = (unsigned int) time(0); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)

			// Mon: /dev/urandom does not work in Windows...
#if defined  _WIN32 || defined _WIN64 // Windows version
			seed = (unsigned int) clock(); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)
#else // Linux version
			// Needed to avoid seed repetition between different runs.
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
#endif
		}

		if(parse_verb)
			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);

		if(Movie.isSet())
		{
			movie_switch = true;
			if(parse_verb)
				printf("Parser input: \"--movie\" = ENABLED\n");
		}

		verbose=Verbose.getValue();
		if(parse_verb)
			printf("Parser input: \"--verb\" = %d\n",verbose);

		if(PdbRef.isSet())
		{
			bench_switch = true; // To allow RMSDs calculation when reference PDB available! (for benchmark use)
			strcpy(file_ref_pdb,((temp=PdbRef.getValue()).c_str())); // Gets Reference PDB file name
			printf("Parser input: Benchmark reference PDB file: %s\n",file_ref_pdb);
		}

		conv_length = Conv_Length.getValue();
		if(parse_verb)
			printf("Parser input: Averaging length (#steps) to reach convergence --conv_length= %d\n",conv_length);

		//		if(MinProb.isSet())
		//		{
		//			min_prob = MinProb.getValue();
		//			if(min_prob<=1 && min_prob>=0)
		//				printf("Parser input: \"--min_prob\" Minimum probability value: %f\n",min_prob);
		//			else
		//			{
		//				printf("Msg(Parser): Please select a valid \"--min_prob\" value!!! [0:1]\nForcing Exit!\n");
		//				exit(1);
		//			}
		//		}

		step = Step.getValue(); // Gets step
		step_begin = step;
		step_end = step;
		if(parse_verb)
			printf("Parser input: Motion --step = %f\n",step);

		if( StepEnd.isSet() )
		{
			rand_step_switch = false;
			step_end = StepEnd.getValue(); // Gets step-end
		}
		else
			step_end = Step.getValue() / 10; // Set default step-end
		if(parse_verb)
			printf("Parser input: Some kind of annealing --step_end = %f\n",step_end);

		if( RandStep.isSet() )
			rand_step_switch = true;

		if( RandExcited.isSet() )
			rand_excited_switch = true;

		if(DeltaSave.isSet())
		{
			delta_save = DeltaSave.getValue(); // Gets Increment to save frames
			if(delta_save<0)
				deltasave_rmsd_switch=false; // non-compatible options
			printf("Parser input: Normalized Score Increment to save frames --delta_save = %f\n",delta_save);
		}

		if( NoWRmsd.isSet() )
		{
			wrmsd_switch = false;
			printf("Parser input: Weighted RMSD PDB alignment disabled\n");
		}
		if( WRmsd.isSet() )
		{
			wrmsd_switch = true;
			min_wrmsd = WRmsd.getValue();
			printf("Parser input: Weighted RMSD PDB alignment enabled, weight= %f\n",min_wrmsd);
		}
		//		if( NoNormEvec.isSet() )
		//		{
		//			norm_evec_switch = false;
		//			printf("Parser input: Eigenvector normalization disabled\n");
		//		}
		if( Rmsds.isSet() )
		{
			more_rmsds = true;
			if(parse_verb)
				printf("Parser input: C-alpha RMSDs enabled\n");
			if(!bench_switch && !pdb_switch) // if PDB-VOL && not benchmark
			{
				printf("Parser Error: Please, include a --pdb_ref file to compute RMSDs...\n");
				exit(3);
			}
		}
		if( PdbsOut.isSet() )
		{
			savemodel_switch = true;
			savefitted_switch = true;
			if(parse_verb)
				printf("Parser input: Initial/final model pdbs output enabled\n");
		}

		if(Rediag.isSet())
		{
			rediag = Rediag.getValue(); // Gets step
			printf("Parser input: Diagonalization threshold (score increment): --rediag = %f\n",rediag);
		}

		if(ShakeRot.isSet())
		{
			shake_rot = ShakeRot.getValue() * (M_PI/180); // radians
			printf("Parser input: Max. Rotational Shake: --shake_rot = %f\n",shake_rot);
		}

		if(ShakeMov.isSet())
		{
			shake_mov = ShakeMov.getValue();
			printf("Parser input: Max. Traslational Shake: --shake_mov = %f\n",shake_mov);
		}

		if(Excited.isSet())
			nex = Excited.getValue();

		if( PoseChain.isSet() )
			sel_chain = 0; // A different chain will be selected each 6D pose ref.
		else
			sel_chain = -1; // the whole molecule will be selected

		refine_delay = PoseDelay.getValue(); // 6D pose refinement delay (number of iterations)
		printf("Parser input: 6D pose refinement number of iterations delay: --6Dref_delay = %d\n",refine_delay);

		fullatom_switch = Full.isSet();
		nomodel = NoModel.isSet();
		fast_switch = FastProj.isSet();
		time_switch = Timer.isSet();
		dump_switch = Dump.isSet();

		//		min_nevs = MinNevs.getValue();
		//		printf("Parser input: Minimum number of eigenvectors: --min_nevs = %d\n",min_nevs);

		//		p_x0 = FuncX0.getValue();
		//		printf("Parser input: Probability method's function X0 factor: --p_x0 = %f\n",p_x0);
		//
		//		p_s = FuncSigma.getValue();
		//		printf("Parser input: Probability method's function sigma factor: --p_s = %f\n",p_s);

	} catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}

void parseOptionsMorph(int argc, char** argv)
{
	std::string temp;
	prog_name = "imorph";
	prog_name2 = "iMORPH";
	prog_url = "iMORPH: Internal coordinates MORPHing tool.";
	printf("%s>\n%s> Welcome to %s %s\n%s>\n",prog_name,prog_name,prog_name2,VERSION_NMAFIT,prog_name);

	CmdLine cmd(prog_name,prog_url, VERSION_NMAFIT );


	try {
		// Define required arguments (not labeled)
		UnlabeledValueArg<std::string> Pdb("pdb1","Initial PDB file.","default","initial_pdb");
		cmd.add( Pdb );
		UnlabeledValueArg<std::string> Ref("pdb2","Target PDB file.","default","target_pdb");
		cmd.add( Ref );

		//		SwitchArg Dump("", "dump","Enable dump stuff (default=disabled).", true);
		//		cmd.add( Dump );

		//		// Gaussian-scoring stuff... EXPERIMENTAL (DON'T DELETE!)
		//      ValueArg<float> Gauss_bound("", "bound_factor","Set bounding box factor (EXPERIMENTAL) (default=disabled).",false,5.0,"float");
		//		cmd.add( Gauss_bound );
		//      ValueArg<float> Gauss("", "gauss","Set resolution for Gaussian overlap scoring (EXPERIMENTAL) (default=disabled).",false,10.0,"float");
		//		cmd.add( Gauss );

		ValueArg<int> Verbose("","verb", "Verbose level (0=low, 1=medium, 2=high) (default=0).",false,0,"int");
		cmd.add( Verbose );
		SwitchArg Server("","server", "Enable \"server mode\" to maximize the automatic selection of parameters. (default=disabled)", false);
		cmd.add( Server );
		SwitchArg DelHeteros("","delete_heteros", "Delete Hetero-atoms, including waters (default=disabled).", true);
		cmd.add( DelHeteros );
		SwitchArg KeepWaters("","keep_waters", "Disables Water molecules deletion (default=disabled).", true);
		cmd.add( KeepWaters );
		SwitchArg KeepHydrogens("","keep_hydrogens", "Disables Hydrogen atoms deletion (default=disabled).", true);
		cmd.add( KeepHydrogens );
		ValueArg<int> AlnLevel("","aln_level", "Clustal alignment level (1, 2 or 3 to match *, *: or *:. ,respectively) (default=2)",false,2,"int");
		cmd.add( AlnLevel );
		ValueArg<double> WrmsdC("","gauss_c", "Distance parameter for gaussian weighting (wrmsd) (default=20).",false,20,"double");
		cmd.add( WrmsdC );
		ValueArg<int> PoseDelay("","6Dref_delay", "Number of iterations before first local 6D pose refinement (default=disabled)",false,-1,"int");
		cmd.add( PoseDelay );
		ValueArg<std::string> AlnFile("a","aln_file", "Clustal's \".aln\" file defining alingment between initial and target atomic structures.\n",false,"alnstring","string");
		cmd.add( AlnFile );
		ValueArg<int> Eigensolver("","eigensolver", "Eigensolver (fastest option will be selected by default):\n"
				"  0= LAPACK/BLAS, fastest for >6% modes [DSPGVX],\n"
				"  1= ARPACK, fastest for <=6% modes [dsdrv1_AP_BP_W_mon] (default),\n"
				"  2= ARPACK-square, [dsdrv1_A_B_W_mon] (experimental).",false,1,"int");
		cmd.add( Eigensolver );
		SwitchArg Timer("", "time","Enable clocks (default=disabled).", true);
		cmd.add( Timer );
		ValueArg<unsigned int> Seed("","seed", "Set the random number generator seed (Mersenne Twister) (default=random-seed from /dev/urandom).",false,386,"unsigned int");
		cmd.add( Seed );
		ValueArg<float> WRmsd("", "wrmsd","Sets the RMSD weighting factor for model alignment. "
				"The range is [0:1], 0=maximum-weighting 1=no-weighting (default=0).",false,0.0,"float");
		cmd.add( WRmsd );
		SwitchArg NoRandWeight("", "norand_weight","Disables random excited modes weighting (default=disabled).", true);
		cmd.add( NoRandWeight );
		SwitchArg SCVWeight("", "noscv","Disable the Scaling Collective Variable method (default=disabled).", true);
		cmd.add( SCVWeight );
		SwitchArg RandExcited("", "rand_excited","Random number of excited modes between --nex and --nex2 (default=disabled).", false);
		cmd.add( RandExcited );
		ValueArg<float> Excited2("","nex2", "Final number of excited eigenvectors, either number [1,nevs] <integer>, or ratio [0,1) <float>. "
				"Number of excited modes will change linearly from --nevs to --nevs2 (default=disabled).",false,0.2,"int");
		cmd.add( Excited2 );
		SwitchArg RandStep("", "rand_step","Random amplitude selection between --step and --step2 (default=disabled).", false);
		cmd.add( RandStep );
		ValueArg<double> StepEnd("","step2", "Final amplitude applied to excited modes (amplitude will decrease linearly from --step to --step2) (default=step/10).",false,1,"double");
		cmd.add( StepEnd );
		ValueArg<int> AddNevs("","addnevs", "Increases --nevs value by --addnevs (default=10).",false,10,"int");
		cmd.add( AddNevs );
		ValueArg<std::string> NM_prob("","prob", "Normal mode selection probabitity. "
				"Each one of the --nex modes will be selected and merged from the --nevs subset according to "
				"the following probabilities (p): (default=var)\n"
				"\tplain: constant probability\n"
				"\tvar: proportional to the i-th mode variance, p(i)= 1/eigenvalue(i)\n"
				"\tline: lineally decreasing probability, p(i)= 1-i/nevs",false,"var","string");
		cmd.add( NM_prob );

		// ONLY FOR FITTING?
		//		// Pose refinement related
		//		ValueArg<float> ShakeRot("","6Dref_rot", "Maximum rotational increment (degrees) for pose refinement (default=1)",false,1,"float");
		//		cmd.add( ShakeRot );
		//		ValueArg<float> ShakeMov("","6Dref_trans", "Maximum translational increment (Angstroms) for pose refinement. (default=2)",false,2,"float");
		//		cmd.add( ShakeMov );
		//		SwitchArg PoseChain("", "6Dref_chain","Independent local 6D pose refinement for each single chain.", false);
		//		cmd.add( PoseChain );
		//		ValueArg<int> PoseDelay("","6Dref_delay", "Number of iterations before first local 6D pose refinement (default=disabled)",false,-1,"int");
		//		cmd.add( PoseDelay );
		ValueArg<int> RefEach("","6Dref", "Number of iterations between local 6D pose refinement (default=200)",false,200,"int"); // Rotational & Traslational refinement each \"refine_each\" iters
		cmd.add( RefEach );

		// Less important parameters (string)
		SwitchArg NoWRmsd("", "nowrmsd","Disables Gaussian weighted RMSD (default=disabled).", false);
		cmd.add( NoWRmsd );
		SwitchArg NoTors("","notors", "Disables extra torsional potential (default=disabled).", true);
		cmd.add( NoTors );
		SwitchArg UnitMass("","unitmass", "Sets unit masses for every atom (default=disabled)", true);
		cmd.add( UnitMass );
		SwitchArg UnitNElec("","unitnelec", "Sets unit electron density for every atom (default=disabled)", true);
		cmd.add( UnitNElec );
		SwitchArg NoModel("", "nomodel","Disables PDB model building. "
				"Warning: introduced PDB model must match the CG selected with the -m option (default=disabled).", false);
		cmd.add( NoModel );
		SwitchArg NoMovie("", "notraj","Disables Multi-PDB trajectory movie output (default=disabled).", false);
		cmd.add( NoMovie );

		//		ValueArg<std::string> PdbRef("","pdb_ref", "Reference PDB file (test only)",false,"","string");
		//		cmd.add( PdbRef ); // Check this if PDB-PDB morphing program...
		SwitchArg Rmsds("","morermsds", "Enables C-alpha RMSD computation (default=disabled).", true);
		cmd.add( Rmsds );
		SwitchArg PdbsOut("","morepdbs", "Saves initial (basename_model.pdb) and final (basename_fitted.pdb) models. "
				"If the -F option is enabled, the fitted CG-model (basename_fitCG.pdb) will be saved (default=disabled).", true);
		cmd.add( PdbsOut );
		ValueArg<double> Rediag("","rediag", "RMSD to trigger NMA (default=0.1A).",false,0.1,"double");
		cmd.add( Rediag );

		ValueArg<std::string> Funcfile("","func", "ASCII file defining the force constant functions to be applied "
				"according to Topology and/or Secondary Structure. "
				"The 5 cols. format is: <SS> <t> <k> <x0> <pow>\n"
				"Where <SS> is the two character pairwise interaction identifier, <t> is the topology, and "
				"<k>,<x0>,<pow> are the corresponding sigmoid function parameters. If --ss "
				"is not specified, the XX pairwise interaction identifier must be introduced. This way, only topologies will be considered. "
				"If <t> is \"-1\", any previously not-matched topology will be considered.",false,"","string");
		cmd.add( Funcfile );

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

		// Medium importance parameters
		ValueArg<float> DeltaSave("","delta_save", "RMSD increment to save a new trajectory frame (default=0.5A). "
				"If a negative integer value is introduced, a new frame will be saved each --delta_save iterations.",false,0.5,"float");
		cmd.add( DeltaSave );
		ValueArg<int> Conv_Length("","conv_win", "Window length (number of iterations) to estimate convergence (default=1000).",false,1000,"int");
		cmd.add( Conv_Length );
		ValueArg<double> Conver("","conv", "Convergence RMSD threshold (default=0.01A).",false,0.01,"double");
		cmd.add( Conver );
		ValueArg<std::string> SSfile("","ss", "Secondary Structure ASCII file with 2 cols.: <n> <char>\n"
				"Where <n> is the corresponding residue index (0,1,...), and <char> is the "
				"single character SS identifier. "
				"By default SS will be computed internally (H=helix, E=strand, C=coil).",false,"","string");
		cmd.add( SSfile );
		ValueArg<std::string> FixSS("S","fixSS", "All dihedral coordinates with a given secondary structure (SS) "
				"will be removed (see --ss). Ex: \"HE\" will fix the dihedrals corresponding to alpha-helices "
				"and beta-sheets.",false,"","string");
		cmd.add( FixSS );
		ValueArg<float> FixRand2("R","fixRand2", "Randomly fixed ratio of Internal Coordinates (default=disabled). Example: 0.7 = 70% of IC will be randomly fixed.",false,0.5,"float");
		cmd.add( FixRand2 );
		ValueArg<float> Excited("e","nex", "Excited modes range, either number [1,nevs] <integer>, or ratio [0,1) <float>. (default=0.1)",false,0.1,"int/float");
		cmd.add( Excited );
		ValueArg<double> Step("s","step", "Initial amplitude applied to the merge displacement (default=5).",false,5,"double");
		cmd.add( Step );
		SwitchArg Chi("x","chi", "Considers first CHI dihedral angle (default=disabled).", true);
		cmd.add( Chi );
		ValueArg<float> Nevs("n","nevs", "Used modes range, either number [1,N] <integer>, or ratio [0,1) <float> (default=0.1). In any case, the value of --addnevs option will be added.",false,0.1,"int/float");
		cmd.add( Nevs );
		ValueArg<int> Iter("i","iter", "Maximum number of iterations (default=100000).",false,100000,"int");
		cmd.add( Iter );

		// Very important parameters (one letter)
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
				"A demo file can be generated using iMode with the --save_fixfile option.",false,"fixstring","string");
		cmd.add( FixFile );
		ValueArg<float> FixRand("r","fixRand", "Randomly fixed ratio of Dihedral Coordinates (default=disabled). Example: 0.7 = 70% of dihedrals will be randomly removed. Rotational/translational coords. always mobile.",false,0.5,"float");
		cmd.add( FixRand );
		SwitchArg Full("F", "full","Enables full-atom output models (default=disabled).", false);
		cmd.add( Full );
		ValueArg<std::string> Name("o","name", "Output files basename (default=imorph).",false,prog_name,"string");
		cmd.add( Name );
		ValueArg<int> Model("m","model", "Coarse-Grained model: 0=CA, 1=C5, 2=Heavy-Atom (default=2).",false,2,"int");
		cmd.add( Model );

		// Parse the command line.
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		// =============================================================================
		strcpy(file_initial,((temp=Pdb.getValue()).c_str())); // Gets Initial PDB file name
		if(parse_verb)
			printf("Parser input: Initial PDB file: %s\n",file_initial);

		// Setting Coarse-Graining models
		// Setting model and chi
		model = Model.getValue();
		if(Chi.isSet())
			type = 2; // phi,chi,psi
		else
			type = 0; // phi,psi

		if(FixFile.isSet())
		{
			strcpy(fix_file,((temp=FixFile.getValue()).c_str())); // Gets Fix-file
			fixmodel = 3;
		}
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

		if(UnitMass.isSet())
		{
			unitmass_switch = true;
			if(parse_verb)
				printf("Parser input: Masses will be set to 1.0.\n");
		}

		if(UnitNElec.isSet())
		{
			unitnelec_switch = true;
			printf("Parser input: Electron density per atom will be set to 1.0.\n");
		}

		if( NoRandWeight.isSet() )
			randweight = false;
		if(parse_verb)
			printf("Parser input: randweight= %s\n",printbool(randweight));

		if( SCVWeight.isSet() )
			scv_weight = false;
		else
			scv_weight = true; // default enabled
		if(parse_verb)
			printf("Parser input: scv_weight = %s\n",printbool(scv_weight));

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

		// Contacting method parameters
		power = k0_Power.getValue();
		cte_k0 = k0_Cte.getValue();
		x0 = k0_X0.getValue();
		cutoff_k0 = k0_Cutoff.getValue();
		cte_k1 = k1_Cte.getValue();
		cutoff_k1 = k1_Cutoff.getValue();
		cutoff_k2 = k2_Cutoff.getValue();

		if(NoTors.isSet())
		{
			notors_switch = true;
			if(parse_verb)
				printf("Parser input: Aditional Torsional Springs DISABLED!\n");
		}

		//		// Selecting Target input (PDB-PDB or PDB-VOL)
		//		if( Morph.isSet() ) // PDB-PDB
		//		{
		pdb_switch = true;
		strcpy(file_ref_pdb,((temp=Ref.getValue()).c_str())); // Gets Target PDB file name
		if(parse_verb)
			printf("Parser input: Target PDB file: %s (PDB mode)\n",file_ref_pdb);
		//		}
		//		else // PDB-VOL
		//		{
		//			pdb_switch = false;
		//			strcpy(file_ref_vol,((temp=Ref.getValue()).c_str())); // Gets Target Map file name
		//			if(parse_verb)
		//				printf("Parser input: Target VOLUME file: %s (VOLUME mode)\n",file_ref_vol);
		//		}

		//		ref_thr = RefThr.getValue();
		//		if(parse_verb)
		//			printf("Parser input: Target Map threshold: %f (Map mode)\n",ref_thr);

		//		model_thr = ModelThr.getValue();
		//		if(parse_verb)
		//			printf("Parser input: Model Map threshold: %f (Map mode)\n",model_thr);

		//		resolution = Res.getValue();
		//		if(parse_verb)
		//			printf("Parser input: Resolution set to: %f\n",resolution);

		strcpy(name,((temp=Name.getValue()).c_str())); // Gets Basename

		if(AlnFile.isSet())
		{
			strcpy(file_aln,((temp=AlnFile.getValue()).c_str())); // Gets Alignment file name
			aln_switch = true;
			if(parse_verb)
				printf("Parser input: Alignment file (%s) for correspondence between initial and target models ENABLED!\n",file_aln);
		}

		gauss_c = WrmsdC.getValue(); // Distance parameter for gaussian weight (wrmsd)
		if(parse_verb)
			printf("Parser input: Distance parameter for gaussian weight --gauss_c = %f\n",gauss_c);

		aln_level = AlnLevel.getValue();
		if(parse_verb)
			printf("Parser input: Clustal alignment level: --aln_level = %d\n",aln_level);

		refine_delay = PoseDelay.getValue(); // 6D pose refinement delay (number of iterations)
		if(parse_verb)
			printf("Parser input: 6D pose refinement iter. delay: --6Dref_delay = %d\n",refine_delay);

		refine_each = RefEach.getValue(); // for both, PDB-PDB and PDB-VOL
		if(parse_verb)
					printf("Parser input: 6D pose refinement interval: --6Dref = %d\n",refine_each);

		if(parse_verb)
			printf("Parser input: Output files Base-Name: %s\n",name);

		convergence = Conver.getValue();
		if(parse_verb)
			printf("Parser input: Convergence cut-off (relative to 1st score): %e\n",convergence);

		// Set NM selection method
		strcpy(text,((temp=NM_prob.getValue()).c_str()));
		if( strcmp(text,"plain") == 0 )
		{
			if(parse_verb)
				printf("Parser input: PLAIN selection probability.\n");
			prob_method = 1;
		}
		else if( strcmp(text,"var") == 0 )
		{
			if(parse_verb)
				printf("Parser input: VAR selection probability.\n");
			prob_method = 2;
		}
		else if( strcmp(text,"line") == 0 )
		{
			prob_method = 4;
			if(parse_verb)
				printf("Parser input: LINE selection probability.\n");
		}
		else if( strcmp(text,"gaussian") == 0 )
		{
			prob_method = 5;
			if(parse_verb)
				printf("Parser input: GAUSSIAN selection probability.\n");
		}
		else if( strcmp(text,"invexp") == 0 )
		{
			prob_method = 6;
			if(parse_verb)
				printf("Parser input: INVEXP selection probability.\n");
		}
		else // Not selected
		{
			printf("\nSorry, you must specify a Normal Mode probability selection method!\nForcing Exit!\n\n");
			exit(1);
		}

		// Number of eigenvectors to be computed 1
		nevec_fact = Nevs.getValue();
		if(nevec_fact <= 0) // checking
		{
			printf("nmafit> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
			exit(1);
		}
		if(parse_verb)
			printf("Parser input: Ratio of Eigenvectors: --nevs = %8.6f\n",nevec_fact);

		nex_fact = Excited.getValue();
		if(nex_fact <= 0) // checking
		{
			printf("nmafit> Error, invalid number excited 1 eigenvectors requested (%f)!\nForcing exit!\n",nex_fact);
			exit(1);
		}
		if(parse_verb)
			printf("Parser input: Ratio of Excited Modes: --nex = %5f\n",nex_fact);

		if(Excited2.isSet())
		{
			nex_fact2 = Excited2.getValue();
			if(nex_fact2 <= 0) // checking
			{
				printf("nmafit> Error, invalid number excited 2 eigenvectors requested (%f)!\nForcing exit!\n",nex_fact2);
				exit(1);
			}
		}
		else
			nex_fact2 = Excited.getValue();
		if(parse_verb)
			printf("Parser input: Ratio of Excited Modes: --nex2 = %5f\n",nex_fact2);

		addnevs = AddNevs.getValue();
		if(parse_verb)
			printf("Parser input: Number of added modes: --addnevs = %d\n",addnevs);

		max_iter = Iter.getValue(); // Gets number of iters
		if(parse_verb)
			printf("Parser input: Maximum Number of Iterations: --iter = %d\n",max_iter);

		if(DelHeteros.isSet())
			delHeteros_switch = true; // Delete heteroatoms
		if(KeepHydrogens.isSet())
			delHydrogens_switch = false; // Keep hydrogens
		if(KeepWaters.isSet())
			delWaters_switch = false; // Keep waters

		//		refine_each = RefEach.getValue(); // for both, PDB-PDB and PDB-VOL
		//		if(parse_verb)
		//			if(pdb_switch) // PDB-PDB
		//				printf("Parser input: PDB-PDB alignment (minRmsd) each X iters: --6Dref = %d\n",refine_each);
		//			else // PDB-VOL
		//				printf("Parser input: Rot.Tras-lational refinement each X iters: --6Dref = %d\n",refine_each);

		//		if(Seed.isSet()) // Fixed seed
		//			seed = (unsigned int) Seed.getValue(); // Gets seed
		//		else // Random seed (time initialization)
		//			seed = (unsigned int) time(0); // int32 seed = (int32)time(0); // random seed (Mersenne.h)
		//		if(parse_verb)
		//			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);
		if(Seed.isSet()) // Fixed seed
			seed = (unsigned int) Seed.getValue(); // Gets seed
		else // Random seed (time initialization)
		{
			// seed = (unsigned int) time(0); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)

			// Mon: /dev/urandom does not work in Windows...
#if defined  _WIN32 || defined _WIN64 // Windows version
			seed = (unsigned int) clock(); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)
#else // Linux version
			// Needed to avoid seed repetition between different runs.
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
#endif
		}

		if(parse_verb)
			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);

		//		if(Movie.isSet())
		//		{
		//			movie_switch = true;
		//			if(parse_verb)
		//				printf("Parser input: \"--movie\" = ENABLED\n");
		//		}

		if(NoMovie.isSet())
		{
			movie_switch = false;
			if(parse_verb)
				printf("Parser input: \"--notraj\" = ENABLED\n");
		}
		else
			movie_switch = true; // movie output by default

		verbose=Verbose.getValue();
		if(parse_verb)
			printf("Parser input: \"--verb\" = %d\n",verbose);

		//		if(PdbRef.isSet())
		//		{
		//			bench_switch = true; // To allow RMSDs calculation when reference PDB available! (for benchmark use)
		//			strcpy(file_ref_pdb,((temp=PdbRef.getValue()).c_str())); // Gets Reference PDB file name
		//			if(parse_verb)
		//				printf("Parser input: Benchmark reference PDB file: %s\n",file_ref_pdb);
		//		}

		conv_length = Conv_Length.getValue();
		if(parse_verb)
			printf("Parser input: Averaging length (#steps) to reach convergence --conv_length= %d\n",conv_length);

		//		if(MinProb.isSet())
		//		{
		//			min_prob = MinProb.getValue();
		//			if(min_prob<=1 && min_prob>=0)
		//			{
		//				if(parse_verb)
		//					printf("Parser input: \"--min_prob\" Minimum probability value: %f\n",min_prob);
		//			}
		//			else
		//			{
		//				printf("Msg(Parser): Please select a valid \"--min_prob\" value!!! [0:1]\nForcing Exit!\n");
		//				exit(1);
		//			}
		//		}

		step = Step.getValue(); // Gets step
		step_begin = step;
		step_end = step;
		if(parse_verb)
			printf("Parser input: Motion --step = %f\n",step);

		if( StepEnd.isSet() )
		{
			rand_step_switch = false;
			step_end = StepEnd.getValue(); // Gets step-end
		}
		else
			step_end = Step.getValue() / 10; // Set default step-end
		if(parse_verb)
			printf("Parser input: Some kind of annealing --step_end = %f\n",step_end);

		if( RandStep.isSet() )
			rand_step_switch = true;

		if( RandExcited.isSet() )
			rand_excited_switch = true;

		if(DeltaSave.isSet())
		{
			delta_save = DeltaSave.getValue(); // Gets Increment to save frames
			if(delta_save<0)
				deltasave_rmsd_switch=false; // non-compatible options
			if(parse_verb)
				printf("Parser input: Normalized Score Increment to save frames --delta_save = %f\n",delta_save);
		}

		if( NoWRmsd.isSet() )
		{
			wrmsd_switch = false;
			if(parse_verb)
				printf("Parser input: Weighted RMSD PDB alignment disabled\n");
		}
		if( WRmsd.isSet() )
		{
			wrmsd_switch = true;
			min_wrmsd = WRmsd.getValue();
			if(parse_verb)
				printf("Parser input: Weighted RMSD PDB alignment enabled, weight= %f\n",min_wrmsd);
		}
		//		if( NoNormEvec.isSet() )
		//		{
		//			norm_evec_switch = false;
		//			if(parse_verb)
		//				printf("Parser input: Eigenvector normalization disabled\n");
		//		}
		if( Rmsds.isSet() )
		{
			more_rmsds = true;
			if(parse_verb)
				printf("Parser input: C-alpha RMSDs enabled\n");
		}
		if( PdbsOut.isSet() )
		{
			savemodel_switch = true;
			savefitted_switch = true;
			if(parse_verb)
				printf("Parser input: Initial/final model pdbs output enabled\n");
		}

		if(Rediag.isSet())
		{
			rediag = Rediag.getValue(); // Gets step
			if(parse_verb)
				printf("Parser input: Diagonalization threshold (score increment): --rediag = %f\n",rediag);
		}

		//		if(ShakeRot.isSet())
		//		{
		//			shake_rot = ShakeRot.getValue() * (M_PI/180); // radians
		//			if(parse_verb)
		//				printf("Parser input: Max. Rotational Shake: --shake_rot = %f\n",shake_rot);
		//		}
		//
		//		if(ShakeMov.isSet())
		//		{
		//			shake_mov = ShakeMov.getValue();
		//			if(parse_verb)
		//				printf("Parser input: Max. Traslational Shake: --shake_mov = %f\n",shake_mov);
		//		}

		if(Excited.isSet())
			nex = Excited.getValue();

		//		if( PoseChain.isSet() )
		//			sel_chain = 0; // A different chain will be selected each 6D pose ref.
		//		else
		//			sel_chain = -1; // the whole molecule will be selected

		//		refine_delay = PoseDelay.getValue(); // 6D pose refinement delay (number of iterations)
		//		if(parse_verb)
		//			printf("Parser input: 6D pose refinement number of iterations delay: --6Dref_delay = %d\n",refine_delay);

		if(Full.isSet())
			fullatom_switch = true;

		if(Timer.isSet())
			time_switch = true;

		eigensolver = Eigensolver.getValue();
		if(Eigensolver.isSet())
		{
			eigensolver_switch = true; // true --> automatic eigensolver choice will be overridden
			printf("Parser> Using eigensolver: %d (automatic choice will be overridden)\n",eigensolver);
		}
		else
			printf("Parser> Fastest eigensolver will be selected (see below)\n");

		server = Server.getValue(); // Server mode

		// Gaussian-scoring stuff... EXPERIMENTAL (DON'T DELETE!)
		//		if(Gauss.isSet())
		//		{
		//			gauss_switch = true;
		//			resolution = Gauss.getValue();
		//		}
		//		gauss_bound = Gauss_bound.getValue();

		//		if(Dump.isSet())
		//			dump_switch = true;

		if(NoModel.isSet())
			nomodel = true;

		//		min_nevs = MinNevs.getValue();
		//		if(parse_verb)
		//			printf("Parser input: Minimum number of eigenvectors: --min_nevs = %d\n",min_nevs);

		//		p_x0 = FuncX0.getValue();
		//		if(parse_verb)
		//			printf("Parser input: Probability method's function X0 factor: --p_x0 = %f\n",p_x0);
		//
		//		p_s = FuncSigma.getValue();
		//		if(parse_verb)
		//			printf("Parser input: Probability method's function sigma factor: --p_s = %f\n",p_s);

	} catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}

void parseOptionsFit(int argc, char** argv)
{
	std::string temp;
	prog_name = "imodfit";
	prog_name2 = "iMODFIT";
	prog_url = "iMODFIT: Internal coordinates MODes FITting tool.";
	printf("%s>\n%s> Welcome to %s %s\n%s>\n",prog_name,prog_name,prog_name2,VERSION_NMAFIT,prog_name);

	CmdLine cmd(prog_name,prog_url, VERSION_NMAFIT );

	try {
		// Define required arguments (not labeled)
		UnlabeledValueArg<std::string> Pdb("pdb","PDB initial model","default","pdb");
		cmd.add( Pdb );
		UnlabeledValueArg<std::string> Ref("map","Target EM map file","default","map");
		cmd.add( Ref );
		UnlabeledValueArg<float> Res("resolution", "Resolution",10,"resolution"); // Sets the resolution (EMAN-like). It should be similar to target map resolution.
		cmd.add( Res );
		UnlabeledValueArg<float> RefThr("cutoff", "EM density map threshold",0.001,"cutoff");
		cmd.add( RefThr );

		SwitchArg Chimera("", "chimera","Enable output for Chimera (dump current movie frame to <base>_fitted.pdb). (default=disabled)", true);
		cmd.add( Chimera );
		ValueArg<float> InterMolec("", "inter_molec","Sets the inter molecular force constant factor (default=disabled)",false,1.0,"float");
		cmd.add( InterMolec );
		SwitchArg Dump("", "dump","Enable dump verbose (default=disabled)", true);
		cmd.add( Dump );
		ValueArg<int> NThreads("","nthreads", "Number of threads used for parallel processing",false,0,"int");
		cmd.add( NThreads );
		ValueArg<int> Eigensolver("","eigensolver", "Eigensolver (fastest option will be selected by default):\n"
				"  0= LAPACK/BLAS, fastest for >6% modes [DSPGVX],\n"
				"  1= ARPACK, fastest for <=6% modes [dsdrv1_AP_BP_W_mon] (default),\n"
				"  2= ARPACK-square, [dsdrv1_A_B_W_mon] (experimental).",false,1,"int");
		cmd.add( Eigensolver );
		SwitchArg Timer("", "time","Enable clocks (default=disabled)", true);
		cmd.add( Timer );
		SwitchArg FastProj("", "fast_proj","Enable fast pdb to map projection method (not-trilinear) (default=disabled)", true);
		cmd.add( FastProj );
		SwitchArg SaveMaps("", "save_maps","Enable to save target, initial and fitted maps. The latter two are projected and filtered from corresponding atomic models. (default=disabled)", true);
		cmd.add( SaveMaps );
		SwitchArg NoModel("", "nomodel","Disables PDB model building. Warning: introduced PDB model should match the selected CG (-m option). (default=disabled)", false);
		cmd.add( NoModel );
		//		ValueArg<int> MinNevs("","min_nevs", "Increases nevs value by min_nevs (default=10)",false,10,"int");
		//		cmd.add( MinNevs );
		SwitchArg SCVWeight("", "noscv","Disable the Scaling Collective Variable method (default=disabled).", true);
		cmd.add( SCVWeight );
		SwitchArg NoRandWeight("", "norand_weight","Disables random excited modes weighting (default=disabled)", true);
		cmd.add( NoRandWeight );
		SwitchArg RandExcited("", "rand_excited","Random number of excited modes between nex and nex2 (default=disabled)", false);
		cmd.add( RandExcited );
		ValueArg<int> AddNevs("","addnevs", "Increases --nevs value by --addnevs (default=10)",false,10,"int");
		cmd.add( AddNevs );
		ValueArg<float> Excited2("","nex2", "Final number of excited eigenvectors, either number [1,nevs] <integer>, or ratio [0,1) <float>. Number of excited modes will change lineally from \"nevs\" to \"nevs2\" (default=disabled)",false,0.2,"int");
		cmd.add( Excited2 );
		SwitchArg RandStep("", "rand_step","Random amplitude selection between step and step2 (default=disabled)", false);
		cmd.add( RandStep );
		ValueArg<double> StepEnd("","step2", "Final amplitude applied to excited modes (amplitude will decrease lineally from \"step\" to \"step2\") (default=step/10)",false,1,"double");
		cmd.add( StepEnd );
		ValueArg<double> Step("","step", "Initial amplitude applied to excited modes (default=5)",false,5,"double");
		cmd.add( Step );

		ValueArg<int> Verbose("","verb", "Verbosity display level (0=low, 1=medium, 2=high) (default=0)",false,0,"int");
		cmd.add( Verbose );
		ValueArg<unsigned int> Seed("","seed", "Pre-define the random number generator SEED (Mersenne Twister) (default=random-seed from /dev/urandom)",false,386,"unsigned int");
		cmd.add( Seed );
		SwitchArg DelHeteros("","delete_heteros", "Delete Hetero-atoms, including waters (default=disabled).", true);
		cmd.add( DelHeteros );
		SwitchArg KeepWaters("","keep_waters", "Disables Water molecules deletion (default=disabled).", true);
		cmd.add( KeepWaters );
		SwitchArg KeepHydrogens("","keep_hydrogens", "Disables Hydrogen atoms deletion (default=disabled).", true);
		cmd.add( KeepHydrogens );
		ValueArg<int> Conv_Length("","conv_win", "Window length (number of iterations) to estimate convergence. (default=1000)",false,1000,"int");
		cmd.add( Conv_Length );
		ValueArg<double> Conver("","conv", "Convergence cross-correlation threshold (default=0.0001)",false,0.0001,"double");
		cmd.add( Conver );
		ValueArg<float> WRmsd("", "wrmsd","Sets the RMSD weighting [0:1] during model alingment (wrmsd=0, means maximum weighting) (default=0)",false,0.0,"float");
		cmd.add( WRmsd );

		// Selection probability related
		ValueArg<std::string> NM_prob("","prob", "Normal Mode Selection Probabitity. "
				"--nex modes will be selected and merged from --nevs subset according to "
				"the following probabilities (p): (default=var)\n"
				"\tplain: equiprobability (--min_prob will be ignored)\n"
				"\tvar: proportional to the mode variance, p(i)= 1/eigenvalue(i)\n"
				"\tline: lineally decreasing probability, p(i)= 1-i/nevs",false,"var","string");
		cmd.add( NM_prob );

		// Pose refinement related
		ValueArg<float> ShakeRot("","6Dref_rot", "Maximum rotational increment (degrees) for pose refinement (default=1)",false,1,"float");
		cmd.add( ShakeRot );
		ValueArg<float> ShakeMov("","6Dref_trans", "Maximum translational increment (Angstroms) for pose refinement. (default=2)",false,2,"float");
		cmd.add( ShakeMov );
		SwitchArg PoseChain("", "6Dref_chain","Independent local 6D pose refinement for each single chain.", false);
		cmd.add( PoseChain );
		ValueArg<int> PoseDelay("","6Dref_delay", "Number of iterations before first local 6D pose refinement (default=disabled)",false,-1,"int");
		cmd.add( PoseDelay );
		ValueArg<int> RefEach("","6Dref", "Number of iterations between local 6D pose refinement (default=200)",false,200,"int"); // Rotational & Traslational refinement each \"refine_each\" iters
		cmd.add( RefEach );

		// Less important parameters (string)
		ValueArg<int> Filter("","filter", "Select filtration method: 1-Fourier, 2-Kernel (By default the fastest will be selected)",false,0,"int");
		cmd.add( Filter );
		ValueArg<float> ModelThr("","cutoff2", "Density cutoff of simulated map (default = 0.001)",false,0.001,"float");
		cmd.add( ModelThr );
		//		SwitchArg NoNormEvec("", "nonormevec","Disables eigenvector \"norm=1\" normalization (default=disabled)", false);
		//		cmd.add( NoNormEvec );
		SwitchArg NoWRmsd("", "nowrmsd","Disables Gaussian Weighted RMSD (default=disabled)", false);
		cmd.add( NoWRmsd );
		SwitchArg NoTors("","notors", "Disables extra torsional potential (default=disabled)", true);
		cmd.add( NoTors );
		SwitchArg UnitMass("","unitmass", "Sets unit masses for every atom (default=disabled)", true);
		cmd.add( UnitMass );
		SwitchArg UnitNElec("","unitnelec", "Sets unit electron density for every atom (default=disabled)", true);
		cmd.add( UnitNElec );
		SwitchArg Rmsds("","morermsds", "Enables C-alpha RMSD computation (default=disabled)", true);
		cmd.add( Rmsds );
		SwitchArg PdbsOut("","morepdbs", "Saves initial (basename_model.pdb) and final (basename_fitted.pdb) models. With -F option, also the fitted CG-model (basename_fitCG.pdb). (default=disabled)", true);
		cmd.add( PdbsOut );

		//		ValueArg<std::string> Funcfile("","ss_func", "Force constants function file based on Topology and Secondary Structure.",false,"","string");
		//		cmd.add( Funcfile );
		//		ValueArg<std::string> SSfile("","ss", "Secondary Structure file with the same format as dssp2ss.pl perl script.",false,"","string");
		//		cmd.add( SSfile );
		ValueArg<std::string> PdbRef("","pdb_ref", "Reference PDB file (test only)",false,"","string");
		cmd.add( PdbRef ); // Check this if PDB-PDB morphing program...

		ValueArg<float> k2_Cutoff("","k2_c","Non-bonding distance cutoff applied to --func option (default=10A).", false, 10,"float");
		cmd.add( k2_Cutoff );
		ValueArg<float> k1_Cte("", "k1_k","Distance cutoff method stiffness constant (default=1.0)",false,1.0,"float");
		cmd.add( k1_Cte );
		ValueArg<float> k1_Cutoff("","k1_c","Distance cutoff method distance cutoff (default=10A)", false, 10,"float");
		cmd.add( k1_Cutoff );
		ValueArg<float> k0_Power("","k0_p", "Inverse Exponential's power term (default=6)",false,6.0,"float");
		cmd.add( k0_Power);
		ValueArg<float> k0_X0("", "k0_x0","Inverse Exponential's inflexion point (default=3.8A)",false,3.8,"float");
		cmd.add( k0_X0 );
		ValueArg<float> k0_Cte("", "k0_k","Inverse Exponential's stiffness constant (default=1.0)",false,1.0,"float");
		cmd.add( k0_Cte );
		ValueArg<float> k0_Cutoff("","k0_c","Inverse Exponential's distance cutoff (default=10A)", false, 10,"float");
		cmd.add( k0_Cutoff );

		ValueArg<std::string> Funcfile("","func", "ASCII file defining the force constant functions to be applied "
				"according to Topology and/or Secondary Structure. "
				"The 5 cols. format is: <SS> <n> <k> <x0> <pow>\n"
				"Where <SS> is the two character pairwise interaction identifier, <n> is the topology, and "
				"<k>,<x0>,<pow> are the corresponding inverse exponential function parameters. If --ss "
				"is not specified, only the XX pairwise interaction identifier will be considered. "
				"If <n> is \"-1\", any previously not-matched topology will be matched.",false,"","string");
		cmd.add( Funcfile );
		ValueArg<std::string> SSfile("","ss", "Secondary Structure ASCII file with 2 cols.: <index> <char>\n"
				"Where <index> is the corresponding residue index (1,2,...), and <char> is the "
				"single character SS identifier. By default SS will be computed internally (H=helix, E=strand, C=coil).",false,"","string");
		cmd.add( SSfile );

		ValueArg<float> DeltaSave("","delta_save", "RMSD increment to save a new trajectory frame (default=0.5A). "
				"If a negative integer value is introduced, a new frame will be saved each --delta_save iterations.",false,0.5,"float");
		cmd.add( DeltaSave );
		ValueArg<double> Rediag("","rediag", "RMSD ratio to trigger diagonalization. (last_diag_RMSD-current_RMSD) > rediag (default=0.1)",false,0.1,"double");
		cmd.add( Rediag );

		// More important parameters (one letter)
		ValueArg<std::string> FixSS("S","fixSS", "All dihedral coordinates belonging to residues matching the indicated SS identifieres will be fixed. "
				"Ex: \"HE\" will fix the dihedrals corresponding to alpha-helices and beta-sheets.",false,"","string");
		cmd.add( FixSS );
		ValueArg<float> FixRand2("R","fixRand2", "Randomly fixed ratio of Internal Coordinates. Example: 0.7 = 70% of IC will be randomly fixed.",false,0.5,"float");
		cmd.add( FixRand2 );
		ValueArg<float> FixRand("r","fixRand", "Randomly fixed ratio of Dihedral Coordinates. Example: 0.7 = 70% of dihedrals will be randomly fixed. (Rotational/Translational coords. always mobile)",false,0.5,"float");
		cmd.add( FixRand );
		ValueArg<std::string> FixFile("f","fixFile", "ASCII file defining the ICs to be fixed with the format:\n"
				"Protein:     \"n phi chi psi\"\n"
				"NAcid:       \"n alpha beta gamma chi epsilon zeta\"\n"
				"Inter-chain: \"n 6D\"\n"
				"Where \"n\" is the residue index (0,1,..) and the coordinate name (phi, psi, etc...) "
				"can be set to 0(fixed) or 1(mobile). "
				"Each one of the 6 inter-chain variables should be specified on separate lines in the "
				"following order: x,y,z,Rx,Ry,Rz. "
				"Note \"n\" is just the sequential residue index (starting with 0) and NOT the PDB's residue index.\n"
				"A dummy file can be generated using iMode with the --save_fixfile option.",false,"fixstring","string");
		cmd.add( FixFile );
		ValueArg<int> Contact("P","potential", "Pairwise interaction potential: (default=0)\n"
				"  0= Sigmoid function (= k/(1+(x/x0)^p), if x < c, else k=0).\n"
				"  1= Tirion's cutoff (= k, if x < c, else k=0).\n"
				"  2= Hinsen's function.\n"
				"  3= Topology & Secondary Structure (--func is mandatory).\n"
				"  4= edNMA formalism (CA-model only).\n"
				"  By default an extra torsional potential will be added.",false,0,"int");
		cmd.add( Contact );
		ValueArg<int> Iter("i","iter", "Maximum number of iterations (default=100000)",false,100000,"int");
		cmd.add( Iter );
		ValueArg<std::string> Name("o","name", "Output files basename. (default=program-name)",false,prog_name,"string");
		cmd.add( Name );
		SwitchArg Full("F", "full","Enables full-atom output models", false);
		cmd.add( Full );
		ValueArg<float> Excited("e","nex", "Excited modes range, either number [1,nevs] <integer>, or ratio [0,1) <float>. (default=0.1)",false,0.1,"int/float");
		cmd.add( Excited );
		ValueArg<float> Nevs("n","nevs", "Used modes range, either number [1,N] <integer>, or ratio [0,1) <float>. (default=0.05) (In any case, the value of --addnevs option will be added to --nevs)",false,0.05,"int/float");
		cmd.add( Nevs );
		SwitchArg Movie("t", "otraj","Outputs a Multi-PDB trajectory movie (basename_movie.pdb)", false);
		cmd.add( Movie );
		SwitchArg Chi("x","chi", "Considers first CHI dihedral angle. (default=disabled)", true);
		cmd.add( Chi );
		ValueArg<int> Model("m","model", "Coarse-Grained model: \n"
				"  0= CA, C-alpha atoms (use this model when your structure has missing side-chain atoms),\n"
				"  1= C5, 3 backbone atoms plus 2 side-chain pseudo-atoms (if present),\n"
				"  2= HA, Heavy atoms (use --keep_hydrogens option for full-atom),\n"
				"  3= NCAC(experimental) (default=2)",false,2,"int");
		cmd.add( Model );
		//		SwitchArg Morph("p", "morph","PDB-PDB morphing (introduce a PDB instead of a target <map>)", true); // , and any value at <resolution> and <cutoff>
		//		cmd.add( Morph ); // Remove it if PDB-PDB morphing program...

		// Parse the command line.
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		// =============================================================================
		strcpy(file_initial,((temp=Pdb.getValue()).c_str())); // Gets Initial PDB file name
		if(parse_verb)
			printf("Parser input: Initial PDB file: %s\n",file_initial);

		// Setting Coarse-Graining models
		// Setting model and chi
		model = Model.getValue();
		if(Chi.isSet())
			type = 2; // phi,chi,psi
		else
			type = 0; // phi,psi

		if(FixFile.isSet())
		{
			strcpy(fix_file,((temp=FixFile.getValue()).c_str())); // Gets Fix-file
			fixmodel = 3;
		}
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

		if(UnitMass.isSet())
		{
			unitmass_switch = true;
			if(parse_verb)
				printf("Parser input: Masses will be set to 1.0.\n");
		}

		if(UnitNElec.isSet())
		{
			unitnelec_switch = true;
			printf("Parser input: Electron density per atom will be set to 1.0.\n");
		}

		if( NoRandWeight.isSet() )
			randweight = false;
		if(parse_verb)
			printf("Parser input: randweight= %s\n",printbool(randweight));

		if( SCVWeight.isSet() )
			scv_weight = false;
		else
			scv_weight = true; // default enabled
		if(parse_verb)
			printf("Parser input: scv_weight = %s\n",printbool(scv_weight));

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
			contacts = 3; // Override "Contact"
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

		// Contacting method parameters
		power = k0_Power.getValue();
		cte_k0 = k0_Cte.getValue();
		x0 = k0_X0.getValue();
		cutoff_k0 = k0_Cutoff.getValue();
		cte_k1 = k1_Cte.getValue();
		cutoff_k1 = k1_Cutoff.getValue();
		cutoff_k2 = k2_Cutoff.getValue();

		if(NoTors.isSet())
		{
			notors_switch = true;
			if(parse_verb)
				printf("Parser input: Aditional Torsional Springs DISABLED!\n");
		}

		// Selecting Target input (PDB-PDB or PDB-VOL)
		//		if( Morph.isSet() ) // PDB-PDB
		//		{
		//			pdb_switch = true;
		//			strcpy(file_ref_pdb,((temp=Ref.getValue()).c_str())); // Gets Target PDB file name
		//			if(parse_verb)
		//				printf("Parser input: PDB Input Reference file: %s (PDB mode)\n",file_ref_pdb);
		//		}
		//		else // PDB-VOL
		//		{
		pdb_switch = false;
		strcpy(file_ref_vol,((temp=Ref.getValue()).c_str())); // Gets Target Map file name
		if(parse_verb)
			printf("Parser input: VOLUME Input Reference file: %s (VOLUME mode)\n",file_ref_vol);
		//		}

		filter = Filter.getValue();
		if(parse_verb)
			printf("Parser input: Filtration method: %d (Map mode)\n",filter);

		ref_thr = RefThr.getValue();
		if(parse_verb)
			printf("Parser input: Target Map threshold: %f (Map mode)\n",ref_thr);

		model_thr = ModelThr.getValue();
		if(parse_verb)
			printf("Parser input: Model Map threshold: %f (Map mode)\n",model_thr);

		resolution = Res.getValue();
		if(parse_verb)
			printf("Parser input: Resolution set to: %f\n",resolution);

		strcpy(name,((temp=Name.getValue()).c_str())); // Gets Basename
		if(parse_verb)
			printf("Parser input: Output files Base-Name: %s\n",name);

		convergence = Conver.getValue();
		if(parse_verb)
			printf("Parser input: Convergence cut-off (relative to 1st score): %e\n",convergence);

		// Set NM selection method
		strcpy(text,((temp=NM_prob.getValue()).c_str()));
		if( strcmp(text,"plain") == 0 )
		{
			if(parse_verb)
				printf("Parser input: PLAIN selection probability.\n");
			prob_method = 1;
		}
		else if( strcmp(text,"var") == 0 )
		{
			if(parse_verb)
				printf("Parser input: VAR selection probability.\n");
			prob_method = 2;
		}
		else if( strcmp(text,"line") == 0 )
		{
			prob_method = 4;
			if(parse_verb)
				printf("Parser input: LINE selection probability.\n");
		}
		else if( strcmp(text,"gaussian") == 0 )
		{
			prob_method = 5;
			if(parse_verb)
				printf("Parser input: GAUSSIAN selection probability.\n");
		}
		else if( strcmp(text,"invexp") == 0 )
		{
			prob_method = 6;
			if(parse_verb)
				printf("Parser input: INVEXP selection probability.\n");
		}
		else // Not selected
		{
			printf("\nSorry, you must specify a Normal Mode probability selection method!\nForcing Exit!\n\n");
			exit(1);
		}

		// Number of eigenvectors to be computed 1
		nevec_fact = Nevs.getValue();
		if(nevec_fact <= 0) // checking
		{
			printf("nmafit> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
			exit(1);
		}
		if(parse_verb)
			printf("Parser input: Ratio of Eigenvectors: --nevs = %8.6f\n",nevec_fact);

		nex_fact = Excited.getValue();
		if(nex_fact <= 0) // checking
		{
			printf("nmafit> Error, invalid number excited 1 eigenvectors requested (%f)!\nForcing exit!\n",nex_fact);
			exit(1);
		}
		if(parse_verb)
			printf("Parser input: Ratio of Excited Modes: --nex = %5f\n",nex_fact);

		if(Excited2.isSet())
		{
			nex_fact2 = Excited2.getValue();
			if(nex_fact2 <= 0) // checking
			{
				printf("nmafit> Error, invalid number excited 2 eigenvectors requested (%f)!\nForcing exit!\n",nex_fact2);
				exit(1);
			}
		}
		else
			nex_fact2 = Excited.getValue();
		if(parse_verb)
			printf("Parser input: Ratio of Excited Modes: --nex2 = %5f\n",nex_fact2);

		addnevs = AddNevs.getValue();
		if(parse_verb)
			printf("Parser input: Number of added modes: --addnevs = %d\n",addnevs);

		max_iter = Iter.getValue(); // Gets number of iters
		if(parse_verb)
			printf("Parser input: Maximum Number of Iterations: --iter = %d\n",max_iter);

		refine_each = RefEach.getValue(); // for both, PDB-PDB and PDB-VOL
		if(parse_verb)
			if(pdb_switch) // PDB-PDB
				printf("Parser input: PDB-PDB alignment (minRmsd) each X iters: --6Dref = %d\n",refine_each);
			else // PDB-VOL
				printf("Parser input: Rot.Tras-lational refinement each X iters: --6Dref = %d\n",refine_each);

		//		if(Seed.isSet()) // Fixed seed
		//			seed = (unsigned int) Seed.getValue(); // Gets seed
		//		else // Random seed (time initialization)
		//			seed = (unsigned int) time(0); // int32 seed = (int32)time(0); // random seed (Mersenne.h)
		//		if(parse_verb)
		//			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);
		if(Seed.isSet()) // Fixed seed
			seed = (unsigned int) Seed.getValue(); // Gets seed
		else // Random seed (time initialization)
		{
			// seed = (unsigned int) time(0); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)

			// Mon: /dev/urandom does not work in Windows...
#if defined  _WIN32 || defined _WIN64 // Windows version
			seed = (unsigned int) clock(); // int32 seed = (int32)time(0);// (Warning, in seconds!) // random seed (Mersenne.h)
#else // Linux version
			// Needed to avoid seed repetition between different runs.
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
#endif
		}

		if(parse_verb)
			printf("Parser input: Mersenne Twister's SEED: --seed = %u\n",seed);

		if(Movie.isSet())
		{
			movie_switch = true;
			if(parse_verb)
				printf("Parser input: \"--movie\" = ENABLED\n");
		}

		verbose=Verbose.getValue();
		if(parse_verb)
			printf("Parser input: \"--verb\" = %d\n",verbose);

		if(PdbRef.isSet())
		{
			bench_switch = true; // To allow RMSDs calculation when reference PDB available! (for benchmark use)
			more_rmsds = true;
			strcpy(file_ref_pdb,((temp=PdbRef.getValue()).c_str())); // Gets Reference PDB file name
			printf("Parser input: Benchmark reference PDB file: %s\n",file_ref_pdb);
		}

		conv_length = Conv_Length.getValue();
		if(parse_verb)
			printf("Parser input: Averaging length (#steps) to reach convergence --conv_length= %d\n",conv_length);

		//		if(MinProb.isSet())
		//		{
		//			min_prob = MinProb.getValue();
		//			if(min_prob<=1 && min_prob>=0)
		//				printf("Parser input: \"--min_prob\" Minimum probability value: %f\n",min_prob);
		//			else
		//			{
		//				printf("Msg(Parser): Please select a valid \"--min_prob\" value!!! [0:1]\nForcing Exit!\n");
		//				exit(1);
		//			}
		//		}

		step = Step.getValue(); // Gets step
		step_begin = step;
		step_end = step;
		if(parse_verb)
			printf("Parser input: Motion --step = %f\n",step);

		if( StepEnd.isSet() )
		{
			rand_step_switch = false;
			step_end = StepEnd.getValue(); // Gets step-end
		}
		else
			step_end = Step.getValue() / 10; // Set default step-end
		if(parse_verb)
			printf("Parser input: Some kind of annealing --step_end = %f\n",step_end);

		if( RandStep.isSet() )
			rand_step_switch = true;

		if( RandExcited.isSet() )
			rand_excited_switch = true;

		if(DeltaSave.isSet())
		{
			delta_save = DeltaSave.getValue(); // Gets Increment to save frames
			if(delta_save<0)
				deltasave_rmsd_switch=false; // non-compatible options
			printf("Parser input: Normalized Score Increment to save frames --delta_save = %f\n",delta_save);
		}

		if( NoWRmsd.isSet() )
		{
			wrmsd_switch = false;
			printf("Parser input: Weighted RMSD PDB alignment disabled\n");
		}
		if( WRmsd.isSet() )
		{
			wrmsd_switch = true;
			min_wrmsd = WRmsd.getValue();
			printf("Parser input: Weighted RMSD PDB alignment enabled, weight= %f\n",min_wrmsd);
		}
		//		if( NoNormEvec.isSet() )
		//		{
		//			norm_evec_switch = false;
		//			printf("Parser input: Eigenvector normalization disabled\n");
		//		}
		if( Rmsds.isSet() )
		{
			more_rmsds = true;
			if(parse_verb)
				printf("Parser input: C-alpha RMSDs enabled\n");
			if(!bench_switch && !pdb_switch) // if PDB-VOL && not benchmark
			{
				printf("Parser Error: Please, include a --pdb_ref file to compute RMSDs...\n");
				exit(3);
			}
		}
		if( PdbsOut.isSet() )
		{
			savemodel_switch = true;
			if(parse_verb)
				printf("Parser input: Initial/final model pdbs output enabled\n");
		}
		savefitted_switch = true; // fitted pdb should be always saved

		if(Rediag.isSet())
		{
			rediag = Rediag.getValue(); // Gets step
			if(parse_verb)
				printf("Parser input: Diagonalization threshold (score increment): --rediag = %f\n",rediag);
		}

		if(ShakeRot.isSet())
		{
			shake_rot = ShakeRot.getValue() * (M_PI/180); // radians
			if(parse_verb)
				printf("Parser input: Max. Rotational Shake: --shake_rot = %f\n",shake_rot);
		}

		if(ShakeMov.isSet())
		{
			shake_mov = ShakeMov.getValue();
			if(parse_verb)
				printf("Parser input: Max. Traslational Shake: --shake_mov = %f\n",shake_mov);
		}

		if(Excited.isSet())
			nex = Excited.getValue();

		if( PoseChain.isSet() )
			sel_chain = 0; // A different chain will be selected each 6D pose ref.
		else
			sel_chain = -1; // the whole molecule will be selected

		refine_delay = PoseDelay.getValue(); // 6D pose refinement delay (number of iterations)
		if(parse_verb)
			printf("Parser input: 6D pose refinement iter. delay: --6Dref_delay = %d\n",refine_delay);

		savemaps_switch = SaveMaps.isSet();
		fullatom_switch = Full.isSet();
		fast_switch = FastProj.isSet();
		nomodel = NoModel.isSet();
		time_switch = Timer.isSet();
		dump_switch = Dump.isSet();
		chimera_switch = Chimera.isSet();

		if(NThreads.isSet())
		{
			nthreads = NThreads.getValue();
			printf("Parser> Using %d threads.\n",nthreads);
		}

		eigensolver = Eigensolver.getValue();
		if(Eigensolver.isSet())
		{
			eigensolver_switch = true; // true --> automatic eigensolver choice will be overridden
			printf("Parser> Using eigensolver: %d (automatic choice will be overridden)\n",eigensolver);
		}

		if(InterMolec.isSet())
		{
			intermolec_factor = InterMolec.getValue();
			intermolec_switch = true;
		}

		delHeteros_switch = DelHeteros.isSet(); // Delete heteroatoms
		delHydrogens_switch = !KeepHydrogens.isSet(); // Keep hydrogens
		delWaters_switch = !KeepWaters.isSet(); // Keep waters

		//		min_nevs = MinNevs.getValue();
		//		printf("Parser input: Minimum number of eigenvectors: --min_nevs = %d\n",min_nevs);

		//		p_x0 = FuncX0.getValue();
		//		printf("Parser input: Probability method's function X0 factor: --p_x0 = %f\n",p_x0);
		//
		//		p_s = FuncSigma.getValue();
		//		printf("Parser input: Probability method's function sigma factor: --p_s = %f\n",p_s);

	} catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}

