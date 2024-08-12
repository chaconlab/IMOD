/************************************************************************
 *                           PCATOOL                                     *
 *************************************************************************
 * Program is part of the ADP package URL: http://chaconlab.org          *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2011)  *
 *************************************************************************
 *                                                                       *
 *   Principal Component Analysis tool                                   *
 *   (from old "nmatool2")
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
#include "libnma/include/libnma_io.h" // Mon's NMA library (Dihedral included)
#include "libnma/include/libnma_cg.h" // Mon's NMA Coarse-Graining library
#include "libnma/include/libnma_diag.h" // Mon's LAPACK Matrix Diagonalization calls
#include <random>

//#define FILE_NAME 500

char file_pdb1[FILE_NAME];
char file_pdb2[FILE_NAME];
char file_covar[FILE_NAME];
char dummy_string[FILE_NAME];
char name[FILE_NAME];
char text[FILE_NAME];
char saved_files[2000]; // string to buffer screen output until the program end.

// char name[FILE_NAME];
char file_ptraj[FILE_NAME];
bool weight = false; // enables mass weighting
bool cart_switch = false; // enables mass un-weighting (after PCA)
bool switch_3BB2R = false; // enables PDB translation into 3BB2R model (5-atom model)
bool format_switch = false; // enables sorting of a PDB file
bool ca_evec = false; // enables CA components extraction from 3BB2R ptraj
bool save_average= false; // save ensemble average with rmsdp in Bfactors
bool save_usedtraj= false; // save used trajectory ensemble
bool save_ref= false; // save used reference
bool save_ref2= false; // save "fine-grained" reference
bool save_rmsd= false; // save text file with "atomic" rmsds
bool save_rmsdlog= false; // save text file with rmsds respecting to the reference
bool rmsdp_switch = false; // enables rmsd profile computation
bool alignto_switch = false; // enables pdb-reference ensemble alignment
bool removeH_switch = false; // Remove hydrogens switch
bool align_switch = false; // It does alignment.
bool pca_switch = false; // Perform PCA
bool proj_switch = false; // Projects the trajectory into the Essential space
bool covar_switch = false; // Enables "just Covariance PCA"
bool save_wcart = false; // mass-weighted cartesian principal components
bool save_covar = false; // save covariance matrix
bool matrix_io_text = false; // false: binary format
bool readpdb_switch = true; // Enables pdb input reading
bool var_switch = false; // Select the range of computed PCs by variance
bool nomodel = false; // true = CG-model building and formating disabled (initial model)
bool ref_average = true; // false = not-using average as reference...
bool verb = false; // Dump some extra info.
bool remove_anchors = false; // Dump some extra info.
bool random_model = false; // Dump some extra info.

float nevec_fact = -1; // % of eigenvectors to be computed (set by parser)
float var_thr=0.0; // Variance threshold (used with var_switch)
int rmsdp_ref = -1; // Align to the "nalign"-th frame from ensemble.
int nevec = -1; // number of eigenvectors to be computed (set by parser)
int modes_saved = 20; // now controled by "nevs" (if there are enougth modes available)
int model = 2; // default CG-model (Heavy Atom, HA)
int initial = 0; // Initial frame to be accounted for
int final = -1; // Final frame to be accounted for
float change_timer = 300; // Time threshold (seconds) to use: timer, instead of Htimer

using namespace TCLAP;
void parseOptions(int argc, char** argv);
char version[]="v1.02";

// Computes the average of an ensemble (allocating memory for "avg")
// mol --> ensemble
// p_avg --> pointer to the average array (of "size" elements)
// initial --> initial molecule index
// final --> final molecule index
void ensemble_average(Macromolecule *mol, double **p_avg, int initial=0, int final=-1);

// Compute the atom-wise RMSD profile from an ensemble
// mol --> input ensemble
// avg --> input average coordiantes array
// p_rmsd --> pointer to the output rmsd profile (if NULL, then allocate memory)
void rmsd_profile(Macromolecule *mol, double *avg, double **p_rmsd, int initial=0, int final=-1);


int main( int argc, char * argv[] )
{
	printf("pcatool>\n");
	// Parsing Input
	parseOptions(argc,argv);
	printf("pcatool>\npcatool> Mon's PCA-tool %s\npcatool>\n",version);

	Macromolecule *molr1,*dummy,*ref,*refr,*aver,*mol;
	pdbIter *iter,*iterR,*iter2,*iter3;
	Molecule *Molec;
	int num_atoms1,num_atoms2,num_molecs,size,ncheck,natoms;
	M4Rot *matrix4_op;
	float matrix4[4][4];
	double *avg=NULL;
	double *rmsd=NULL;
	double *eigval=NULL;
	double *eigvect=NULL;
	double *mass=NULL;
	double *covar = NULL;
	double totvar=0.0;
	double *dev = NULL; // Deviation matrix (Q)
	float rmsdi,rmsdf,rmsd_min;
	int model_min;
	int saved_files_len; // store the "saved_files" array current length
	int accounted_mols; // number of accounted for molecules
	Tcoor pos;
	Htimer ht_timer; // timer
	timer t_timer; // timer

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
	sprintf(dummy_string,"%s.log",name);
	if( !(f_com=(FILE *)fopen(dummy_string,"w") ) )
	{
		printf("Sorry, unable to open LOG FILE: %s\n",dummy_string);
		exit(1);
	}
	sprintf(saved_files,"pcatool> Log file:                            %35s\n",dummy_string);
	for(int i=0; i<argc; i++)
		fprintf(f_com,"%s ",argv[i]);
	fclose(f_com);

	// Aligning the ensemble
	if(readpdb_switch)
	{
		// Initialize aminoacids and nucleotids
		init_aminoacids();

		// reading pdb 1
		printf( "pcatool> Reading INPUT Multi-PDB file\n" );
		molr1=new Macromolecule( file_pdb1 );
		molr1->readPDB(file_pdb1);

		// Removing Hydrogens
		if(removeH_switch)
		{
			molr1->delete_hydrogens(); // Remove Hydrogens
			printf( "pcatool> Hydrogens were removed.\n" );
		}

		// Showing input information
		if(verb)
			molr1->info(stdout);

		// Initializing Molecule iterator
		iter = (pdbIter *) new pdbIter(molr1);
		num_molecs = iter->num_molecule();
		num_atoms1 = molr1->get_num_atoms();

		// Get anchor residue indices from first loop in the Multi-PDB file
		pdbIter *iter_seg;
		Segment *seg;
		iter->pos_segment=0;
		seg = iter->get_segment();
		iter_seg = new pdbIter(seg);
		iter_seg->pos_fragment=0;
		int ri = (iter_seg->get_fragment())->getIdNumber(); // get Nt anchor residue number (PDB) <-- First residue in PDB
		iter_seg->pos_fragment=iter_seg->num_fragment()-1;
		int rf = (iter_seg->get_fragment())->getIdNumber(); // get Ct anchor residue number (PDB), i.e. the last that moves... <-- Last residue in PDB
		delete iter_seg;
		delete iter;
		printf("pcatool> Found %d molecules in %s (total atoms: %d) residues from %d to %d\n",num_molecs,file_pdb1,num_atoms1, ri, rf);





		if(final < 0)
			final = num_molecs-1; // #frames, if not specified...
		if(final >= num_molecs)
		{
			printf("Maximum number of frames exceeded, using the maximum available (%d)\n", num_molecs);
			final = num_molecs-1; // #frames, if not specified...
		}

		accounted_mols = final - initial + 1; // number of accounted for molecules
		printf("pcatool> Considering molecules from %d to %d (%d)\n",initial+1,final+1,accounted_mols);

		// Sorting/Formatting
		if(format_switch && !nomodel)
		{
			if(verb)
				molr1->info(stdout);
			printf( "pcatool> Sorting/Formatting input Multi-PDB %s\n",file_pdb1 );
			molr1->format_residues();
		}

		// Setting Coarse-Graining model for the whole trajectory
		switch(model)
		{
		case 0: // CA-only + (NH and CO)-terminal model
			printf( "pcatool> Coarse-Graining model: CA-model\n");
			// WARNING: input model should have, at least, N, CA and C atoms!
			//
			// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
			// Warning, atoms are not copied, they're just pointers to the original atoms.
			// setmass = true --> adds masses to Occupancy and Bfactor, otherwise left unchanged.
			// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
			//                   (unit mass is divided between NH or CO and their CAs)
			// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
			//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
			if(!nomodel) // Makes model (if it's needed)
				mol = cg_CA( molr1, true, false ); // masses will be computed
			else
				mol = molr1; // readed pdb will be used as it is.
			break;

		case 4: // CA-only
			printf( "pcatool> Coarse-Graining model: CA-only\n");
			// if CA selection
			mol = molr1->select( calpha2 );
			// Sets masses to a CA-only-model
			// nomass = true --> masses = 1.0
			// nomass = false --> each CA will weight its residue mass
			mass_CAonly(mol, false);
			break;

		case 3: // N,CA,C-model
			printf( "pcatool> Coarse-Graining model: N,CA,C-model\n");

			// N,CA,C selection
			mol = molr1->select( ncac2 );
			if(!nomodel) // Set masses (if it's needed)
				// mass_NCAC( mol, false );
				mass_NCAC( mol, false, true, false, false );

			break;

		case 1: // 3BB2R model
			printf( "pcatool> Coarse-Graining model: C5\n");
			mol = new Macromolecule(molr1);
			if(!nomodel) // Makes 3BB2R model (if it's needed)
			{
				// CREATES a 3BB2R reduced model
				//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
				//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
				//     There are a few exceptions: for Ala, Gly and Pro,
				//     and for the 1st and last residues.
				printf("pcatool> Creating 3BB2R reduced model:\n");
				cg_3BBR2(mol);
			}
			break;

		case 2: // Full-Atom
			printf( "pcatool> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
			mol = molr1;
			if(!nomodel) // Sets masses (if it's needed)
				mass_FA(mol, false); // sets masses
			break;
		} // Now "mol" represents the coarse-grained trajectory...

		natoms = mol->get_num_atoms(); // ensemble's number of atoms
		printf( "pcatool> Selected model number of (pseudo)atoms: %d\n", natoms );

		Condition *mobile;
		Conditions *mobile2 = new Conditions();
		mobile = new Condition(-1,-1,-1,-1,-1,-1,ri+1,rf-1,-1,-1); // get residues from Nt+1 to Ct-1, i.e. only those mobile...
		mobile2->add(mobile);
		iterR = (pdbIter *) new pdbIter(molr1); // iter to get models from "read" PDB

		if (remove_anchors) {
			Macromolecule *mold;

			mold=mol->select_cpy(mobile2);
			//        delete mol;
			mol=mold;
		}
		//



		mol->writePDB("mierda.pdb"); // ref is the average

		iter = (pdbIter *) new pdbIter(mol); // iter to get models from "coarse-grained" PDB



		if (random_model) {


			Tcoor pos2;
			Atom *at, *at2;
			std::srand ( unsigned ( std::time(0) ) );

			for(iter->pos_molecule = 0; iter->pos_molecule < num_molecs; iter->next_molecule())
			{
				// Get current molecule
				Molec = iter->get_molecule();
				iter2 = new pdbIter(Molec);


				int num_atoms = iter2->num_atom();
				int indx[num_atoms]; int indx2[num_atoms/3];

				for(iter2->pos_atom=0; !iter2->gend_atom(); iter2->next_atom() ) // iter atoms
					indx[iter2->pos_atom]=iter2->pos_atom;


				for(int f=0; f<num_atoms/3; f++ ) {
					indx2[f]=indx[f*3];
//					printf ("%d ",indx2[f]);

				}
//				printf ("\n");

				//random_shuffle(&indx[0],&indx[num_atoms]);

				random_shuffle(&indx2[1],&indx2[num_atoms/3-1]);

//				for(int f=0; f<num_atoms/3; f++ ) {
//								printf ("%d ",indx2[f]);
//
//							}
//				printf ("\n");
//				getchar();

				for(int f=0; f<num_atoms/3; f++ )
				{

					//printf ("%d atoms \n",f);
					for(int j=0; j<3; j++ ) {

					iter2->pos_atom=f*3+j;
					at = (Atom *) iter2->get_atom();
					at->getPosition(pos);
					//printf ("%d %f %f %f\n",f+j, pos[0], pos[1], pos[2]);

					iter2->pos_atom=indx2[f]+j;
					at2 = (Atom *) iter2->get_atom();
					at2->getPosition(pos2);

					//printf ("%d %f %f %f\n",indx2[f]+j, pos2[0], pos2[1],pos2[2]);
					//getchar();
					at->setPosition(pos2);
					at2->setPosition(pos);
					}
				}

				delete iter2;
			}
		}

		mol->writePDB("mierda2.pdb"); // ref is the average



		if(rmsdp_switch || align_switch)
		{
			if(rmsdp_ref >= num_molecs) // Some checking...
			{
				printf("pcatool> Please, introduce a valid reference (<%d): %d\n",num_molecs,rmsdp_ref);
				exit(2);
			}

			// Selecting the reference structure
			if(alignto_switch) // Reference provided by user
			{
				// Getting reference structure from <pdb2>
				printf( "pcatool> Reading second PDB file to be used as reference structure\n" );
				refr = new Macromolecule( file_pdb2 );
				refr->readPDB(file_pdb2);
				// Showing info
				// ref->info(stdout);

				// Removing Hydrogens
				if(removeH_switch)
				{
					refr->delete_hydrogens(); // Remove Hydrogens
					printf( "pcatool> Hydrogens were removed from reference structure.\n" );
				}

				// Setting Coarse-Graining model
				switch(model)
				{
				case 0: // CA-only + (NH and CO)-terminal model
					printf( "pcatool> Changing reference CG-model into: CA-model\n");
					// WARNING: input model should have, at least, N, CA and C atoms!
					//
					// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
					// Warning, atoms are not copied, they're just pointers to the original atoms.
					// setmass = true --> adds masses to Occupancy and Bfactor, otherwise left unchanged.
					// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
					//                   (unit mass is divided between NH or CO and their CAs)
					// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
					//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
					if(!nomodel) // Makes model (if it's needed)
						ref = cg_CA( refr, true, false ); // masses will be computed
					else
						ref = cg_CA( refr, false, false ); // input pdb masses will be used
					break;

				case 4: // CA-only
					printf( "pcatool> Changing reference CG-model into: CA-only\n");
					// if CA selection
					ref = refr->select( calpha2 );
					// Sets masses to a CA-only-model
					// nomass = true --> masses = 1.0
					// nomass = false --> each CA will weight its residue mass
					mass_CAonly(ref, false);
					break;

				case 3: // N,CA,C-model
					printf( "pcatool> Changing reference CG-model into: N,CA,C-model\n");
					// N,CA,C selection
					ref = refr->select( ncac2 );
					if(!nomodel) // Set masses (if it's needed)
						mass_NCAC( ref, false );
					break;

				case 1: // 3BB2R model
					printf( "pcatool> Changing reference CG-model into: 3BB2R\n");
					ref = new Macromolecule(refr);
					if(!nomodel) // Makes 3BB2R model (if it's needed)
					{
						// CREATES a 3BB2R reduced model
						//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
						//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
						//     There are a few exceptions: for Ala, Gly and Pro,
						//     and for the 1st and last residues.
						printf("pcatool> Creating reference 3BB2R reduced model:\n");
						cg_3BBR2(ref);
					}
					break;

				case 2: // Full-Atom
					printf( "pcatool> Changing reference CG-model into: Heavy-Atom model (no CG)\n");
					ref = refr;
					if(!nomodel) // Sets masses (if it's needed)
						mass_FA(ref, false); // sets masses
					break;
				}
			}
			else // Reference obtained from the ensemble
			{
				// COMPUTING THE FIRST ENSEMBLE AVERAGE...
				fprintf(stdout,"pcatool> Computing first ensemble average (without any alignment).\n");
				ensemble_average(mol, &avg, initial, final); // Computes the average of an ensemble

				// Getting first frame structure from ensemble to compute average
				iter->pos_molecule = initial;
				Molec = iter->get_molecule(); // get reference to the first molecule in the ensemble
				ref = new Macromolecule(); // build empty Macromolecule
				ref->add(Molec); // now "ref" has the first frame (Molec is just a reference to the molecule)
				aver = new Macromolecule(ref); // Copying the full Macromolecule (mandatory!!!).
				ref->removeAll(); // remove list memory (not contained element, i.e. the ensemble)
				num_atoms2 = aver->num_atoms(); // Get the number of atoms from the first Molecule in the ensemble
				aver->copy_coordinates(avg);

				// Compute the atom-wise RMSD profile from an ensemble
				fprintf(stdout,"pcatool> Computing first RMSD-profile from the intact ensemble.\n");
				rmsd_profile(mol, avg, &rmsd, initial, final);

				// Saving initial average PDB
				if(save_average)
				{
					aver->exchange_Pdbfact(rmsd);
					//					free(rmsd);
					sprintf(dummy_string,"%s_avgi.pdb",name);
					fprintf(stdout,"pcatool> Saving initial average structure (with atom-wise-RMSD in Bfactor) %s\n",dummy_string);
					aver->writePDB(dummy_string); // ref is the average
					printf("pcatool> Written first average: %s\n",dummy_string);
					sprintf(text,"pcatool> Initial average PDB:                 %35s\n",dummy_string);
					saved_files_len = strlen(saved_files);
					strcpy(saved_files+saved_files_len,text);
				}

				if(ref_average) // Reference from the initial average
				{
					printf("pcatool> The reference structure will be the average one.\n");
					ref = aver;
					//					ref = new Macromolecule(aver); // Copying the full Macromolecule (mandatory!!!).
				}
				else // Reference from some frame...
				{
					if(rmsdp_ref < 0) // Searching the closest to the average structure (upon alignment)
					{
						printf("pcatool> Searching the closest to the average structure (upon alignment)\n");

						// Computing RMSDs respecting to the initial average (without moving) to get the closest model.
						rmsd_min = 999999; // a very high RMSD
						for(iter->pos_molecule = initial; iter->pos_molecule <= final; iter->next_molecule())
						{
							// Get current molecule
							Molec = iter->get_molecule();
							dummy = new Macromolecule();
							dummy->add(Molec);

							// Computing RMSDs vs. average model. (without moving any of them)
							rmsdi = aver->rmsd(dummy);
							rmsdf = aver->minRmsd(dummy,matrix4);
							// Don't moving the structures at all!
							// matrix4_op = new M4Rot(matrix4);
							// dummy->applyAtoms(matrix4_op);
							// delete(matrix4_op);
							if(verb)
								printf("pcatool> RMSD form frame %3d to initial average. rmsdi= %f  rmsdf= %f\n"
										,iter->pos_molecule+1,rmsdi,rmsdf);
							if(rmsdf < rmsd_min)
							{
								rmsd_min = rmsdf;
								model_min = iter->pos_molecule; // get the closest to average model
							}
							dummy->removeAll(); // remove list memory (not contained elements, i.e. ensemble)
						}
						rmsdp_ref = model_min;
						printf("pcatool> The closest model to the initial average is %d (rmsd= %f A)\n",rmsdp_ref+1,rmsd_min);
						delete aver;
					}

					printf("pcatool> Using the %d-th structure as reference\n",rmsdp_ref+1);

					// Getting reference structure from "coarse-grained" ensemble
					iter->pos_molecule=rmsdp_ref;
					Molec = iter->get_molecule();
					ref = new Macromolecule();
					ref->add(Molec);
					// ref->writePDB("mierda.pdb");

				}




				// Saving reference either from initial average PDB or to the closest to it (up to this point nothing is moved)
				if(save_ref)
				{
					sprintf(dummy_string,"%s_ref.pdb",name);
					printf("pcatool> Saving reference PDB: %s\n",dummy_string);
					ref->writePDB(dummy_string);
					sprintf(text,"pcatool> Reference PDB:                       %35s\n",dummy_string);
					saved_files_len = strlen(saved_files);
					strcpy(saved_files+saved_files_len,text);
				}

				if(save_ref2)
				{
					// Getting reference structure from read ensemble
					iterR->pos_molecule=rmsdp_ref;
					Molec = iterR->get_molecule();
					dummy = new Macromolecule();
					dummy->add(Molec);
					// Saving "finest-grained" reference structure
					sprintf(dummy_string,"%s_ref2.pdb",name);
					printf("pcatool> Saving \"finest-grained\" reference PDB: %s\n",dummy_string);
					dummy->writePDB(dummy_string);
					sprintf(text,"pcatool> Reference \"fine-grained\" PDB:      %35s\n",dummy_string);
					saved_files_len = strlen(saved_files);
					strcpy(saved_files+saved_files_len,text);
					dummy->removeAll(); // remove list memory (not contained elements, i.e. ensemble)
				}
			} // The reference structure should have been chosen at this point
			num_atoms2 = ref->num_atoms();
			size = num_atoms2 * 3; // number of variables

			// Saving log-file with the RMSDs to the reference structure
			FILE *f_rmsdlog;
			if(save_rmsdlog)
			{
				printf("pcatool> Output alignment RMSD profile\n");
				sprintf(dummy_string,"%s_rmsd.txt",name);
				f_rmsdlog = fopen( dummy_string,"w" );
				fprintf(f_rmsdlog,"%8s %5s %9s %9s\n","#  frame","ref.","RMSDi","RMSDf");
				sprintf(text,"pcatool> Final RMSD-profile (from reference): %35s\n",dummy_string);
				saved_files_len = strlen(saved_files);
				strcpy(saved_files+saved_files_len,text);
			}



			// Aligning or saving log file with each frame RMSD respecting to the reference

			if(align_switch)
				printf("pcatool> Aligning the input trajectory to the selected structure frame %d\n",rmsdp_ref+1);
			else
				printf("pcatool> Not-Aligning input trajectory to selected structure, just computing RMSD %d\n",rmsdp_ref+1 );
			double armsd=0;
			// Alignment to the reference (closest to initial average PDB)

			double *sigrmsd;

			sigrmsd = ( double * ) malloc( num_molecs * sizeof( double ));
			int nmole=0;
			for(iter->pos_molecule = initial; iter->pos_molecule <= final; iter->next_molecule())
			{
				// Get current molecule
				Molec = iter->get_molecule();
				dummy = new Macromolecule();
				dummy->add(Molec);

				// Some checking
				num_atoms1 = dummy->get_num_atoms();
				if(num_atoms1 != num_atoms2)
				{
					printf("pcatool> I'm sorry! Number of atoms mismatch: <pdb1>= %d  <pdb2>= %d\nForcing exit!\n",num_atoms1,num_atoms2);
					exit(1);
				}



				// Computing RMSD without any alignment
				rmsdi = ref->rmsd(dummy);

				// This is needed because "ref" and "dummy" are the same atoms when dummy==ref
				// (If you do "minRmsd" with the same molecule, a bad result appears.)
				if(align_switch) { // Only applying transformation if explicitly requested!
					if(iter->pos_molecule != rmsdp_ref)
					{
						// Align respecting the selected reference PDB model. (It's fast and works nice)
						rmsdf = ref->minRmsd(dummy,matrix4);
						matrix4_op = new M4Rot(matrix4);
						dummy->applyAtoms(matrix4_op);
						delete(matrix4_op);

					}
					else
						rmsdf = 0.0;
				} else
					rmsdf = rmsdi;

				armsd+=rmsdf;
				sigrmsd[nmole]=rmsdf;
				nmole++;


				if(verb)

					printf("nmovdm> Aligned frame %3d to reference frame %3d. rmsdi= %f  rmsdf= %f\n"
							,iter->pos_molecule+1,rmsdp_ref+1,rmsdi,rmsdf);

				if(save_rmsdlog)
					fprintf(f_rmsdlog,"%8d %5d %9.6f %9.6f\n",iter->pos_molecule+1,rmsdp_ref+1,rmsdi,rmsdf);

				dummy->removeAll(); // remove list memory (not contained elements, i.e. the ensemble)
			}
			if(save_rmsdlog)
				fclose(f_rmsdlog);
			double averrmsd=armsd/num_molecs;
			double sigmarmsd=0;
			for(int i = 0; i < num_molecs; i++)
				sigmarmsd += pow(sigrmsd[i]-averrmsd,2);

			printf("pcatool> Average rmsd from ref: %.3f %.3f\n",averrmsd, sqrt(sigmarmsd/num_molecs));

			free(sigrmsd);



		}
		else // If both alignment and rmsd-profile switches are disabled, "ref" should contain a valid structure anyways
		{
			// Get initial frame (later its coordinates will be the average ones)
			printf("pcatool> Using the %d-th structure as reference\n",rmsdp_ref+1);
			if (rmsdp_ref<0) rmsdp_ref=0;
			iter->pos_molecule=rmsdp_ref;
			Molec = iter->get_molecule();
			ref = new Macromolecule();
			ref->add(Molec);
			num_atoms2 = ref->num_atoms();
			size = num_atoms2 * 3; // number of variables

			if(save_ref)
			{
				sprintf(dummy_string,"%s_ref.pdb",name);
				printf("pcatool> Saving reference PDB: %s\n",dummy_string);
				ref->writePDB(dummy_string);
				sprintf(text,"pcatool> Reference PDB:                       %35s\n",dummy_string);
				saved_files_len = strlen(saved_files);
				strcpy(saved_files+saved_files_len,text);
			}



		}
		aver = new Macromolecule(ref); // Create a full copy from "ref"
		ref->removeAll(); // remove list memory (not contained elements)

		// COMPUTING THE FINAL ENSEMBLE AVERAGE...
		fprintf(stdout,"pcatool> Computing final ensemble average.\n");
		ensemble_average(mol, &avg, initial, final); // Computes the average of an ensemble
		aver->copy_coordinates(avg);

		// Computing "atomic" RMSD-profile
		// see: ec.(4) from Yang,...Bahar. Structure (2007)
		// pdbIter *nuevo = (pdbIter *) new pdbIter(mol);
		if(readpdb_switch && (save_rmsd || save_average || pca_switch))
		{
			// Compute the atom-wise RMSD profile from an ensemble
			fprintf(stdout,"pcatool> Computing final RMSD-profile from the ensemble used in PCA.\n");
			rmsd_profile(mol, avg, &rmsd, initial, final);
		}

		// Saving average PDB
		if(readpdb_switch && save_average)
		{
			aver->exchange_Pdbfact(rmsd);
			sprintf(dummy_string,"%s_avg.pdb",name);
			fprintf(stdout,"pcatool> Saving final average structure (with atom-wise-RMSD in Bfactor) %s\n",dummy_string);
			aver->writePDB(dummy_string);
			sprintf(text,"pcatool> Final average PDB:                   %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			delete aver;
		}

		// Save Multi-PDB
		if(readpdb_switch && save_usedtraj)
		{
			sprintf(dummy_string,"%s_traj.pdb",name);
			fprintf(stdout,"pcatool> Saving the used trajectory in %s\n",dummy_string);
			FILE *f_multi;
			f_multi = fopen(dummy_string,"w"); // delete old file
			fclose(f_multi);
			for(iter->pos_molecule = initial; iter->pos_molecule <= final; iter->next_molecule())
			{
				Molec = iter->get_molecule();
				dummy = new Macromolecule();
				dummy->add(Molec);
				dummy->writeMPDB(dummy_string,iter->pos_molecule+1);
				dummy->removeAll();
			}
			sprintf(text,"pcatool> Aligned trajectory Multi-PDB:        %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
		}

		if(readpdb_switch)
		{
			// Getting masses
			if(!(mass = ( double * ) malloc( num_atoms2 * sizeof( double ) )))
			{
				printf("pcatool> Sorry, unable to allocate masses array %d bytes\nExiting!\n",(int) sizeof(double) * num_atoms2);
				exit(1);
			}
			if(weight) // If mass-weighted PCA
			{
				printf("pcatool> Getting masses from First MODEL's occupancy.\n");
				// Masses are needed for mass-weighting
				iter->pos_molecule=0; // The masses will be taken from first ensemble model.
				Molec = iter->get_molecule();
				iter2 = new pdbIter(Molec);

				// Getting & computing sqrt(masses)
				for( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
					mass[iter2->pos_atom] = sqrt((iter2->get_atom())->getPdbocc()); // get mass from occupancy
			}
			else
			{
				printf("pcatool> No mass-weighting.\n");
				// If "just coordinates" PCA, all masses set to unit.
				for(int i=0; i<num_atoms2; i++)
					mass[i] = 1.0;
			}
		}

		// ****************************************************************************************
		// * COMPUTING COVARIANCE
		// ****************************************************************************************
		if(readpdb_switch && (pca_switch || save_covar))
		{
			// Allocating/Initializing Deviation matrix (Q)
			printf("pcatool> Allocating the Deviation matrix (Q) bytes=%d\n",(int) sizeof(double)*size*accounted_mols);
			if( !(dev = (double *) malloc(sizeof(double) * size*accounted_mols) ) )
			{
				printf("pcatool> Sorry, unable to allocate Deviation matrix (Q): %d bytes\nExiting!\n",(int) sizeof(double) * size*accounted_mols);
				exit(1);
			}
			for(int i=0; i< size*accounted_mols; i++)
				dev[i] = 0.0; // initialization

			// Allocating/Initializing Covariance matrix (C)
			printf("pcatool> Allocating the Covariance matrix (C) bytes=%d\n",(int) sizeof(double) * size*(size+1)/2);
			if( !(covar = (double *) malloc(sizeof(double) * size*(size+1)/2)) )
			{
				printf("pcatool> Sorry, unable to allocate Covariance matrix (C): %d bytes\nExiting!\n",(int) sizeof(double) * size*(size+1)/2);
				exit(1);
			}
			for(int i=0; i< size*(size+1)/2; i++)
				covar[i] = 0.0; // initialization

			// Computing the transpose of the Deviation matrix (Q^t)
			ht_timer.restart(); // timer
			t_timer.restart(); // timer
			printf("pcatool> Computing the Deviation matrix (Q) ");
			fflush(stdout);
			int offset=0;
			int nmol = 0;
			for(iter->pos_molecule = initial; iter->pos_molecule <= final; iter->next_molecule())
			{
				// Get current molecule
				Molec = iter->get_molecule();
				iter2 = new pdbIter(Molec);
				offset = size * nmol; // molecule start

				for(iter2->pos_atom=0; !iter2->gend_atom(); iter2->next_atom() ) // iter atoms
				{
					(iter2->get_atom())->getPosition(pos);
					dev[offset + iter2->pos_atom*3]     = (pos[0]-avg[iter2->pos_atom*3])*mass[iter2->pos_atom];
					dev[offset + iter2->pos_atom*3 + 1] = (pos[1]-avg[iter2->pos_atom*3 + 1])*mass[iter2->pos_atom];
					dev[offset + iter2->pos_atom*3 + 2] = (pos[2]-avg[iter2->pos_atom*3 + 2])*mass[iter2->pos_atom];
				}

				nmol++; // counting "in range" molecules...
				delete iter2;
			}
			if(t_timer.elapsed() > change_timer)
				printf("%s\n",t_timer.print_time());
			else
				printf("%s\n",ht_timer.print_time_sec());

			// Computing Covariance matrix (C)
			ht_timer.restart(); // timer
			t_timer.restart(); // timer
			printf("pcatool> Computing the Covariance matrix (C) ");
			fflush(stdout);
			for(int i=0; i<size; i++) // iter rows. (coordinate)
				for(int j=i; j<size; j++) // iter cols. (coordinates, including diagonal)
				{
					offset = i + j*(j+1)/2;
					// "i" and "j" specify the coordinate offset
					for(int k=0; k<accounted_mols*size; k+=size) // iter molecules
						covar[ offset ] += dev[k+i] * dev[k+j]; // coord. "i" vs. "j" (for every molecule)
					covar[ offset ] /= accounted_mols;
				}
			if(t_timer.elapsed() > change_timer)
				printf("%s\n",t_timer.print_time());
			else
				printf("%s\n",ht_timer.print_time_sec());

			// Free Deviation matrix (Q)
			if(!proj_switch) // If "dev" matrix not needed anymore (It may be required for projection into essential space)
			{
				free(dev);
				printf("pcatool> Freed Deviation matrix (Q)\n");
			}

			// Saving a triangular packed matrix
			if(save_covar)
			{
				ht_timer.restart(); // timer
				t_timer.restart(); // timer
				if(matrix_io_text) // Text format
				{
					sprintf(dummy_string,"%s_covar.txt",name);
					save_matrix(covar, size, dummy_string);
					printf("pcatool> Saving covariance matrix (C=%d^2) (text): %s ",size,dummy_string);
					sprintf(text,"pcatool> Covariance matrix (text):            %35s\n",dummy_string);
				}
				else // Binary format
				{
					sprintf(dummy_string,"%s_covar.bin",name); // ptraj name
					printf("pcatool> Saving covariance matrix (C=%d^2) (binary): %s ",size,dummy_string);
					save_matrixB(covar, size, dummy_string);
					sprintf(text,"pcatool> Covariance matrix (binary):          %35s\n",dummy_string);
				}
				saved_files_len = strlen(saved_files);
				strcpy(saved_files+saved_files_len,text);
				if(t_timer.elapsed() > change_timer)
					printf("%s\n",t_timer.print_time());
				else
					printf("%s\n",ht_timer.print_time_sec());
			}
		} // END COVAR COMPUTATION
	}
	else
	{
		printf( "nmadm> Skipped input PDB reading!\n" );
	}

	// Reading input covariance matrix (C)
	if( covar_switch )
	{
		ht_timer.restart(); // timer
		t_timer.restart(); // timer
		if(matrix_io_text) // Text format
		{
			// Reading a triangular packed matrix (into memory)
			// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
			printf("pcatool> Reading covariance (C) from text file: %s ",file_covar);
			read_matrix(&covar,&size,file_covar);
		}
		else
		{
			// Reading a triangular packed matrix (into memory) from BINARY file
			// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
			printf("pcatool> Reading covariance (C) from binary file: %s ",file_covar);
			read_matrixB(&covar,&size,file_covar);
		}
		if(t_timer.elapsed() > change_timer)
			printf("%s\n",t_timer.print_time());
		else
			printf("%s\n",ht_timer.print_time_sec());
	}
	else // Deleting all models (they are not needed anymore!)
	{
		// Setting Coarse-Graining model
		switch(model)
		{
		case 0: // CA-only + (NH and CO)-terminal model
			mol->removeAll();
			delete molr1;
			break;

		case 4: // CA-only
			mol->removeAll();
			delete molr1;
			break;

		case 3: // N,CA,C-model
			mol->removeAll();
			delete molr1;
			break;

		case 1: // 3BB2R model
			delete molr1;
			delete mol;
			break;

		case 2: // Full-Atom
			delete molr1;
			break;
		}
		printf("pcatool> Freed input data.\n");
	}
	printf("pcatool> Covariance matrix size= %d x %d\n",size,size);

	// Compute the total variance directly from the covariance matrix
	if(pca_switch || save_covar || covar_switch)
	{
		for(int i=0; i<size; i++)
			totvar += covar[i + i*(i+1)/2]; // diagonal in packed storage
		printf("pcatool> Total variance computed from covariance matrix, tr(C)= %10f\n",totvar);
	}

	if(	pca_switch || covar_switch )
	{
		// Selecting principal components range to be computed
		if(var_switch)
		{
			// Buffering covariance matrix, because DSPEVX destroys it!
			double *covar2;
			printf("pcatool> Allocating Covariance matrix Buffer (C2) bytes=%d\n",(int) sizeof(double) * size*(size+1)/2);
			if( !(covar2 = (double *) malloc(sizeof(double) * size*(size+1)/2)) )
			{
				// IF NOT ENOUGHT MEMORY, WE SHOULD TRY TO SAVE A FILE!
				printf("pcatool> Sorry, unable to allocate Covariance matrix Buffer (C2): %d bytes\nExiting!\n",(int) sizeof(double) * size*(size+1)/2);
				exit(1);
			}
			for(int i=0; i< size*(size+1)/2; i++)
				covar2[i] = covar[i]; // buffering

			ht_timer.restart(); // timer
			t_timer.restart(); // timer
			printf("pcatool> Computing ALL (%d) eigenvalues with DSPEVX: ",size);
			fflush(stdout);

			// Compute all eigenvalues (no eigenvectors)
			diag_dspevx(covar2,&eigval, &eigvect, size,1,size,0.0,0.0,'N','A',NULL);
			if(t_timer.elapsed() > change_timer)
				printf("%s\n",t_timer.print_time());
			else
				printf("%s\n",ht_timer.print_time_sec());

			// Free covariance matrix buffer (C2)
			free(covar2);

			printf("pcatool> Reversing eigenvalues ordering.\n");
			reverse_ptraj(eigval,NULL,size,size);

			float dummy=0.0;
			int cont=0;
			do {
				dummy += (float)eigval[cont]/totvar;
				if( dummy > var_thr )
					nevec = cont;
				cont++;
			} while(cont<size && nevec<0);
			if(nevec<0)
				nevec=size;
			else
				nevec++;

			printf("pcatool> The first %d principal components explain the %f%% of variance.\n",nevec,100*dummy);




		}
		else
		{
			// Number of eigenvectors to be computed (we need to know "size" first)
			if(nevec_fact >= 1.0) // number of modes
			{
				nevec = (int) nevec_fact;
				printf( "pcatool> Range of computed principal components: 1-%d\n", nevec);
			}
			else
			{
				if (nevec_fact>0.95) nevec = size;
				else nevec = (int) (nevec_fact * size);
				printf( "pcatool> Range of computed principal components: 1-%d (%.0f%)\n", nevec, nevec_fact*100);
			}
			// Checking
			if(nevec > size)
			{
				printf("pcatool> Sorry, more eigenvectors requested (%d) than available (%d), forcing maximum.\n",nevec,size);
				nevec = size;
			}
			else if(nevec <= 0) // checking
			{
				printf("pcatool> Error, invalid number of eigenvectors requested %d (%f)!\nForcing exit!\n",nevec,nevec_fact);
				exit(1);
			}
		}

		// ****************************************************************************************
		// * PERFORMING PCA
		// ****************************************************************************************
		// You only save what you compute!
		modes_saved = nevec;
		if(verb > 1)
			printf( "pcatool> Number of principal components to be saved: %d\n", modes_saved);

		//*  DSPEVX computes SELECTED eigenvalues and, optionally, eigenvectors
		//*  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
		//*  can be selected by specifying either a range of values or a range of
		//*  indices for the desired eigenvalues.
		//   (SMALL MEMORY REQUIREMENTS) (default compute eigenvalues and eigenvectors)
		ht_timer.restart(); // timer
		t_timer.restart(); // timer
		printf("pcatool> Computing %d eigenvectors/values of %d using DSPEVX: ",nevec,size);
		fflush(stdout);
		diag_dspevx(covar,&eigval, &eigvect, size, size-nevec+1, size, 'V');
		if(t_timer.elapsed() > change_timer)
			printf("%s\n",t_timer.print_time());
		else
			printf("%s\n",ht_timer.print_time_sec());
		free(covar);

		printf("pcatool> Reversing eigenpairs ordering.\n");
		reverse_ptraj(eigval,eigvect,size,modes_saved);

		int max=modes_saved;

		if(max > nevec)
			max = nevec;
		printf("pcatool> Showing the %d first eigenvalues:\n",max);
		printf("pcatool>\npcatool>   %4s %10s %10s\n","MODE","EIGENVALUE","%VARIANCE");


		float dummy=0.0;
		for(int i=0; i<max; i++) {
			dummy += (float)eigval[i]/totvar;
			printf("pcatool>   %4d  %14e  %8.6f\n",i+1,eigval[i],dummy);
		}
		printf("pcatool>\n");

		if(weight)
		{
			sprintf(dummy_string,"%s_wpca.evec",name); // ptraj name
			printf("pcatool> Saving %d principal components in \"ptraj\" format (mass-weighted): %s\n",modes_saved,dummy_string);
			sprintf(text,"pcatool> MW Principal Components (ptraj):     %35s\n",dummy_string);
		}
		else
		{
			sprintf(dummy_string,"%s.evec",name); // ptraj name
			printf("pcatool> Saving %d principal components in \"ptraj\" format: %s\n",modes_saved,dummy_string);
			sprintf(text,"pcatool> Principal Components (ptraj):        %35s\n",dummy_string);
		}
		saved_files_len = strlen(saved_files);
		strcpy(saved_files+saved_files_len,text);
		save_ptraj_modes(dummy_string, size, 0, modes_saved, eigval, eigvect, false);

		if(cart_switch)
		{
			// Mass Weight/Un-weight eigenvectors.
			weight_ptraj(eigvect,mass,size,modes_saved,false);

			sprintf(dummy_string,"%s_cpca.evec",name); // ptraj name
			printf("pcatool> Saving %d principal components in \"ptraj\" format (un-weighted cartesian): %s\n",modes_saved,dummy_string);
			sprintf(text,"pcatool> Un-MW Principal Components (ptraj):  %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);
			save_ptraj_modes(dummy_string, size, 0, modes_saved, eigval, eigvect, false);
		}

		// Projecting the computed essential space into the deviations (this is equivalent to pcasuite's trajectory compression)
		if(proj_switch)
		{
			// Opening trajectory projection file (it will be dumped directly into a file on-the-fly)
			FILE *f_proj;
			sprintf(dummy_string,"%s.proj",name);
			f_proj = fopen( dummy_string,"w" );

			// Append file output summary
			sprintf(text,"pcatool> Trajectory projection file:          %35s\n",dummy_string);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files+saved_files_len,text);

			// Writing header
			fprintf(f_proj,"%12s","# Frame, PC:");
			for(int i=1; i<=modes_saved; i++)
				fprintf(f_proj," %11d",i);
			fprintf(f_proj,"\n");

			// Projecting deviation matrix (Q) into eigenvectors
			ht_timer.restart(); // timer
			t_timer.restart(); // timer
			printf("pcatool> Projecting deviation matrix (Q) into the computed essential space: ");
			fflush(stdout);
			for(int k=0; k<accounted_mols*size; k+=size) // iter molecules
			{
				fprintf(f_proj,"%7d     ",initial + (k/size)+1);
				for(int n=0; n<modes_saved; n++)
				{
					double proj = 0.0; // frame deviation vs. current principal component projection
					for(int i=0; i<size; i++) // iter coordinates
						proj += dev[k+i] * eigvect[n*size + i]; // coord. "i" for every molecule "k"
					fprintf(f_proj," %11f",proj);
				}
				fprintf(f_proj,"\n");
			}
			if(t_timer.elapsed() > change_timer)
				printf("%s\n",t_timer.print_time());
			else
				printf("%s\n",ht_timer.print_time_sec());

			// Free Deviation matrix (Q)
			free(dev);
			printf("pcatool> Freed Deviation matrix (Q)\n");
		}

	} // END COVAR PCA

	//	delete molr1;
	if(readpdb_switch)
		delete iter;
	printf("pcatool>\npcatool> SAVED FILES:\n%s",saved_files);
	printf("pcatool>\npcatool> Bye!\n");
	return 0;
}

/*==============================================================================================*/
void parseOptions(int argc, char** argv)
{
	std::string temp;
	CmdLine cmd("pcatool","PCA-tool", version );

	try {
		// Define required arguments not labeled
		UnlabeledValueArg<std::string> Pdb1("pdb","Input Multi-PDB file","default","pdb");
		cmd.add( Pdb1 );


		// Define labeled arguments
		ValueArg<int> Verb("","verb", "Verbose level (default=0, no verbose)",false,0,"int");
		cmd.add( Verb );

		ValueArg<std::string> Name("o","name", "Sets the output files base-name.",false,"pca","string");
		cmd.add( Name );

		SwitchArg SaveCovar("","save_covar", "Saves covariance matrix in binary packed storage format. (default=disabled)", true);
		cmd.add( SaveCovar );

		SwitchArg SaveWCart("", "save_wcart","Save Mass-weighted Cartesian pricipal components as <basename_wpca.evec> (default=disabled)", false);
		cmd.add( SaveWCart );

		SwitchArg SaveAverage("", "save_avg","Used together with --rmsd_profile, saves average PDB", false);
		cmd.add( SaveAverage );

		SwitchArg SaveReFinest("", "save_refinest","Saves the \"finest-grained\" reference PDB", false);
		cmd.add( SaveReFinest );

		SwitchArg SaveRef("", "save_ref","Saves the reference PDB", false);
		cmd.add( SaveRef );

		SwitchArg SaveUsedtraj("", "save_usedtraj","Saves the aligned ensemble (Multi-PDB)", false);
		cmd.add( SaveUsedtraj );

		SwitchArg NoModel("", "nomodel","Disables PDB model building. Warning1: introduced PDB model must match the selected CG (-m option), and the masses should be already set in Occupancy PDB-field. (default=disabled)", false);
		cmd.add( NoModel );

		SwitchArg Cartesian("", "cartesian","Enables mass un-weighting after PCA (default=false)", false);
		cmd.add( Cartesian );

		SwitchArg Weight("", "weight","Enables Mass-weighting cartesian displacement vector (default=false)", false);
		cmd.add( Weight );

		SwitchArg Format("", "format","Sort properly the atoms acording to any model  (<pdb2> will be the output) (you should introduce something into <ptraj>, it will be ignored)", false);
		cmd.add( Format );

		SwitchArg RemoveH("", "noH","Removing Hydrogen atoms from output files and computations.", false);
		cmd.add( RemoveH );


		SwitchArg RemoveA("", "remove_anchors","Removing anchors loop", false);
		cmd.add( RemoveA );

		SwitchArg Random_model("", "random_model","Random model", false);
		cmd.add( Random_model );


		SwitchArg RMSDlog("", "rmsd_log","Save text file with the RMSD of each frame respecting to the reference structure.", false);
		cmd.add( RMSDlog );

		SwitchArg RMSDprofile("", "rmsd_profile","Compute the atomic RMSD profile. Plain text file output <basename_rmsd.txt>, and average PDB with RMSD in Bfactors <basename_avg.pdb>.", false);
		cmd.add( RMSDprofile );

		SwitchArg Proj("p","proj", "Projects the trajectory into the Essential space and generates a text file <basename.proj> with the projections of one frame per row. (default=disabled)", true);
		cmd.add( Proj );

		ValueArg<int> NAlign("","nalign", "Number of frame which will be used as alignment reference (from 1 to #models). "
				"If integer, using n-th frame as reference (1,2,..,N). "
				"If <1 the ensemble will be aligned to the closest frame respecting to the average. (default=disabled)",false,0,"int");
		cmd.add( NAlign );

		ValueArg<std::string> Alignto("","alignto", "PDB which will be used as reference during ensemble alignment."
				"If integer, using n-th frame as reference (1,2,..,N). "
				"If [rmsd_profile]<1 the ensemble will be aligned to the closest frame respecting to the average. (default=0)"
				,false,"nmatool","string");
		cmd.add( Alignto );

		ValueArg<std::string> Covar("","covar", "Reads a covariance matrix from file and performs PCA",false,"pcatool","string");
		cmd.add( Covar );

		ValueArg<float> RVar("r","rvar", "Select the number of principal components by Ratio of variance explained. (example: 0.9 = 90%% of var.) (default=disabled)",false,0.9,"int/float");
		cmd.add( RVar );

		ValueArg<float> Nevs("n","nevs", "Used principal components range, either number [1,N] <integer>, or ratio [0,1) <float>. (default=10)",false,10,"int/float");
		cmd.add( Nevs );

		SwitchArg Align("", "align","Align the each input trajectory frame to the selected reference.", false);
		cmd.add( Align );

		ValueArg<int> FrameF("f","final", "Final frame (default= #frames, last)",false,-1,"int");
		cmd.add( FrameF );
		ValueArg<int> FrameI("i","initial", "Initial frame (default= 1, first)",false,1,"int");
		cmd.add( FrameI );

		SwitchArg PCA("", "pca","Enable PCA computation from the input trajectory.", false);
		cmd.add( PCA );

		ValueArg<int> Model("m","model", "Input Coarse-Graining model: 0=CA-IC, 1=3BB2R, 2=Full-Atom, 3=NCAC, 4=CA-only. Note input model should match input data! (default=2)",false,2,"int");
		cmd.add( Model );

		// Parse the command line. -------------------------------------------------------
		cmd.parse(argc,argv);

		// Getting the command line arguments.
		strcpy(file_pdb1,((temp=Pdb1.getValue()).c_str())); // Gets PDB file name
		strcpy(name,((temp=Name.getValue()).c_str())); // Gets Base-Name

		if(PCA.isSet())
		{
			pca_switch = true;
			printf("pcatool> Parser: Perform Principal Component Analysis.\n");
		}

		if(NAlign.isSet())
		{
			//	align_switch = true;
			save_ref = true;
			alignto_switch = false;
			ref_average = false;
			rmsdp_ref = NAlign.getValue();
			if (rmsdp_ref==-1) {
				rmsdp_switch = true;
				printf("pcatool> Select reference close to average.\n");
			}
			else {
				rmsdp_switch = true;
				rmsdp_ref--;
				printf("pcatool> Select reference %d-th frame.\n",rmsdp_ref);
			}
		}

		if(Alignto.isSet())
		{
			strcpy(file_pdb2,((temp=Alignto.getValue()).c_str())); // now the file name is in <pdb2>
			align_switch = true;
			alignto_switch = true;
			ref_average = false;
			printf("pcatool> Parser: Aligning to %s\n",file_pdb2);
		}

		if(Align.isSet())
		{
			align_switch = true;
			printf("pcatool> Parser: Align the Multi-PDB trajectory into selected reference\n");
		}

		if(RMSDprofile.isSet())
		{
			rmsdp_switch = true;
			printf("pcatool> Parser: Compute the RMSD profile and save it into .pdb and .txt\n");
		}

		if(RMSDlog.isSet())
		{
			save_rmsdlog = true;
			printf("pcatool> Parser: Save a text file with the RMSDs respecting to the reference.\n");
		}

		if(RemoveH.isSet())
		{
			removeH_switch = true;
			printf("pcatool> Parser: Remove hydrogens\n");
		}

		if(RemoveA.isSet()) {
			remove_anchors = true;
			save_ref2 = true;
			printf("pcatool> Parser: Remove anchors\n");
		}
		if(Random_model.isSet()) {
			random_model = true;
			printf("pcatool> Parser: Random Model\n");
		}

		if(Proj.isSet())
		{
			proj_switch = true;
			printf("pcatool> Parser: Projects the trajectory into the Essential space.\n");
		}

		if(SaveUsedtraj.isSet())
		{
			save_usedtraj = true;
			printf("pcatool> Parser: Saving the used trajectory\n");
		}

		if(SaveRef.isSet())
		{
			save_ref = true;
			printf("pcatool> Parser: Saving reference PDB\n");
		}

		if(SaveReFinest.isSet())
		{
			save_ref2 = true;
			printf("pcatool> Parser: Saving \"finest-grained\" reference PDB\n");
		}

		if(SaveAverage.isSet())
		{
			save_average = true;
			printf("pcatool> Parser: Saving average PDB/s.\n");
		}

		if(Weight.isSet())
		{
			weight = true;
			printf("pcatool> Parser: Using mass-weighted coordinates ( multiply by sqrt(mi), mi= atomic mass )\n");
		}

		if(Cartesian.isSet())
		{
			cart_switch = true;
			printf("pcatool> Parser: In the end, it un-weights the coordinates ( divide by sqrt(mi) )\n");
		}

		if(Format.isSet())
		{
			format_switch = true;
			printf("pcatool> Parser: Formating PDB properly.\n");
		}

		if(SaveWCart.isSet())
		{
			save_wcart = true;
			printf("pcatool> Parser: Save mass-weighted cartesian principal components.\n");
		}

		if(SaveCovar.isSet())
		{
			save_covar = true;
			printf("pcatool> Parser: Saving covariance matrix (C) in binary packed storage format.\n");
		}

		// Number of eigenvectors to be computed
		nevec_fact = Nevs.getValue();
		if(nevec_fact <= 0) // checking
		{
			printf("pcatool> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
			exit(1);
		}

		// Range of eigenvectors to be computed
		var_thr = RVar.getValue();
		var_switch = RVar.isSet();
		if(var_thr <= 0 || var_thr > 1.0) // checking
		{
			printf("pcatool> Error, invalid Ratio of variance requested (%f)! It should be in range (0,1]\nForcing exit!\n",var_thr);
			exit(1);
		}

		if(FrameI.isSet())
		{
			initial = FrameI.getValue();
			printf("pcatool> Parser: Initial frame= %d\n",initial);
			initial--;
		}
		if(FrameF.isSet())
		{
			final = FrameF.getValue();
			printf("pcatool> Parser: Final frame= %d\n",final);
			final--;
		}

		if(Verb.isSet())
		{
			verb = Verb.getValue();
			printf("pcatool> Parser: Verbose level= %d\n",verb);
		}

		if(Covar.isSet())
		{
			covar_switch = true;
			weight = false;
			strcpy(file_covar,((temp=Covar.getValue()).c_str())); // Gets Covariance file name
			printf("pcatool> Parser: Using a pre-computed binary covariance matrix: %s\n",file_covar);
			printf("pcatool>         Warning, all previous parser options will be overridden!\n");
			pca_switch = false;
			readpdb_switch = false;
		}

		// Setting model
		model = Model.getValue();

		if(NoModel.isSet())
			nomodel = true;

	} catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}

// Computes the average of an ensemble (allocating memory for "avg")
// mol --> ensemble
// p_avg --> pointer to the average array (of "size" elements)
// initial --> initial molecule index
// final --> final molecule index
void ensemble_average(Macromolecule *mol, double **p_avg, int initial, int final)
{
	double *avg = *p_avg;
	Molecule *Molec;
	pdbIter *iter2;
	Tcoor pos;
	int num_atoms,num_molecs,ncheck,size,nmols;

	pdbIter *iter = (pdbIter *) new pdbIter(mol); // this must be done here! (mandatory)

	num_molecs = iter->num_molecule();
	if(final<0)
		final = num_molecs-1;

	// Computing initial average
	nmols = 0; // counting molecules
	for(iter->pos_molecule = initial; iter->pos_molecule <= final; iter->next_molecule())
	{
		// Get current molecule
		Molec = iter->get_molecule();
		iter2 = new pdbIter(Molec);

		// Some initializations
		if(iter->pos_molecule == initial)
		{
			num_atoms = iter2->num_atom();
			size = num_atoms * 3; // setting the number of Cartesian coordinates from the first molecule
			if(avg==NULL)
			{
				avg = (double *) malloc( sizeof(double) * size ); // allocating memory (if necessary)
				*p_avg = avg; // outputs the allocated address!
			}
			for(int i=0; i<size; i++)
				avg[i] = 0.0; // initialization
		}

		// Screening Molec's atoms...
		ncheck=0; // to check current molecule number of atoms
		for(iter2->pos_atom=0; !iter2->gend_atom(); iter2->next_atom() )
		{
			(iter2->get_atom())->getPosition(pos);
			avg[ncheck*3]     += pos[0]; // accumulate
			avg[ncheck*3 + 1] += pos[1];
			avg[ncheck*3 + 2] += pos[2];
			ncheck++;
		}
		if(ncheck != num_atoms) // Checking whether the correct number of atoms exist in the full trajectory
		{
			printf("ensemble_average> I'm sorry, number of atoms mismatch found!  first= %d  current(frame=%d)= %d\n",num_atoms,iter->pos_molecule+1,ncheck);
			exit(1);
		}
		delete iter2;
		nmols++; // counting molecules
	}

	for(int i=0; i<num_atoms; i++)
	{
		avg[i*3] /= nmols;
		avg[i*3 + 1] /= nmols;
		avg[i*3 + 2] /= nmols;
	}

	delete iter;
}

// Compute the atom-wise RMSD profile from an ensemble
// mol --> input ensemble
// avg --> input average coordiantes array
// p_rmsd --> pointer to the output rmsd profile (if NULL, then allocate memory)
void rmsd_profile(Macromolecule *mol, double *avg, double **p_rmsd, int initial, int final)
{
	Molecule *Molec;
	pdbIter *iter2;
	Tcoor pos;
	double *rmsd = *p_rmsd;
	int num_molecs,num_atoms,nmols;

	pdbIter *iter = (pdbIter *) new pdbIter(mol); // this must be done here! (mandatory)

	num_molecs = iter->num_molecule();
	if(final<0)
		final = num_molecs-1;

	nmols=0;
	for(iter->pos_molecule = initial; iter->pos_molecule <= final; iter->next_molecule())
	{
		// Get current molecule
		Molec = iter->get_molecule();
		iter2 = new pdbIter(Molec);

		if(iter->pos_molecule == initial)
		{
			num_atoms = iter2->num_atom();
			if(rmsd==NULL)
			{
				rmsd = (double *) malloc( sizeof(double) * num_atoms);
				*p_rmsd = rmsd; // outputs the allocated address!
			}
			for(int i=0; i<num_atoms; i++)
				rmsd[i] = 0.0;
		}

		for(iter2->pos_atom=0; !iter2->gend_atom(); iter2->next_atom() )
		{
			(iter2->get_atom())->getPosition(pos);
			rmsd[iter2->pos_atom] += pow(pos[0]-avg[iter2->pos_atom*3],2) +
					pow(pos[1]-avg[iter2->pos_atom*3 + 1],2) +
					pow(pos[2]-avg[iter2->pos_atom*3 + 2],2); // accumulate
			//						if(iter->pos_molecule == initial)
			//						{
			//							printf("pos=(%f,%f,%f) avg=(%f,%f,%f)",pos[0],pos[1],pos[2],avg[iter2->pos_atom*3],avg[iter2->pos_atom*3 + 1],avg[iter2->pos_atom*3 + 2]); // rmsd[iter2->pos_atom]);
			//						}
		}
		delete iter2;
		nmols++; // counting molecules
	}

	for(int i=0; i<num_atoms; i++)
	{
		rmsd[i] = sqrt( rmsd[i]/nmols );
		//					fprintf(stderr,"rmsd[%d]= %lf (accounted_mols= %d)\n",i,rmsd[i],accounted_mols);
	}

	delete iter;
}
