/************************************************************************
*                           NMACM                                       *
*************************************************************************
* Program is part of the ADP package URL: http://sbg.cib.csic.es        *
* (c) Jose Ramon Lopez-Blanco, Jose Ignacio Garzon and Pablo Chacon     *
* Structural Bioinformatics Group (2004-2009)                           *
*************************************************************************
*                                                                       *
*   Tool for NMA computation in Cartesian Coordinate Space (CCS)        *
*   a) Coarse-Graining models: CA, 3BB2R(5-atom), Full-Atom             *
*   b) Contacting rules: cutoff, inverse-exponential, K-list            *
*   c) Mass-weighted CCS modes computation                              *
*                                                                       *
*************************************************************************
* This program is free software; you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation; either version 2 of the License, or     *
* (at your option) any later version.                                   *
************************************************************************/

#include <stdio.h> // Needed to define EOF
#include "cmdl/CmdLine.h"
//#include "libnma/nma.h"
#include "libnma/include/libnma_io.h" // Mon's NMA Input-Output library
#include "libnma/include/libnma_cg.h" // Mon's NMA Coarse-Graining library
#include "libnma/include/libnma_hessian.h" // Mon's NMA Hessian Matrix library
#include "libnma/include/libnma_diag.h" // Mon's LAPACK Matrix Diagonalization calls
#include "libnma/include/libnma_misc.h" // Mon's NMA related library
#include "libnma/include/libnma_def.h" // Pablo's Deformability library

/*==============================================================================================*/
char version[]="v1.05";
// Input variables
char file_pdb[FILE_NAME]; // pdb input
char Kfile[FILE_NAME]; // Kfile input
char name[FILE_NAME];
char text[FILE_NAME];
char dummy_string[FILE_NAME];

int contacts;
double cutoff;
bool no_mass=false;
bool null=false;
float nevec_fact=0;
float cte=1;
//double power=6;
//double x0=3.8;
bool norm_modes=true;
int modes_saved=56;
int model=0; // coarse graining model
int nevec=0;
int verb = 0; // Enables different levels of dirty debugging output!
bool save_cart=true; // false --> cartesian modes
bool Kout_switch=false; // outputs contacts
bool debug_3BB2R = false; // Allows 3BB2R PDB input (initial model)
bool saveformat_switch = false; // = true --> save formated input PDB
bool savemodel_switch = true; // = true --> save 3BB2R-model input PDB
//bool nevecfact_switch = false; // = true --> compute a fraction of available eigenvectors
bool deform_switch=false; // = true --> enable deformability computations

// INPUT PARAMETERS (CUSTOMIZABLE by parser)
float cte_k0 = 1.0;
float cte_k1 = 1.0;
float x0 = 3.8;
float power = 6.0;
float cutoff_k0 = 10;
float cutoff_k1 = 10;
bool save_wcart=false; // true --> mass-weighted cartesian modes

/*==============================================================================================*/
using namespace TCLAP;
void parseOptions(int argc, char** argv);

/*==============================================================================================*/
int main( int argc, char * argv[] )
{
	int i, j, cont, index;
	Htimer ht_timer; // timer

	// Parsing input
	parseOptions(argc,argv);

	// Initialize aminoacids and nucleotids
	init_aminoacids();

	// CA  Conditions
	Condition *calpha = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	calpha->add( " CA  " );
	Conditions *calpha2 = new Conditions();
	calpha2->add( calpha );

	// NCAC  Conditions ( N-, CA-, C- selection)
	Condition *ncac = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	ncac->add( " N   " );
	ncac->add( " CA  " );
	ncac->add( " C   " );
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
	for(i=0; i<argc; i++)
		fprintf(f_com,"%s ",argv[i]);
	fprintf(f_com,"\n");

	// Reading pdb
	int num_res=0;
	int num_atoms=0;
	Macromolecule *molr = new Macromolecule( "pdb" );
	molr->readPDB( file_pdb );

	// Formating PDB first
	printf( "nmacm>  Formatting residues names\n" );
	molr->format_residues(); // formating allways
	if(saveformat_switch)
	{
		sprintf(dummy_string,"%s_format.pdb",name);
		molr->writePDB( dummy_string );
	}

    molr->info(stdout);

    // Some checking...
	num_atoms = molr->get_num_atoms();
	printf( "nmacm>\nnmacm> Atoms read -> %d\n", num_atoms );
	if(num_atoms==0)
	{
		printf( "nmacm>\nnmacm> Error %s seems to be a not valid PDB\n", file_pdb );
		exit(0);
	}

	// Setting Coarse-Graining model
	Macromolecule *mol;
	switch(model)
	{
	case 0: // CA-model: CA + (NH and CO)-terminal model
		printf( "nmacm> Coarse-Graining model: CA-model\n");
		fprintf( f_com, "# nmacm> Coarse-Graining model: CA-model\n");
		// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
		// Warning, atoms are not copied, they're just pointers to the original atoms.
		// setmass = true --> adds masses to Occupancy and Bfactor, otherwise not left unchanged.
		// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
		//                   (unit mass is divided between NH or CO and their CAs)
		// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
		//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
		mol = cg_CA( molr, true, no_mass ); // masses will be computed
		break;

	case 4: // CA-only
		printf( "nmacm> Coarse-Graining model: CA-only\n");
		fprintf( f_com, "# nmacm> Coarse-Graining model: CA-only\n");
		// if CA selection
		mol = molr->select( calpha2 );
		mass_CAonly(mol, no_mass);
		break;

	case 3: // N,CA,C-model
	{
		printf( "nmacm> Coarse-Graining model: N,CA,C-model (Experimental)\n");

		// N,CA,C selection
		mol = molr->select( ncac2 );
		mass_NCAC( mol, no_mass, true ); // add masses
		break;
	}

	case 1: // 3BB2R model
		printf( "nmacm> Coarse-Graining model: 3BB2R\n");
		fprintf( f_com, "# nmacm> Coarse-Graining model: 3BB2R\n");
		mol = molr;
		if(!debug_3BB2R) // Makes 3BB2R model (if it's needed)
		{
			// CREATES a 3BB2R reduced model
			//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
			//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
			//     There are a few exceptions: for Ala, Gly and Pro,
			//     and for the 1st and last residues.
			printf("nmacm> Creating 3BB2R reduced model:\n");
//			num_res = mol->reducedmodel_3BBR2(no_mass, no_mass); // masses and charges on/off the same
			cg_3BBR2(mol, no_mass, no_mass); // masses and charges on/off the same
		}
		break;

	case 2: // Full-Atom
		printf( "nmacm> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
		fprintf( f_com, "# nmacm> Coarse-Graining model: Full-Atom (no coarse-graining)\n");
		mol = molr;
		mass_FA(mol, no_mass);
		break;
	}
	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();
	printf( "nmacm> Selected model number of residues: %d\n", num_res );
	printf( "nmacm> Selected model number of (pseudo)atoms: %d\n", num_atoms );
	fprintf( f_com, "# nmacm> Selected model number of residues: %d\n", num_res );
	fprintf( f_com, "# nmacm> Selected model number of (pseudo)atoms: %d\n", num_atoms );

	// Saving current model
	if(savemodel_switch)
	{
		sprintf(dummy_string,"%s_model.pdb",name);
		mol->writePDB( dummy_string );
	}

	// Setting model-dependent masses (needed below...)
	pdbIter *iter2 = new pdbIter( mol );
	pdbIter *iter2_res = new pdbIter( mol );
	Atom *atom;
	Residue *res;
	double *mass;

	printf("nmacm> Setting (pseudo)atomic masses:\n");
	if(no_mass && !(model==0 || model==1) ) // the CA-model has first NH and CA, and last CO and CA masses different.
	{
		printf("nmacm> No mass-weighting selected\n");
		fprintf(f_com,"# nmacm> No mass-weighting selected\n");
		mass = NULL;
	}
	else
	{
		// Allocating mass array memory
		if(!(mass = ( double * ) malloc( num_atoms * sizeof( double ) )))
		{
			printf("nmacm> Masses array memory allocation failed!\n");
			exit(1);
		}

		// Model dependent masses
		switch(model)
		{
//		case 4: // CA-only
//			printf( "nmacm> CA-only coarse-graining --> CA's will weight their residue mass.\n");
//			fprintf( f_com, "# nmacm> CA-only coarse-graining --> CA's will weight their residue mass.\n");
//			for ( iter2_res->pos_fragment = 0; !iter2_res->gend_fragment(); iter2_res->next_fragment() )
//			{
//				res = ( Residue * ) iter2_res->get_fragment();
//				mass[iter2_res->pos_fragment] = AA[ resnum_from_resname( res->getName() ) ].mass; // Residue mass
//				if( verb>0 )
//					printf("residue %4d  id= %3s  mass= %10.6f\n",iter2_res->pos_fragment,res->getName(), mass[iter2_res->pos_fragment]);
//			}
//			break;

		// Masses already set in occupancies!
		case 0: // CA-model
		case 1: // 3BB2R model
		case 2: // Full-Atom
		case 3: // NCAC model
		case 4: // CA-only
			printf( "nmacm> Atoms will weight their occupancy mass.\n");
			fprintf( f_com, "# nmacm> Atoms will weight their occupancy mass.\n");
			for ( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
			{
				atom = ( Atom * ) iter2->get_atom();
				mass[iter2->pos_atom] = atom->getPdbocc();
				if( verb>0 )
					printf("atom %4d  id= %3s  mass= %f10.6\n",iter2->pos_atom,atom->getName(),atom->getPdbocc() );
			}
			break;

//		case 2: // Full-Atom
//			printf( "nmacm> Full-Atom model --> Atoms will weight their atomic mass.\n");
//			fprintf( f_com, "# nmacm> Full-Atom model --> Atoms will weight their atomic mass.\n");
//			for ( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
//			{
//				atom = ( Atom * ) iter2->get_atom();
//				mass[iter2->pos_atom] = atom->getElement()->weight;
//				if( verb>0 )
//					printf("atom %4d  id= %3s  mass= %10.6f\n",iter2->pos_atom,atom->getName(),mass[iter2->pos_atom]);
//			}
//			break;
		}
	}
	delete iter2;
	delete iter2_res;

	// compute pairwise distance matrix
	double *dist_matrix;
	mol->distanceMatrix( &dist_matrix );

	// Checking min and max pair distances (Pablo's code)
	// (checking whether some atoms are too far from each other)
	double maxdist, mindist, currdist, currmindist;
	mindist = 1e20;
	maxdist = -1;
	for( i = 0; i < num_atoms; i++ )
	{
		currmindist = 1e20;
		for( j = 0; j < num_atoms; j++ )
		{
			if ( i == j )
				continue;
			else if ( i < j )
				index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2;
			else
				index = num_atoms * j + i - 1 - j * ( j + 3 ) / 2;
			currdist = dist_matrix[index];
			if ( currdist < currmindist ) currmindist = currdist;
		}
		if ( currmindist < mindist ) mindist = currmindist;
		if ( currmindist > maxdist ) maxdist = currmindist;
	}
	printf( "nmacm> Range of closest neighbor distances: %5.2f - %5.2f.\n", mindist, maxdist );
	fprintf( f_com, "# nmacm> Range of closest neighbor distances: %5.2f - %5.2f.\n", mindist, maxdist );

	// Initialize contact matrix
	double *cont_matrix;
	if( !(cont_matrix = ( double * ) malloc( sizeof( double ) * num_atoms * ( num_atoms - 1 ) / 2 ) ) )
	{
		printf("nmacm> Masses array memory allocation failed!\n");
		exit(1);
	}
	for ( i = 0; i <  num_atoms * ( num_atoms - 1 ) / 2 ; i++ )
		cont_matrix[i] = 0.0;

	// Opening output Force-Constants file
	FILE *f_Kout;
	if(Kout_switch)
	{
		sprintf(dummy_string,"%s_Kfile.dat",name);
		if( !(f_Kout = (FILE *)fopen(dummy_string,"w") ) )
		{
			printf("Sorry, unable to open Output contacts file: %s\nCheck input!!!\n",dummy_string);
			exit(1);
		}
	}

	// Setting up the force constants
	int nipa;
	switch(contacts)
	{
	case 0:
		//
		// POWER DISTANCE METHOD
		//
		printf( "nmacm> Power distance method (Inverse Exponential)\n");
		printf( "nmacm> \tk = %f\n", cte_k0);
		printf( "nmacm> \tcutoff %f\n", cutoff_k0 );
		printf( "nmacm> \tx0: %f\n", x0 );
		printf( "nmacm> \tpower: %f\n", power );
		fprintf( f_com, "# nmacm> Power distance method (Inverse Exponential)\n");
		fprintf( f_com, "# nmacm> \tk = %f\n", cte_k0);
		fprintf( f_com, "# nmacm> \tcutoff %f\n", cutoff_k0 );
		fprintf( f_com, "# nmacm> \tx0: %f\n", x0 );
		fprintf( f_com, "# nmacm> \tpower: %f\n", power );

		cont=0;
		for( i = 0; i < num_atoms; i++ )
			for( j = i+1; j < num_atoms; j++ ) //  (i < j, allways!)
			{
				index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ;
				if (dist_matrix[index] <= cutoff_k0)
				{
					cont_matrix[index] = inv_exp(cte_k0,dist_matrix[index],x0,power); // setting Force Constants
					if(Kout_switch)
						fprintf(f_Kout,"%d %d %E\n",i+1,j+1,cont_matrix[index]);
					cont++;
				}
			}
		break;

	case 1:
		//
		// DISTANCE CUTOFF METHOD
		//
		printf( "nmacm> Distance Cutoff method\n");
		printf( "nmacm> \tk = %f\n", cte_k1);
		printf( "nmacm> \tcutoff = %f\n", cutoff_k1);
		fprintf( f_com, "# nmacm> Distance Cutoff method\n");
		fprintf( f_com, "# nmacm> \tk = %f\n", cte_k1);
		fprintf( f_com, "# nmacm> \tcutoff = %f\n", cutoff_k1);
		if( cutoff_k1 < maxdist )
		{
			printf( "nmacm> Cutoff must be at least %5.2f to ensure that all atoms have at least one connected neighbor\n", maxdist );
			printf( "nmacm> (Every atom should be connected by at least 3 springs to avoid extra null-models, so increase the cutoff a little bit more!)\n");
			exit( 1 );
		}

		cont=0;
		for( i = 0; i < num_atoms; i++ )
			for( j = i+1; j < num_atoms; j++ ) //  (i < j, allways!)
			{
				index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ;
				if ( dist_matrix[index] <= cutoff_k1)
				{
					cont_matrix[index] = cte_k1;
					if(Kout_switch) // Writting Force-Constants
						fprintf(f_Kout,"%d %d %E\n",i+1,j+1,cont_matrix[index]);
					cont++;
				}
			}
		break;

	case 2:
		//
		// READING CONTACTS form Force Constants file (Klist)
		//
		printf("nmacm> Reading Force Constants file %s\n",Kfile);
		fprintf(f_com,"# nmacm> Reading Force Constants file %s\n",Kfile);
		FILE *f_Ks;
		if( !(f_Ks = (FILE *)fopen(Kfile,"r") ) )
		{
			printf("Sorry, unable to open Force Constants file: %s\nCheck input!!!\n",Kfile);
			exit(1);
		}

		cont=0;
		while( fscanf(f_Ks,"%d %d %f",&i,&j,&cte) != EOF )
		{
			i--; // input indices run from 1,...,N
			j--;
			if( j > i )
			{
				if( verb>0 )
					printf("i= %4d  j= %4d  cte= %10.6f ADDED!\n",i,j,cte);
				cont_matrix[ num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ] = (double) cte;
				if(Kout_switch)
					fprintf(f_Kout,"%d %d %E\n",i+1,j+1,cte);
				cont++;
			}
		}
		printf("nmacm> %d contacts readed!\n",cont);
		fprintf(f_com,"# nmacm> %d contacts readed!\n",cont);
		fclose( f_Ks );
		break;

	case 3: // Laura's "Mixed" model
		// Making Interacting Pair of (non-virtual) Atoms (ipas)
		cont_matrix_mix(mol, cont_matrix, num_atoms, &nipa);
		printf("nmadm>\nnmadm> \"Mixed\" method (%d nipas)\nnmadm>\n",nipa);

		if(Kout_switch)
		for( i = 0; i < num_atoms; i++ )
			for( j = i+1; j < num_atoms; j++ ) //  (i < j, allways!)
			{
				index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ;
				if( cont_matrix[num_atoms * i + j - 1 - i * ( i + 3 ) / 2] != 0.0 )
					fprintf(f_Kout,"%d %d %E\n",i+1,j+1,cont_matrix[index]);
			}

		break;

	case 4: // Yang,...,Jernigan's "pfENM". PNAS. 2009
		// Making Interacting Pair of (non-virtual) Atoms (ipas)
		cont_matrix_pfENM(mol, cont_matrix, num_atoms, &nipa);
		printf("nmadm>\nnmadm> \"pfENM\" method (%d nipas)\nnmadm>\n",nipa);

		if(Kout_switch)
		for( i = 0; i < num_atoms; i++ )
			for( j = i+1; j < num_atoms; j++ ) //  (i < j, allways!)
			{
				index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ;
				if( cont_matrix[num_atoms * i + j - 1 - i * ( i + 3 ) / 2] != 0.0 )
					fprintf(f_Kout,"%d %d %E\n",i+1,j+1,cont_matrix[index]);
			}

		break;

	case 5: // Yang,...,Jernigan's "pfENM". PNAS. 2009
		// Making Interacting Pair of (non-virtual) Atoms (ipas)
		cont_matrix_Kovacs(mol, cont_matrix, num_atoms, &nipa);
		printf("nmadm>\nnmadm> \"pfENM\" method (%d nipas)\nnmadm>\n",nipa);

		if(Kout_switch)
		for( i = 0; i < num_atoms; i++ )
			for( j = i+1; j < num_atoms; j++ ) //  (i < j, allways!)
			{
				index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ;
				if( cont_matrix[num_atoms * i + j - 1 - i * ( i + 3 ) / 2] != 0.0 )
					fprintf(f_Kout,"%d %d %E\n",i+1,j+1,cont_matrix[index]);
			}

		break;

	default:
		printf("Please, introduce a valid Conection method to continue!!!\n\nForcing exit!\n\n");
		exit(1);
		break;
	}

	// Closing output Force-Constants file
	if(Kout_switch)
	{
		fclose(f_Kout);
		printf( "nmacm> Contacts found: %d --> Written to %s\n", cont, dummy_string);
		fprintf( f_com,"# nmacm> Contacts found: %d --> Written to %s\n", cont, dummy_string);
	}
	else
		printf( "nmacm> Contacts found: %d\n", cont );

	// Getting coordinates single row vector
	float *coord;
	mol->coordMatrix( & coord );

	// HESSIAN MATRIX COMPUTATION
	double *hess_matrix=NULL; // forces automatic memmory allocation
	int size = num_atoms * 3;

	if(size<=0)
	{
		printf("nmadm> Sorry invalid size (%d), forcing exit!\n",size);
		exit(1);
	}

	// Number of eigenvectors to be computed (we need to know "size" first)
    if(nevec_fact >= 1.0) // number of modes
    {
    	nevec = (int) nevec_fact;
		printf( "nmadm> Number of requested modes: %d\n", nevec);
    }
    else
    {
    	nevec = (int) (nevec_fact * size);
    	printf( "nmadm> Number of requested modes: %d (%.0f%)\n", nevec, nevec_fact*100);
    }
    // Checking
	if(nevec > size)
	{
		printf("nmadm> Sorry, more eigenvectors requested (%d) than available (%d), forcing maximum.\n",nevec,size);
		nevec = size;
	}
	else if(nevec <= 0) // checking
    {
    	printf("nmadm> Error, invalid number of eigenvectors requested %d (%f)!\nForcing exit!\n",nevec,nevec_fact);
    	exit(1);
    }

	// You only save what you compute!
	modes_saved = nevec;
	printf( "nmadm> Number of Saved modes: %d\n", modes_saved);

	// Setting the number of modes to be saved
	if(!null) // saving only not-null modes
	{
		if(nevec < size-6)
		{
			nevec += 6; // to save "nevec" not-null vectors
			printf( "nmacm> Number of computed modes (null-included): %d (%.0f%)\n", nevec, ((float)nevec/size)*100);
		}
		else
		{
			printf( "nmacm> More not-null modes requested (%d) than available (%d), saving the (not-null) maximum available.\n",modes_saved,size-6);
			nevec = size;
			modes_saved = size-6; // then set number of modes saved to the maximum available
		}
	}

	ht_timer.restart(); // Hessian computation timer
	double r[3], rsqu, dummy;
	int m, n;
	if(no_mass && !(model==0 || model==1) )
		printf("nmacm>\nnmacm> Computing Hessian Matrix\n");
	else
		printf("nmacm>\nnmacm> Computing Mass-Weighted Hessian Matrix\n");
	hessianC(&hess_matrix,mass,coord,dist_matrix,cont_matrix,size);
	free( cont_matrix );
	if(!deform_switch) // if no deformability computations, then not needed dist_matrix!
		free( dist_matrix);

	if(verb > 2)
	{
		printf("Hessian Matrix (H) hess_matrix:\n");
		for ( int i = 0; i < size; i++ )
		{
			for ( int j = 0; j < size; j++ )
			{
				if(j>=i)
					printf("%9.6f ",hess_matrix[i + j*(j+1)/2]);
				else
					printf("%9.6f ",0.0);
			}
			printf("\n");
		}
	}

	printf( "nmacm> Number of Saved modes: %d\n", modes_saved);
	printf( "nmacm> Hessian computation time: %s\n",ht_timer.print_time_sec());
	fprintf( f_com, "# nmacm> Hessian computation time: %s\n",ht_timer.print_time_sec());

	// HESSIAN DIAGONALIZATION
	// Diagonalization input
	double *eigval=NULL;
	double *eigvect=NULL;
	// Allocating eigenvalues (if requested by NULL)
	if( !(eigval = (double *) malloc( sizeof(double) * nevec) ) )
	{
		printf("Msg(diag_dspevx): I'm sorry, Eigenvalues memory allocation failed!\n"
				"Forcing exit!\n");
		exit(1);
	}
	for( i=0; i<nevec; i++)
		eigval[i] = 0.0;

	// Allocating eigenvectors (if requested by NULL)
	if( !(eigvect = (double *) malloc( sizeof(double) * size * nevec ) ) )
	{
		printf("Msg(diag_dspevx): I'm sorry, Eigenvectors memory allocation failed!\n"
				"Forcing exit!\n");
		exit(1);
	}
	for( i=0; i<nevec*size; i++)
		eigvect[i] = 0.0;

	ht_timer.restart();
//	printf("nmacm\nnmacm> Hessian Diagonalizing with MKL's dspevx (%d modes requested)\n",nevec);
	printf("nmacm\nnmacm> Hessian Diagonalizing with DSPEVX (%d modes requested)\n",nevec);
	fflush(stdout);
	diag_dspevx(hess_matrix, &eigval, &eigvect, size, 1, nevec);
	printf("nmacm> Diagonalization time: %s\n",ht_timer.print_time_sec());
	fprintf(f_com,"# nmacm> Diagonalization time: %s\n",ht_timer.print_time_sec());

	// Showing output...
	int max=15;
	if(max>nevec)
		max=nevec;
	printf("nmacm> Showing the first %d eigenvalues:\n",max);
	printf("nmacm>\nnmacm> %4s %12s\n","MODE","EIGENVALUE");
	for(i=0; i<max; i++)
		printf("nmacm> %4d %12.5e\n",i+1, eigval[i]);

	// Saving Mass-Weighted cartesian modes
	if(save_wcart) // mass != NULL --> some masses are != 1.0
	{
		sprintf(dummy_string,"%s_wcart.evec",name);
		if(!null)
			save_ptraj_modes(dummy_string, size, 6, modes_saved+6, eigval, eigvect, norm_modes); // Saves & Normalizes
		else
			save_ptraj_modes(dummy_string, size, 0, modes_saved, eigval, eigvect, norm_modes); // Saves & Normalizes
		printf( "nmacm> Saving %d modes in ptraj file %s\n",modes_saved,dummy_string);
		fprintf( f_com, "# nmacm> Saving %d modes in ptraj file %s\n",modes_saved,dummy_string);
	}

	// switch from mass-weighted coords to standard cartesian coords.
	if(save_cart)
	{
		// Un-mass weighting eigenvectors
		if(mass != NULL) // mass != NULL --> some masses are != 1.0
		{
			for(i=0; i<nevec; i++)
				for(j=0; j<num_atoms; j++)
					for(m=0; m<3; m++)
						eigvect[i*size+3*j+m] /= sqrt(mass[j]);
		}

		// Saving modes
		sprintf(dummy_string,"%s_cart.evec",name);
		if(!null)
			save_ptraj_modes(dummy_string, size, 6, modes_saved+6, eigval, eigvect, norm_modes); // Saves & Normalizes
		else
			save_ptraj_modes(dummy_string, size, 0, modes_saved, eigval, eigvect, norm_modes); // Saves & Normalizes
		printf( "nmacm> Saving %d modes in ptraj file %s\n",modes_saved,dummy_string);
		fprintf( f_com, "# nmacm> Saving %d modes in ptraj file %s\n",modes_saved,dummy_string);
	}

	// COMPUTING DEFORMABILITY
	if(deform_switch)
	{
		double *deform, *mob;
		double *bf;
		char *list;

		bf = (double *) malloc(num_atoms * sizeof(double));
		mol->get_Pdbfact(bf);

		ht_timer.restart();
		printf("nmacm> Computing Deformability: ");
		fflush(stdout);
		compute_def(num_atoms, coord, dist_matrix, eigval, eigvect, nevec, 6, &deform, &mob);
		printf("%s\n",ht_timer.print_time_sec());

		//saving mobility in PDB
		norm_defmob (num_atoms, deform);
		mol->exchange_Pdbfact(deform);
		sprintf(dummy_string,"%s_def.pdb",name);
		printf("nmacm> Saving deformability in: %s\n",dummy_string);
		mol->writePDB(dummy_string);

		//saving mobilitiy in PDB
		norm_defmob(num_atoms,mob);
		mol->exchange_Pdbfact(mob);
		sprintf(dummy_string,"%s_mob.pdb",name);
		printf("nmadm> Saving mobility in: %s\n",dummy_string);
		mol->writePDB(dummy_string);

		free( dist_matrix );
		free(deform);
		free(mob);
	}

	if(!no_mass)
		free(mass);
	free(eigvect);
	free(eigval);

	fclose(f_com);
	printf("nmacm> Saved COMMAND & LOG file: %s\n",text);
	printf("nmacm>\nnmacm> The End\nnmacm>\nnmacm>\n");
}

void parseOptions(int argc, char** argv)
{
	std::string temp;
	CmdLine cmd("nmacm","NMACM: Normal Mode Analysis in Cartesian Coordinate Space (CCS)\n"
			"a) Coarse-Graining models: CA-model, 3BB2R(5-atom), Full-Atom, NCAC-model, CA-only.\n"
			"b) Contacting rules: cutoff, inverse-exponential, K-list.\n"
			"c) Mass-weighted CCS modes computation.\n"
			"\n ", version );

	try {

		// Define required no labeled arguments
		//
		UnlabeledValueArg<std::string> pdb("pdb","pdb input file","default","pdb");
		cmd.add( pdb );

		// Define labeled arguments
		//
        ValueArg<int> Verb("","verb", "Verbose level (default=0, no verbose)",false,0,"int");
        cmd.add( Verb );

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

		SwitchArg NoNorm("", "nonorm","Disables (norm=1) eigenvector normalization. (Note this does not affect vector direction)", true);
		cmd.add( NoNorm );
		SwitchArg NoMass("","nomass", "Disables mass weighting (default=mass-weighted)", true);
		cmd.add( NoMass );
		SwitchArg SaveWCart("", "save_wcart","Save Mass-weighted Cartesian modes(default=disabled)", false);
		cmd.add( SaveWCart );
		SwitchArg KOut("", "save_Kfile","Save atom-pairwise force constants file (to be used with -K option) (default=disabled)", false);
		cmd.add( KOut );
		SwitchArg Null("","null", "Save null-modes as well (default=disabled)", false);
		cmd.add( Null );

	    ValueArg<std::string> KFile("K","Kfile", "Force constants input file (Format 3 cols. plain-text: <i-atom> <j-atom> <K>)\n"
	    		"It can be made with --save_Kfile option.",false,"Kfile","string");
        cmd.add( KFile );
        SwitchArg Def("d","deform", "Turn on deformability calculations", true);
        cmd.add( Def );
        ValueArg<std::string> Name("o","name", "Output files basename.",false,"nmadm_out","string");
        cmd.add( Name );
        ValueArg<float> Nevs("n","nevs", "Number of eigenvectors (and eigenvalues) to be computed and saved. (n>=1 --> integer, n<1 --> ratio from maximum available) (default=0.2)",false,50,"int");
        cmd.add( Nevs );
        ValueArg<int> Contact("P","potential", "Contact method: (default=0)\n"
        		"  0=inverse exponential (= k/(1+(x/x0)^p), if x < c, else k=0)\n"
        		"  1=cutoff (= k, if x < c, else k=0)\n"
        		"  2=Kfile (--Kfile is mandatory)\n"
        		"  3=MIX-model\n"
        		"  4=pfENM\n"
        		"  5=Kovacs inverse exponential",false,0,"int");
        cmd.add( Contact );
		ValueArg<int> Model("m","model", "Coarse-Graining model: 0=CA, 1=3BB2R, 2=Full-Atom, 3=NCAC(experimental), 4=CA-only. (default=2)",false,2,"int");
        cmd.add( Model );

        // Parse the command line.
		cmd.parse(argc,argv);

		strcpy(file_pdb,((temp=pdb.getValue()).c_str()));
        strcpy(name,((temp=Name.getValue()).c_str())); // Gets Basename
        strcpy(Kfile,((temp=KFile.getValue()).c_str())); // Gets Force constants file-name
		if(KOut.isSet()) Kout_switch = true;

		// this should not be here...
		FILE *f;
		if ( (f=fopen(file_pdb, "r"))==NULL)
		{
			printf( "\n  Error->Cannot open file '%s'\n\n", file_pdb);
			exit(1);
		} fclose(f);

		// Contacting method parameters
        power = k0_Power.getValue();
		cte_k0 = k0_Cte.getValue();
		x0 = k0_X0.getValue();
		cutoff_k0 = k0_Cutoff.getValue();
		cte_k1 = k1_Cte.getValue();
		cutoff_k1 = k1_Cutoff.getValue();

    	// Number of eigenvectors to be computed
        nevec_fact = Nevs.getValue();
        if(nevec_fact <= 0) // checking
        {
        	printf("nmadm> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",nevec_fact);
        	exit(1);
        }
        //        modes_saved = ModesSaved.getValue();

		contacts = Contact.getValue(); // Contact method
		if(contacts == 2 && !KFile.isSet() ) // checking
		{
        	printf("nmadm> Error, you should include a KFile (see: --Kfile)!\nForcing exit!\n");
        	exit(1);
		}
		if( KFile.isSet() )
			contacts = 2;

		if (NoMass.isSet())
			no_mass = true;

		if (Null.isSet())
			null=true;

        if (Def.isSet())
        	deform_switch = true;

	    if(NoNorm.isSet())
	      	norm_modes = false;

    	model = Model.getValue();
    	if(model == 0 && no_mass) // CA-only and no-mass weighting
	    	save_cart = false; // not needed (wcart = cart)

	    if(SaveWCart.isSet())
	    	save_wcart = true;

		if( Verb.isSet() )
		{
			verb = Verb.getValue();
			printf("Parser input: Verbose level: --verb = %d\n",verb);
		}

	} catch ( ArgException& e )
	{ std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
	//cmd.~CmdLine();
}
