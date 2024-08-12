/************************************************************************
*                           NMATOOL                                     *
*************************************************************************
* Program is part of the ADP package URL: http://sbg.cib.csic.es        *
* (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
* Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
*************************************************************************
*                                                                       *
*   Tool for NMA stuff                                                  *
*   (from old "irbtool")
*                                                                       *
*************************************************************************
* This program is free software; you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation; either version 2 of the License, or     *
* (at your option) any later version.                                   *
************************************************************************/

#include "libnma/include/libnma_io.h" // Mon's NMA library (Dihedral included)
#include "cmdl/CmdLine.h"
//#include "libnma/nma.h" // Mon's NMA library (Dihedral included)

#define FILE_NAME 300

char file_pdb1[FILE_NAME];
char file_pdb2[FILE_NAME];
// char name[FILE_NAME];
char file_ptraj[FILE_NAME];
bool weight = false; // enables mass weighting
bool switch_3BB2R = false; // enables PDB translation into 3BB2R model (5-atom model)
bool ca_evec = false; // enables CA components extraction from 3BB2R ptraj

using namespace TCLAP;
void parseOptions(int argc, char** argv);

int main( int argc, char * argv[] )
{
	// Parsing Input
	parseOptions(argc,argv);
	printf("Mon's IRB-tool\n");

	int num_atoms1,num_atoms2;

	// create residues from templates
	init_aminoacids();

	// reading pdb 1
	printf( "nmadm>  1) Reading PDB 1 file\n" );
	Macromolecule * molr1 = new Macromolecule( file_pdb1 );
	molr1->readPDB(file_pdb1);
//	molr1->check_info_PDB();
	molr1->info(stdout);
	num_atoms1 = molr1->get_num_atoms();

	// Extracting CA eigenvector components for CA coarse-graining model comparation
	if(ca_evec)
	{
		double *invec,*inval,*evec_CA;
		int numatoms, numvectors, numcomps;
		// Reads an entire ptraj file allocating memory
//		int read_ptraj( char * file, double **evect, double **evals, int *numatoms, int *n_vectors, int *n_components  )
		read_ptraj(file_pdb2, &invec, &inval, &numatoms, &numvectors, &numcomps);
		printf( "irbtool>  %d atoms, %d eigenvectors and %d components found in %s\n",numatoms, numvectors, numcomps,file_pdb2 );

		int num_res = molr1->get_num_fragments();
		// Cartesian Eigenvector ---> CA-cartesian eigenvector
		int size_CA = num_res*3;
		if( !(evec_CA=(double *) malloc( size_CA * sizeof(double) * numvectors ) ) )
		{
			printf("Error, unable to allocate CA Cartesian NM's memory!\n");
			exit(1);
		}

		// Extracting CA eigenvector components for CA coarse-graining model comparation
		//Iterador para recorrer atomos
		pdbIter *iter = new pdbIter( molr1 );
		pdbIter *iter1 = new pdbIter( molr1 );
		int atom_index=0;
		int res_index;
		Residue *res;
		for ( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		{
			res_index=iter->pos_fragment; // shorter
			res = (Residue *) iter->get_fragment();
			iter1 = new pdbIter( res );
			for ( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() ) // screen atoms
			{
				if(	iter1->pos_atom == 1 ) // CA position index
				{
					for(int n=0; n<numvectors; n++)
					{
						evec_CA[n*size_CA+res_index*3] = invec[n*numcomps+atom_index*3];
						evec_CA[n*size_CA+res_index*3+1] = invec[n*numcomps+atom_index*3+1];
						evec_CA[n*size_CA+res_index*3+2] = invec[n*numcomps+atom_index*3+2];
					}
				}
				atom_index++;
			}
		}

		// Saving & Normalizing CARTESIAN Eigenvectors with "ptraj" format (CA-model)
		printf( "irbtool> Saving and normalizing %d modes in ptraj file %s\n",numvectors,file_ptraj);
		// Saves & Normalizes eigenpairs according to "ptraj" format
		save_ptraj_modes(file_ptraj, size_CA, 0, numvectors, inval, evec_CA, true); // Saves & Normalizes

		exit(0);
	}

	if( switch_3BB2R )
	{
		printf( "irbtool>  Set-up pseudo atomic model: 3BB + 2R\n" );
		molr1->format_residues();
		int num_res = molr1->reducedmodel_3BBR2();
		printf( "irbtool> Saving %s 3BB2R model (%d res.) into %s\n",file_pdb1,num_res,file_pdb2);
		molr1->writePDB(file_pdb2);
		printf( "irbtool> Saving succeeded.\n");
		exit(0);
	}

	printf( "nmadm>  2) Reading PDB 2 file\n" );
	Macromolecule * molr2 = new Macromolecule( file_pdb2 );
	molr2->readPDB(file_pdb2);
//	molr2->check_info_PDB();
	molr2->info(stdout);
	num_atoms2 = molr2->get_num_atoms();

	// Some checking
	if(num_atoms1 != num_atoms2)
	{
		printf("I'm sorry; please check input... different number of atoms found! (1:%d 2:%d)\n",num_atoms1,num_atoms2);
		exit(1);
	}
	else
		printf("OK, both PDB-files number of atoms match! (1:%d 2:%d)\n",num_atoms1,num_atoms2);

	double *difvec;
	difvec = (double *) malloc( sizeof(double) * num_atoms1 * 3 );

	// Getting masses
	double *mass;
	mass = ( double * ) malloc( num_atoms1 * sizeof( double ) );
	pdbIter *iter1 = new pdbIter( molr1 );
	Atom *atom;
	for ( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() )
	{
		atom = ( Atom * ) iter1->get_atom();
//		printf("atom %d  id= %s  mass= %f\n",iter1->pos_atom,atom->getName(),atom->getPdbocc() );
		mass[iter1->pos_atom] = atom->getPdbocc();
	}

	Tcoor pos1,pos2;
//	pdbIter *iter1 = new pdbIter( molr1 );
	pdbIter *iter2 = new pdbIter( molr2 );
	int index=0;
	double massw;

	if(weight) // mass-weighting
	{ // See: (ec. 59) from Brooks et al. J.Comput.Biol 16,12, 1522-42 (1995)
		printf("Mass-weighting the cartesian displacement vector.\n");
		for( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() ) // screen atoms
		{
			iter2->next_atom();
			( iter1->get_atom() )->getPosition( pos1 ); // Load pos1
			( iter2->get_atom() )->getPosition( pos2 ); // Load pos1
			massw = sqrt( mass[iter1->pos_atom] ); // mass-weighting...
			difvec[ iter1->pos_atom * 3 ] = (pos2[0] - pos1[0]) * massw;
			difvec[ iter1->pos_atom * 3 + 1] = (pos2[1] - pos1[1]) * massw;
			difvec[ iter1->pos_atom * 3 + 2] = (pos2[2] - pos1[2]) * massw;
		}
	}
	else // standard cartesian displacement vector
		for( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() ) // screen atoms
		{
			iter2->next_atom();
			( iter1->get_atom() )->getPosition( pos1 ); // Load pos1
			( iter2->get_atom() )->getPosition( pos2 ); // Load pos1
			difvec[ iter1->pos_atom * 3 ] = pos2[0] - pos1[0];
			difvec[ iter1->pos_atom * 3 + 1] = pos2[1] - pos1[1];
			difvec[ iter1->pos_atom * 3 + 2] = pos2[2] - pos1[2];
		}

	double eigval=0.0;
	printf("Saving ptraj difference vector: %s\n",file_ptraj);
	save_ptraj_modes(file_ptraj, 3*num_atoms1, 0, 1, &eigval, difvec, true); // Saves & Normalizes

	return 0;
}

/*==============================================================================================*/
void parseOptions(int argc, char** argv)
{
        std::string temp;
        CmdLine cmd("irbtool","IRB-tool (ADP-Platform)", "1.03" );

   try {
       // Define required arguments not labeled
       UnlabeledValueArg<std::string> Pdb1("pdb1","Initial PDB input file","default","pdb1");
       cmd.add( Pdb1 );

       // Define required arguments not labeled
       UnlabeledValueArg<std::string> Pdb2("pdb2","Final PDB input file","default","pdb2");
       cmd.add( Pdb2 );

       // Define required arguments not labeled
       UnlabeledValueArg<std::string> Ptraj("ptraj","output ptraj file with displacement vector","default","ptraj");
       cmd.add( Ptraj );

        // Define labeled arguments
	    ValueArg<std::string> Name("","name", "Sets the output files base-name.",false,"nmadm_ptraj","string");
        cmd.add( Name );

//		SwitchArg Verb("", "verb","Enables Hessian visualization and comparation with Naive version", true);
//		cmd.add( Verb );

		SwitchArg Weight("", "weight","Enables Mass-weighting cartesian displacement vector (default= false)", false);
		cmd.add( Weight );

		SwitchArg Convert3BB2R("", "3BB2R","Standard PDB conversion into 5-atoms model 3BB2R (<pdb2> will be the output) (you should introduce something into <ptraj>, it will be ignored)", false);
		cmd.add( Convert3BB2R );

		SwitchArg CAevec("", "ca_evec","Given a 3BB2R model (in <pdb1>) and a 3BB2R input ptraj file (in <pdb2>), the CA-components will be saved into a CA ptraj file (in <ptraj>). The new ptraj will be normalized. The remaining options will be ignored!", false);
		cmd.add( CAevec );

		// Parse the command line.
        cmd.parse(argc,argv);

		// Getting the command line arguments.
        strcpy(file_pdb1,((temp=Pdb1.getValue()).c_str())); // Gets PDB file name
        strcpy(file_pdb2,((temp=Pdb2.getValue()).c_str())); // Gets PDB file name
        strcpy(file_ptraj,((temp=Ptraj.getValue()).c_str())); // Gets PDB file name

	    if(Weight.isSet())
	    {
	    	weight = true;
	    	printf("irbtool> Cartesian displacement vector will be mass weighted ( sqrt(mi) )!\n");
	    }

	    if(Convert3BB2R.isSet())
	    {
	    	switch_3BB2R = true;
	    	printf("irbtool> Converting standard PDB into a 3BB2R model PDB (5-pseudo-atom model)!\n");
	    }

	    if(CAevec.isSet())
	    {
	    	ca_evec = true;
	    	printf("irbtool> Extracting CA-components form a 3BB2R ptraj file.\n");
	    }

        } catch ( ArgException& e )
        { std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }
}
