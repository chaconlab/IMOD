/*************************************************************************
 *                 libnmafit's LIBRARY: libnmafit.cpp                    *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Fitting/Molphing-related functions.                                 *
 *   (needs "libnma")                                                    *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include "nmafit.h"

extern CRandomMersenne *rg; // Mersenne Twister global object
//CRandomMersenne * rg; // Mersenne Twister global object

// Projects PDB to VOLUME & Filters it apropriately (Pablo's emt_pdb2vol code)
// set model_3BB2R = false to project with Masses
// set model_3BB2R = true to project with Charges (default)
void pdb2map(Macromolecule *in_pdb, vlVolume **volume, float resolution, float grid_size, int filter_method, bool model_3BB2R)
{
	bool debug=false;
	float sigma_factor = 3; // from "emt_pdb2vol"
    vlVolume *vol,*vol2;
	int pad[3];
	float r_half, r_cut, sigm;

	sigm = resolution/2.0;  // sigma
    r_half = sigm*sqrt(log(2.0))/sqrt(3.0/2.0); // r_half (Takes into account 3D...)
   	r_cut = sqrt(3.0)*sigm; // r_cut

    if (filter_method==1)
    {
     if ( r_half*0.8/grid_size < 1.0 )
     { // this is an empirical error bound
		fprintf(stdout,"pdb2map> Increase the sampling (>=%.2f) and try again\n>",r_half*0.8);
		exit(0);
	 }
   	}
/*
  	    else // Theoric way...
         if ( r_half/grid_size < 1.0 )  {
    		fprintf(stdout," Increase the sampling (>=%.2f) and try again\n>",r_half);
    		exit(0);
  		 }
*/
	in_pdb->geoBox();

	if(debug)
	{
    	cout<<"pdb2map> Number of read atoms: "<<in_pdb->get_num_atoms()<<endl;
    	cout<<"pdb2map> PDB electronic density "<<in_pdb->eDensity()<<endl;
	}

   switch (filter_method)
   {
    case 0:
    	// Gaussian Expansion
	    if(debug) cout<<"pdb2map> Gaussian expansion"<<endl;
	    vol2=in_pdb->pdb2map_real(resolution,grid_size);
	    // Saving files
//	    vol2 = FOPS::crop( vol2, 1E-12 );
    break;

    case 1:
    	// Fourier filter
	    if(debug) cout<<"pdb2map> Fourier filter"<<endl;
		if(model_3BB2R)
	    	vol=in_pdb->fillVolumeNE_3BB2R(grid_size);
	    else
	    	vol=in_pdb->fillVolumeNE(grid_size);
	    //cout<<"pdb2map> Electronic density "<<FOPS::calc_total(vol)<<endl;
	    pad[0]= (int)( (vol->dim().x()+resolution/grid_size+2) /4.0);
	    pad[1]= (int)( (vol->dim().y()+resolution/grid_size+2) /4.0);
	    pad[2]= (int)( (vol->dim().z()+resolution/grid_size+2) /4.0);
	    vol2=FOPS::padVolume(vol,vlDim(pad[0],pad[1],pad[2]),false);

		FOPS::GaussFilter(vol2, resolution);

      //FOPS::ButterFilter(vol2, resolution, 2, 0);
      //FOPS::ButterBandFilter(vol2, 1/width, resol, 2);
      //vol2 = FOPS::crop( vol2, 1E-12 );
    break;

    case 2:
	     // Kernel filter
	     if(debug) cout<<"pdb2map> Kernel filter"<<endl;
		 if(model_3BB2R)
	    	 vol=in_pdb->fillVolumeNE_3BB2R(grid_size);
	     else
	    	 vol=in_pdb->fillVolumeNE(grid_size);
	     float *kernel;
	     int dim_vox;
	     //cout<<"pdb2map> Electronic Density "<<FOPS::calc_total(vol)<<endl;
	     FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,vol->units().x(),resolution, sigma_factor);
	     vol = FOPS::padVolume( vol , vlDim( (dim_vox-1)/2, (dim_vox-1)/2,  (dim_vox-1)/2));
	     vol2 = FOPS::convoluteK(vol,kernel,dim_vox);
	     //FOPS::crop( vol2, 1E-12);
	     float cx_hi, cy_hi, cz_hi;
	     FOPS::center_masses( vol2, & cx_hi, & cy_hi, & cz_hi );
	     if(debug) cout<<"pdb2map> Mass center "<<" "<< cx_hi <<" "<< cy_hi <<" "<<cz_hi<<endl;
    break;
  }

   if(debug)
   {
	   cout<<"pdb2map> Electron density "<<FOPS::calc_total(vol2)<<endl;
	   cout<<"pdb2map> Volume size "<< vol2->dim() << endl;
	   cout<<"pdb2map> Filtering succeeded."<<endl;
	   cout<<"pdb2map> "<< endl;
   }

	*volume = vol2;
}

// INVERSE PARABOLIC INTERPOLATION
// Returns the expected minimum displacement coordinate (abscissa, delta_x)
// Minimum = orig.(fx) - parab_interpol()
double parab_iterpol(double fx, double fa, double fb, double offset)
{ // Numerical Recipes, pag 402.
  // 10.2 Parabolic Interpolation and Brent's Method in One Dimension
	double fxfa,fxfb,d,q,dzpl;
	if(fx<fa && fx<fb)
	{
	    // parabolic fit weights Numerical recipies 402
	    fxfa=fx-fa;
	    fxfb=fx-fb;

		d = (offset*offset) * (fxfb-fxfa);
		q = 2.0 * offset *(fxfb + fxfa);

		if( fabs(q) > 1.0E-10)
		{
			dzpl = -(d / q); // expected minimum position (parabolic approx.)
			return(dzpl);
		}
		else
		{
			printf("Msg(parab_interpol): Sorry, unstable interpolation, denominator close to Zero!\n");
			exit(1);
		}
	}
	else
	  if(fa<fb)
	    return(-offset);
	  else
	    return(offset);
}

// Rotates and traslates a Macromolecule, checking whether its necessary!
// Set "n_chain" <0 to rotate the whole molecule (no chain selection!)
void rotrans(Macromolecule *mol, float x, float y, float z, float a, float b, float c, int n_chain)
{
	bool debug=false;
	EulerRot *rotEuler;
	Tcoor pos;
	float mov[3] = { x, y, z };
	float rot[3] = { a, b, c };
	bool rotation, traslation;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen chains
	pdbIter *iter_chain;

	rotation = (rot[0]!=0 || rot[1]!=0 || rot[2]!=0);
	traslation = (mov[0]!=0 || mov[1]!=0 || mov[2]!=0);

	if(debug) printf( "Msg(rotrans): Chain %d => %5.3f %5.3f %5.3f (Rot,º) %5.3f %5.3f %5.3f (Trans,A)\n"
					, n_chain, rot[0]*(180/M_PI), rot[1]*(180/M_PI), rot[2]*(180/M_PI)
					, mov[0], mov[1], mov[2] );

//	for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
//	{
//		atom = ( iter_com->get_atom() );
//		mta = atom->getPdbocc(); // Load mass...
//		mtot += mta;
//		atom->getPosition(pos);
//		/* Sum(mass*coord) before putting the CoM at 0 */
//		r[0] += mta * pos[0];
//		r[1] += mta * pos[1];
//		r[2] += mta * pos[2];
//	}
//	r[0] /= mtot;
//	r[1] /= mtot;
//	r[2] /= mtot;
//	if(debug) printf( "Msg(move_dihedralM): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, r[0], r[1], r[2] );
//	// shift CM of pdb_model to 0,0,0
//	for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
//	{
//		atom = ( iter_com->get_atom() );
//		atom->getPosition(pos);
//		pos[0] -= r[0];
//		pos[1] -= r[1];
//		pos[2] -= r[2];
//		atom->setPosition(pos);
//	}

//	float r[3] = { 0.0, 0.0, 0.0 };
//	float mta,mtot=0;
	Macromolecule *molec = new Macromolecule();
	Protein *prot= new Protein();
////	Chain *chain;
////	Atom *atom;
//	prot->add(chain);
//	molchain->add(protchain);
//	molchain->e
//	molchain->add()
//
//	protchain->removeAll();
//	delete molchain;

	if(n_chain >= 0) // if: it has more than one chain
	{
		iter->pos_chain = n_chain;
		prot->add( iter->get_chain() ); // adds chain to the protein (reference)
		molec->add( prot ); // adds protein to the macromolecule

		if( rotation && traslation )
		{ // ROT + TRAS
			if(debug) printf("Msg(rotrans): Rotation + Traslation Detected!\n");
			rotEuler = new EulerRot( rot[0], rot[1], rot[2] ); // radians
			molec->geoCenter( pos ); // gets Center Of Masses from PDB
			for(int n=0; n<3; n++)
				pos[n] = -pos[n];
			molec->moveAll( pos ); // brings PDB to its COM
			molec->applyAtoms( rotEuler ); // rotates
			for(int n=0; n<3; n++)
				pos[n] = -pos[n] + mov[n]; // back movement + SHAKE
			molec->moveAll( pos ); // places back the PDB after rotation
		}
		else
		{
			if( traslation )
			{ // TRAS
				if(debug) printf("Msg(rotrans): Traslation Detected Only!\n");
				molec->moveAll( mov ); // traslates PDB
			}
			else if( rotation )
			{ // ROT
				if(debug) printf("Msg(rotrans): Rotation Detected Only!\n");
				rotEuler = new EulerRot( rot[0], rot[1], rot[2] ); // radians
				molec->geoCenter( pos ); // gets Center Of Masses from PDB
				for(int n=0; n<3; n++) pos[n] = -pos[n];
				molec->moveAll( pos ); // brings PDB to its COM
				molec->applyAtoms( rotEuler ); // rotates
				for(int n=0; n<3; n++) pos[n] = -pos[n]; // back movement + SHAKE
				molec->moveAll( pos ); // places back the PDB after rotation
			}
			else
			{ // NONE
				if(debug) printf("Msg(rotrans): Null movement detected! Doing Nothing!\n");
			}
		}
		prot->removeAll(); // deletes references only!
		delete molec;

	}
	else // < 0  ==> then old method (full rigid body rotation)
	{
		if( rotation && traslation )
		{ // ROT + TRAS
			if(debug) printf("Msg(rotrans): Rotation + Traslation Detected!\n");
			rotEuler = new EulerRot( rot[0], rot[1], rot[2] ); // radians
			mol->geoCenter( pos ); // gets Center Of Masses from PDB
// BUG 26/5/2009 (We weren't placing the molecule on the COM before rotation!!!
			for(int n=0; n<3; n++) pos[n] = -pos[n];
			mol->moveAll( pos ); // brings PDB to its COM
			mol->applyAtoms( rotEuler ); // rotates
			for(int n=0; n<3; n++) pos[n] = -pos[n] + mov[n]; // back movement + SHAKE
			mol->moveAll( pos ); // places back the PDB after rotation
		}
		else
		{
			if( traslation )
			{ // TRAS
				if(debug) printf("Msg(rotrans): Traslation Detected Only!\n");
				mol->moveAll( mov ); // traslates PDB
			}
			else if( rotation )
			{ // ROT
				if(debug) printf("Msg(rotrans): Rotation Detected Only!\n");
				rotEuler = new EulerRot( rot[0], rot[1], rot[2] ); // radians
				mol->geoCenter( pos ); // gets Center Of Masses from PDB
				// BUG 26/5/2009 (We weren't placing the molecule on the COM before rotation!!!
				for(int n=0; n<3; n++) pos[n] = -pos[n];
				mol->moveAll( pos ); // brings PDB to its COM
				mol->applyAtoms( rotEuler ); // rotates
				for(int n=0; n<3; n++) pos[n] = -pos[n]; // back movement + SHAKE
				mol->moveAll( pos ); // places back the PDB after rotation
			}
			else
			{ // NONE
				if(debug) printf("Msg(rotrans): Null movement detected! Doing Nothing!\n");
			}
		}
	}
	delete( iter );
}

// Rotates and traslates a Macromolecule, checking whether its necessary!
void rotrans_old(Macromolecule *mol, float x, float y, float z, float a, float b, float c)
{
	bool debug=true;
	EulerRot *rotEuler;
	Tcoor pos;
	float mov[3] = { x, y, z };
	float rot[3] = { a, b, c };
	bool rotation, traslation;

	rotation = (rot[0]!=0 || rot[1]!=0 || rot[2]!=0);
	traslation = (mov[0]!=0 || mov[1]!=0 || mov[2]!=0);

	if(debug)	printf( "Msg(rotrans): %5.3f %5.3f %5.3f (Rot,º) %5.3f %5.3f %5.3f (Trans,A)\n"
			, rot[0]*(180/M_PI), rot[1]*(180/M_PI), rot[2]*(180/M_PI)
			, mov[0], mov[1], mov[2] );

	if( rotation && traslation )
	{ // ROT + TRAS
		if(debug) printf("Msg(rotrans): Rotation + Traslation Detected!\n");
		rotEuler = new EulerRot( rot[0], rot[1], rot[2] ); // radians
		mol->geoCenter( pos ); // gets Center Of Masses from PDB
// BUG 26/5/2009 (We are not placing the molecule on the COM before rotation!!!
		mol->moveAll( pos ); // brings PDB to its COM
		mol->applyAtoms( rotEuler ); // rotates
		for(int n=0; n<3; n++) pos[n] = -pos[n] + mov[n]; // back movement + SHAKE
		mol->moveAll( pos ); // places back the PDB after rotation
	}
	else
	{
		if( traslation )
		{ // TRAS
			if(debug) printf("Msg(rotrans): Traslation Detected Only!\n");
			mol->moveAll( mov ); // traslates PDB
		}
		else if( rotation )
		{ // ROT
			if(debug) printf("Msg(rotrans): Rotation Detected Only!\n");
			rotEuler = new EulerRot( rot[0], rot[1], rot[2] ); // radians
			mol->geoCenter( pos ); // gets Center Of Masses from PDB
// BUG 26/5/2009 (We are not placing the molecule on the COM before rotation!!!
			mol->moveAll( pos ); // brings PDB to its COM
			mol->applyAtoms( rotEuler ); // rotates
			for(int n=0; n<3; n++) pos[n] = -pos[n]; // back movement + SHAKE
			mol->moveAll( pos ); // places back the PDB after rotation
		}
		else
		{ // NONE
			if(debug) printf("Msg(rotrans): Null movement detected! Doing Nothing!\n");
		}
	}
}

void build_tree(tree *current, int level, int max_level, float *prof, int *p_index, int max_index)
{
	int index = *p_index;
	// Some initializations...
	current->available = true; // allways available during initialization
	current->sumr = 0; // Sum values will be updated afterwards, see: init_tree()
	current->suml = 0;
	if( level < max_level )
	{
		if(index < max_index)
		{
			current->left = (tree *)malloc( sizeof(tree) ); // Allocate child's memory
			current->left->parent = current; // Tells child who is its parent
			build_tree(current->left,level+1,max_level,prof,p_index,max_index); // Propagates tree
		}
		else
		{  // por si acaso...
			current->left = NULL;
			return;
		}
		if(index < max_index)
		{
			current->right = (tree *)malloc( sizeof(tree) );
			current->right->parent = current;
			build_tree(current->right,level+1,max_level,prof,p_index,max_index);
		}
		else
		{
			current->right = NULL;
			return;
		}
	}
	else // if we are in the last level (adding leaf)
	{
		// Setting leaf values
		current->left = NULL;
		current->right = NULL;
		current->index = index;
		current->suml = prof[index];
		current->sumr = prof[index];
		(*p_index)++; // counts the number of leaves added so far (global-like var.)
	}
}

void init_tree(tree *current, int level, int max_level, float *prof, int *p_index, int max_index)
{
	int index = *p_index;
	if( level < max_level )
	{
		// Tree direct-propagation (from root to leaves)
		if(index < max_index)
			init_tree(current->left,level+1,max_level,prof,p_index,max_index);
		else
			return;

		if(index < max_index)
			init_tree(current->right,level+1,max_level,prof,p_index,max_index);
		else
			return;
	}
	else // if we are in the last level (adding leaf)
	{
		tree *p_now,*p_down;
		float weight;
		weight = current->suml;
		p_down = current;

		// Back-propagation (from leaf to the root node)
		for(int i = (max_level-1); i >= 0; i--)
		{
			// Going down
			p_now = p_down;
			p_down = p_now->parent;

			// Computing partial sums
			if(p_down->left == p_now) // Checking left/right branch
				p_down->suml += weight; // Adding leaf left-branch contribution
			else
				p_down->sumr += weight; // Adding leaf right-branch contribution
		}
		(*p_index)++; // counts the number of leaves added so far (global-like var.)
	}
}

// Weigted Random Sampling - WithOut Replacement (WRS-WOR)
// Based on the Wong & Easton (1980)'s method taken from Olken & Rotem (1995)
// (allocates memory to store the selected modes array, if *modes==NULL)
// profile (float)
void select_mode_WOR(float *profile, int size, int num, mode **modes, float cutoff)
{
	mode *p_modes;
	tree root; // allocating root
	tree *current = &root;
	int rank,rank_i,item_i;

	if(*modes == NULL) // Allocate memory
		if( !(p_modes=(mode *)malloc(sizeof(mode) * num)) )
		{
			printf("Unable to allocate memory (select_mode_WOR)\n");
			exit(1);
		}
	else // Don't allocate memory
		p_modes = *modes;

	// Determining tree rank
	rank = (int) ceil( log((float)size)/log((float)2) ); // binary logarithm (base 2)
	rank_i = rank; // rank index

	// Building up tree
	item_i = 0;
	build_tree(current,0,rank,profile,&item_i,size);

	// Intializing weight sums
	item_i = 0;
	init_tree(current,0,rank,profile,&item_i,size);

	// Sampling modes

}

// It selects "num" modes according to a given probability profile,
// and it sets their weights too!
// "Without replacement" was done DELETING each selected item from list.
// (allocates memory to store the selected modes array, if *modes==NULL)
// profile --> (float *)
// eigval --> (double *) Eigenvalues array for SCV weighting...
//void select_mode_prob(float *profile, int size, int num, mode **modes, float cutoff, bool randweight)
void select_mode_prob(float *profile, double *eigval, int size, int num, mode **modes, float cutoff, bool randweight)
{
	bool debug=false;
	float dice;
	int i;
	float *prob,*prof;
	mode *p_modes;
//	float cutoff=1; // cutoff to take into account probabilities
	bool noselected;
	int j,n_mode,n_mode_i,chosen=0;
	double weight;
	int *index_modes;
	int curr_size=size;
//	double KbT = 0.00198717 * 300; // Kboltz * Tempereature(K)

	if( !(prob=(float *)malloc(sizeof(float) * size)) )
	{
		printf("Unable to allocate memory (select_mode_prob)\n");
		exit(1);
	}
	if( !(prof=(float *)malloc(sizeof(float) * size)) )
	{
		printf("Unable to allocate memory (select_mode_prob)\n");
		exit(1);
	}
	if( !(index_modes=(int *)malloc(sizeof(int) * size)) )
	{
		printf("Unable to allocate memory (select_mode_prob)\n");
		exit(1);
	}

	if(*modes == NULL)
	{
		if( !(p_modes=(mode *)malloc(sizeof(mode) * num)) )
		{
			printf("Unable to allocate memory (select_mode_prob)\n");
			exit(1);
		}
	}
	else
	{
		p_modes = *modes;
		for(int i=0; i<num; i++)
		{
			p_modes[i].n = -1;
			p_modes[i].weight = 0.0;
		}
	}

	// Buffer Prob. Profile Initialization
	for(i=0; i<size; i++)
	{
		prof[i] = profile[i]; // "prof" initialization
		index_modes[i] = i; // indices initialization
	}

	// First Cumulative probability
	prob[0]=prof[0];
	for(i=1; i<size-chosen; i++)
		prob[i] = prob[i-1] + prof[i];

	// Cumulative Probability Normalization
//	if(debug) printf("Normalized --> prob[i] (i=0;i<size;i++)\n");
//	for(i=0; i<size; i++)
//	{
////		prob[i] /= prob[size-1]; // normalizes the cumulative prob.
//		if(debug) printf("%7.4f ",prob[i]);
//	}
//	if(debug) printf("\n");

	for(int n=0;n<num;n++)
	{
		// NORMAL MODE SELECTION
//		dice = (float) rg->Random() * prob[curr_size-1]; // Playing dice! [0:1) (Mersenne)
		dice = rg->Random() * prob[curr_size-1]; // Playing dice! [0:1) (Mersenne)
		// Selecting mode
		for(i=0; i<curr_size; i++) // screens probability profile
			if( dice <= prob[i] )
			{
				n_mode = index_modes[i]; // selecting (mode index)
				n_mode_i = i; // current mode index
				i = size; // forces loop exiting
			}
		curr_size--;
		// Deleting selected mode from profiles
		for(i=n_mode_i; i<curr_size; i++) // starting from selected mode index
		{
			prof[i] = prof[i+1]; // updates probability
			index_modes[i] = index_modes[i+1]; // updates indices
		}
		// Cumulative probability
		prob[0] = prof[0];
		for(i=n_mode_i+1; i<curr_size; i++)
			prob[i] = prob[i-1] + prof[i]; // updates cumulated prob.
//		for(i=0; i<curr_size; i++)
//			prob[i] /= prob[curr_size-1]; // normalizes the cumulative prob.

		// Storing selected mode
		p_modes[n].n = n_mode; // selecting (mode index)
//		double Random();
//		 * Gives a floating point random number in the interval 0 <= x < 1.
		if(randweight) // random within [-1.0,+1.0)
			weight = (double) 2 * rg->Random() - 1; // Playing dice! (-1,+1) (Mersenne)
		else // just random sign
			if(rg->Random() < 0.5)
				weight = 1.0;
			else
				weight = -1.0;

//		p_modes[n].weight = weight; // weighting

		// Apply the Scaling Colective Variable (SCV) Method (Go et al. 1985) to morphing protocol...
		if(eigval == NULL)
			p_modes[n].weight = weight; // weighting current mode like in SCV (KbT defined in this function)
		else
//			p_modes[n].weight = weight * sqrt(KbT/eigval[n_mode]); // weighting current mode like in SCV (KbT defined in this function)
			p_modes[n].weight = weight * sqrt(eigval[0]/eigval[n_mode]); // weighting current mode like in SCV (ENM independent)

		if(debug)
			printf("Msg(select_modes): # %3d: dice= %f  ==> Mode %4d ---> weight= %8.5f\n",n,dice,p_modes[n].n+1,p_modes[n].weight);
	}

	free(prob);
	free(prof);
	free(index_modes);
	*modes=p_modes; // outputs NM selection
	if(debug)
		printf("Msg(select_modes):  %d modes selected!\n",num);
}


// It selects "num" modes according to a given probability profile.
// "Without replacement" was done SCREENING the selected list to see whether new item matches.
// (allocates memory to store the selected modes array, if *modes==NULL)
// profile --> (float)
void select_mode_prob2(float *profile, int size, int num, mode **modes, float cutoff)
{
	bool debug=false;
	float dice;
	int nev; // chosen mode
	float *prob;
	mode *p_modes;
//	float cutoff=1; // cutoff to take into account probabilities
	bool noselected;
	int j,n_mode;
	double weight;

	if( !(prob=(float *)malloc(sizeof(float) * size)) )
	{
		printf("Unable to allocate memory (select_mode_prob)\n");
		exit(1);
	}

	if(*modes == NULL)
	{
		if( !(p_modes=(mode *)malloc(sizeof(mode) * num)) )
		{
			printf("Unable to allocate memory (select_mode_prob)\n");
			exit(1);
		}
	}
	else
	{
		p_modes = *modes;
		for(int i=0; i<num; i++)
		{
			p_modes[i].n = -1;
			p_modes[i].weight = 0.0;
		}
	}

	// Converting "profile" to an accumulative probability profile!
	if(debug) printf("Msg(select_modes): Converting to an accumulative probability! (normalized)\n");
//	prob[0]=0;
	for(int i=0; i<size; i++)
	{
		if(profile[i] > cutoff)
			if(i>0)
				prob[i] = prob[i-1] + profile[i] - cutoff; // accummulates probabilities
			else
				prob[i] = profile[i] - cutoff; // accummulates probabilities (1st time)
		else // zero probability
			if(i>0)
				prob[i] = prob[i-1];
			else
				prob[i] = 0;
	}
	for(int i=0; i<size; i++)
	{
		prob[i] /= prob[size-1]; // normalizes the accumulative prob.
		if(debug) printf("%7.4f ",prob[i]);
	}
	if(debug) printf("\n");

	for(int n=0;n<num;n++)
	{
		// NORMAL MODE SELECTION
//		dice = (float) rand()/RAND_MAX; // Playing dice! [0:1]
		dice = rg->Random(); // Playing dice! [0:1) (Mersenne)
		n_mode = 0; // we dont test the first one, so this it's needed!
					// if during last iter it doesn't match, then "nev=0" (1st evect)
		weight = 1;
		for(int i=0; i<size-1; i++)
			if( dice > prob[i] && dice <= prob[i+1] )
			{
				n_mode = i+1; // selecting (mode index)
				i = size; // forces loop exiting
			}

		// checks whether it has been already selected!
		noselected=true;
		j=0;
		while(noselected && j < n)
		{
			noselected = p_modes[j].n != n_mode; // noselected = false, when the mode was already selected
			j++;								  // ...then exits loop!
		}
		if(noselected) // If mode was not selected before, then...
		{
			p_modes[n].n = n_mode; // selecting (mode index)
//			p_modes[n].weight = (double) profile[i]; // weighiting
//			weight = (double) 2 * rand()/RAND_MAX - 1; // Playing dice! [-1,+1]
			weight = (double) 2 * rg->Random() - 1; // Playing dice! (-1,+1) (Mersenne)
			p_modes[n].weight = weight; // weighiting
			if(debug) printf("Msg(select_modes): dice= %f  ==> Mode %4d ---> weight= %8.5f\n",dice,p_modes[n].n+1,p_modes[n].weight);
		}
		else
		{
			if(debug) printf("Msg(select_modes): dice= %f  ==> Mode %4d ---> Already selected!\n",dice,n_mode+1);
			n--; // it tries again!
		}

	}

	free(prob);
	*modes=p_modes; // outputs NM selection
	if(debug) printf("Msg(select_modes):  %d modes selected!\n",num);
}

// It selects the "num" biggest "profile value" modes.
// (allocates memory to store the selected modes array, if *modes==NULL)
void select_mode_big(float *profile, int size, int num, mode **modes)
{
	bool debug=true;
	mode *p_modes;

	if(*modes == NULL)
	{
		if( !(p_modes=(mode *)malloc(sizeof(mode) * num)) )
		{
			printf("Unable to allocate memory (select_modes)\n");
			exit(1);
		}
	}
	else
		p_modes = *modes;

	int *indexes=NULL;
	sort_profile(profile, &indexes, size);

	for(int n=0;n<num;n++)
	{
		// NORMAL MODE SELECTION
		p_modes[n].n = indexes[n]; // selecting
//		p_modes[n].weight = (double)profile[ indexes[n] ]; // weighiting
		p_modes[n].weight = 1; // weighiting
		if(debug) printf("Msg(select_modes): Mode %4d ---> weight= %f\n",p_modes[n].n,p_modes[n].weight);
	}

	*modes=p_modes; // outputs NM selection
	if(debug) printf("Msg(select_modes):  %d modes selected!\n",num);
	free(indexes);
}

// It selects a peak according to a given probability profile.
int select_index(double *profile, int size, float cutoff)
{
	bool debug=true;
	float dice;
	double *prob;
	int nev;

	if( !(prob=(double *)malloc(sizeof(double) * size)) )
	{
		printf("Unable to allocate memory (select_index)\n");
		exit(1);
	}

	// Converting "profile" to an accumulative probability profile!
	if(debug) printf("Msg(select_index): Converting to an accumulative probability! (normalized)\n");
	prob[0]=0;
	for(int i=1; i<size; i++)
	{
		if( profile[i] > cutoff )
			prob[i] = prob[i-1] + profile[i] - cutoff; // it accummulates probabilities
		else
			prob[i] = prob[i-1]; // it accummulates probabilities
	}
	for(int i=0; i<size; i++)
	{
		prob[i] /= prob[size-1]; // normalizes the accumulative prob.
//		if(debug) printf("prob[%4d] = %f\n",i,prob[i]);
	}

	// SELECTION
//	dice = (double) rand()/RAND_MAX; // Playing dice!
	dice = rg->Random(); // Playing dice! [0:1) (Mersenne)
	nev = size-1; // we dont test the last one, so this it's needed!
						 // (if during last iter doesn't match, then "nev=last evect")
	for(int i=0; i<size-1; i++)
		if( dice > prob[i] && dice < prob[i+1] )
		{
			nev = i; // storing index
			i = size; // forcing exit
		}
	free(prob);
	if(debug) printf("Msg(select_index): Array index selected = %d\n",nev);
	return(nev);
}

//// It Merges "nev" modes, selected in "nm_props" (FAST)
//// (warning: not allocates memory!)
//// type = 0 --> phi,psi,...
//// type = 1 --> phi,psi,chi,...
//// type = 2 --> phi,chi,psi...
//void merge_modes(Macromolecule *mol, NM *uu, tri *props, double *hess_matrix, int size, mode *nm_props, int nev, int type)
//{
//	bool debug=false;
//	double weight;
//	int num_res = mol->get_num_fragments();
//	int j,cont;
//	double merged[size];
//	double *p_mode=NULL;
//	pdbIter *iter;
//	Residue *res;
//	double zero = ZERO; // -99999.0;
//
//	// Merged mode initialization
//	for( int k = 0; k < size; k++ )
//		merged[k] = 0.0;
//
//	if(debug) printf("Msg(merge_modes): Merging %d modes\n",nev);
//
//// Big BUG FIX
//	for(int k=0; k<num_res; k++)
//		uu[k].dphi = uu[k].dpsi = uu[k].dchi = 0.0;
//
//	for(int i=0; i<nev; i++)
//	{
//		if(debug) printf("Msg(merge_modes):\t%2d --> Mode: %4d  Weight: %6.4f\n",i+1,nm_props[i].n+1,nm_props[i].weight);
//
//		// Extracts Dihedral-Angle components from the Hessian Matrix (raw modes)
//		// Allocates the mode array, only for the first time.
////		dihedral_comps(mol,hess_matrix,nm_props[i].n,size,props,&p_mode);
//		dihedral_comps(hess_matrix,nm_props[i].n,size,&p_mode);
//
//		// Updating merged mode
//		weight = nm_props[i].weight;
//		for( int k=0; k<size; k++ )
//			merged[k] += p_mode[k] * weight; // also weighting
//	}
//
//	cont = 0;
//	// Outputs a "NM" struct
//	iter = (pdbIter *) new pdbIter(mol);
//	if(debug) printf("Msg(dihedral_comps): Finding Phi & Psi indices\n");
//	// Sorting TYPE 1 ( phi,psi,chi,.. ) DEFAULT
//	if(type == 1)
//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			j=iter->pos_fragment; // shorter
//			res = ( Residue * ) iter->get_fragment();
//
//			// PHI
//			if( j != 0 && (strcmp(res->getName(), "PRO") != 0) ) // != PRO
//			{
//				uu[j].dphi = merged[cont]; // storing PHI amplitude
//				cont++;
//			}
//			else
//				uu[j].dphi = zero; // loading PHI
//
//			// PSI
//			if(j!=num_res-1)
//			{
//				uu[j].dpsi = merged[cont]; // storing PSI amplitude
//				cont++;
//			}
//			else
//				uu[j].dpsi = zero; // loading PSI
//
//			// CHI
//			if(props[j].nan==3 || (props[j].nan==2 && (j==0 || j==num_res-1)))
//			{
//				uu[j].dchi = merged[cont]; // storing PSI amplitude
//				cont++;
//			}
//			else
//				uu[j].dchi = zero; // loading CHI
//		}
//	// Sorting TYPE 2 ( phi,chi,psi... )
//	if(type == 2)
//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			j=iter->pos_fragment; // shorter
//			res = ( Residue * ) iter->get_fragment();
//
//			// PHI
//			if( j != 0 && (strcmp(res->getName(), "PRO") != 0) ) // != PRO
//			{
//				uu[j].dphi = merged[cont]; // storing PHI amplitude
//				cont++;
//			}
//			else
//				uu[j].dphi = zero; // loading PHI
//
//			// CHI
//			if(props[j].nan==3 || (props[j].nan==2 && (j==0 || j==num_res-1)))
//			{
//				uu[j].dchi = merged[cont]; // storing PSI amplitude
//				cont++;
//			}
//			else
//				uu[j].dchi = zero; // loading CHI
//
//			// PSI
//			if(j!=num_res-1)
//			{
//				uu[j].dpsi = merged[cont]; // storing PSI amplitude
//				cont++;
//			}
//			else
//				uu[j].dpsi = zero; // loading PSI
//		}
//
//	// Sorting TYPE 0 ( phi,psi... )
//	if(type == 0)
//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			j=iter->pos_fragment; // shorter
//			res = ( Residue * ) iter->get_fragment();
//
//			// PHI
//			if( j != 0 && (strcmp(res->getName(), "PRO") != 0) ) // != PRO
//			{
//				uu[j].dphi = merged[cont]; // storing PHI amplitude
//				cont++;
//			}
//			else
//				uu[j].dphi = zero; // loading PHI
//
//			// CHI
//			if(props[j].nan==3 || (props[j].nan==2 && (j==0 || j==num_res-1)))
//			{
//				uu[j].dchi = zero; // reset CHI (it exists, but should not be accounted for)
//			}
//
//			// PSI
//			if(j!=num_res-1)
//			{
//				uu[j].dpsi = merged[cont]; // storing PSI amplitude
//				cont++;
//			}
//			else
//				uu[j].dpsi = zero; // loading PSI
//		}
//
//	delete iter;
//	free( p_mode );
//}

// It Merges "nev" modes, selected in "nm_props" (FAST)
// (warning: not allocates memory!)
void merge_modes(double *uu, double *hess_matrix, int size, mode *nm_props, int nev)
{
	bool debug=false;
	double weight;
	double *mode;

	// Merged mode memory allocation
	if( !(mode = (double *) malloc( sizeof(double) * size ) ) )
	{
		printf("Msg(merge_modes): Memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	// initializations
	for( int k = 0; k < size; k++ )
		uu[k] = 0.0;

	if(debug) printf("Msg(merge_modes): Merging %d modes\n",nev);

	for(int i=0; i<nev; i++)
	{
		if(debug) printf("Msg(merge_modes):\t%2d --> Mode: %4d  Weight: %6.4f\n",i+1,nm_props[i].n+1,nm_props[i].weight);

		// Extracts Dihedral-Angle components from the Hessian Matrix (raw modes)
		// Allocates the mode array, only for the first time.
//		dihedral_comps(mol,hess_matrix,nm_props[i].n,size,props,&p_mode,num_res);
//		get_mode(mode, hess_matrix, size, nm_props[i].n);

		// Extracts Dihedral-Angle components directly from the Hessian Matrix raw modes
		// Allocates the mode array, only for the first time.
		// (Please, free it when not needed)
		dihedral_comps(hess_matrix,nm_props[i].n,size,&mode);

		// Updating merged mode
		weight = nm_props[i].weight;
		for( int k=0; k<size; k++ )
			uu[k] += mode[k] * weight; // also weighting
	}
	free( mode );
}

// Randomly wheights the already present NM weight.
void rand_mode_weight(mode *modes, int num)
{
	bool debug=true;
	double weight;

	for(int n=0;n<num;n++)
	{
//		weight = (double) 2 * rand()/RAND_MAX - 1; // [-1,+1]
		weight = (double) 2 * rg->Random() - 1; // Playing dice! (-1,+1) (Mersenne)
		modes[n].weight *= weight; // weighiting
		if(debug) printf("Msg(rand_mode_weight): Mode index %4d ---> weight= %f\n",modes[n].n,modes[n].weight);
	}
}

//// Computes dihedral angle diferences between two macromolecules
//// The result will be stored inside "NM" struct
//void dihedral_diff(Macromolecule *mol1, Macromolecule *mol2, NM *uu)
//{
//	float *dihedrals1,*dihedrals2;
//	int dh_comps=7; // (phi,psi,omega,4x Chi)
//	int num_res = mol1->get_num_fragments();
//	int num_res2 = mol2->get_num_fragments();
//
//	if(num_res != num_res2)
//	{
//		printf("Msg(dihedral_diff): Sorry, residue numbers must agree!!!\nForcing Exit\n");
//		exit(1);
//	}
//
//	// Dihedral angles calculation
//	mol1->all_dihedrals( &dihedrals1 );
//	mol2->all_dihedrals( &dihedrals2 );
//
//	for(int i=0;i<num_res;i++)
//	{
//		for(int n=0;n<dh_comps;n++)
//		{
//			printf("%9.6f ",dihedrals1[i*dh_comps+n] - dihedrals2[i*dh_comps+n]);
//		}
//		printf("\n");
//	}
//}

// Sorts indexes from a profile based on its profile values (big-->small)
// (allocates memory to store the selected modes array, if *indexes==NULL)
void sort_profile(float *profile, int **indexes, int size)
{
	int *p_indexes;
	bool changed=true;
	float dummy[size];
	float tempf;
	int tempi;

	if( *indexes == NULL )
		if( !(p_indexes = (int *)malloc( sizeof(int) * size ) ) )
		{
			printf("Msg(sort_profile): Unable to allocate memory!\n");
			exit(1);
		}

	// Indexes array initialization
	for(int i=0;i<size;i++)
	{
		p_indexes[i] = i;
		dummy[i] = profile[i];
	}

	// Sorting
	while(changed)
	{
		changed=false;
		for(int i=0;i<size-1;i++)
		{
			if( dummy[i] < dummy[i+1] )
			{ // swapping
				tempf = dummy[i+1];
				dummy[i+1] = dummy[i];
				dummy[i] = tempf;
				tempi = p_indexes[i+1];
				p_indexes[i+1] = p_indexes[i];
				p_indexes[i] = tempi;
				changed=true;
			}
		}
	}
	*indexes = p_indexes; // output
}

// Overwrites the "mol_ref" atomic coordinates over "mol" ones.
void replace_pos(Macromolecule *mol, Macromolecule *mol_ref)
{
	int num,num_ref;
	pdbIter *iter,*iter_ref;
	Tcoor pos_ref;

	iter = new pdbIter(mol);
	iter_ref = new pdbIter(mol_ref);

	// Some checking
	num = iter->num_fragment();
	num_ref = iter_ref->num_fragment();
	if(num != num_ref)
	{
		printf("Msg(replace_pos): Number of atoms mismatch! Forcing Exit!\n");
		exit(1);
	}

	for( iter->pos_atom=0, iter_ref->pos_atom=0 ; !iter->gend_atom(); iter->next_atom(), iter_ref->next_atom() )
	{
		(iter_ref->get_atom())->getPosition(pos_ref); // get reference position
		(iter->get_atom())->setPosition(pos_ref); // set position
	}
}

// Projects a Macromolecule (mol) into a volume (pvol) without unnecessary memory allocations.
// The convolution Kernel and Two padded volumes must be already allocated (pvol & pdummy).
// ("dim_vox" is only needed by "Gausian filter", otherwise it will be ignored)
// nthreads --> Number of threads used in parallel. (=0 for non-parallel)

#ifdef USE_PTHREAD // Enables PThread parallel routines
void project_pdb(Macromolecule *mol,vlVolume **pvol,vlVolume **pdummy,float *kernel,int dim_vox,int method, float res, bool fast,
				int nthreads, convoluteK_data *threads_data)
#else
void project_pdb(Macromolecule *mol,vlVolume **pvol,vlVolume **pdummy,float *kernel,int dim_vox,int method, float res, bool fast)
#endif
{
	bool debug = false; // = true --> Enables screen output and timers
	vlVolume *temp;
//	Htimer ht_timer; // timer
	timerReal rt_timer1,rt_timer2; // Rapid timers!
	float timer1,timer2;

//	if(debug)
//		ht_timer.restart();

#ifdef USE_PTHREAD // Enables PThread parallel routines
	if(debug)
		printf ("Msg(project_pdb):  Projecting Model into Map (nthreads=%d, method=%d) ",nthreads,method);
#endif

	(*pvol)->clear(0);

	// PROJECTION AND FILTRATION
	if(debug)
		rt_timer1.startTimer();
//	if(nthreads > 0 && threads_data_proj != NULL ) // If parallel...
//	{
////		fprintf(stderr,"before project_RealNE_3BB2R_par...\n");
//		mol->project_RealNE_3BB2R_par( *pvol, nthreads, threads_data_proj );
////		fprintf(stderr,"after project_RealNE_3BB2R_par...\n");
//	}
//	else
	mol->project_RealNE_3BB2R( *pvol, fast ); // projects without moving or creating volume either

	if(debug)
	{
		rt_timer1.stopTimer();
		timer1 = (float) rt_timer1.getElapsedTime();
		rt_timer2.startTimer();
	}

	switch(method)
	{
		case 1:
		    FOPS::GaussFilter( *pvol, res);
		    break;

		case 2:
			if(debug)
			{
//				printf("(%s)\n",ht_timer.print_time());
//				ht_timer.restart();
				printf ("Msg(project_pdb):  Kernel convolution ");
				fflush(stdout);
			}

#ifdef USE_PTHREAD  // If PThread  enabled
			if(nthreads > 0 && threads_data != NULL ) // If parallel...
			{
////				fprintf(stderr,"Ejecutando: convoluteK_nopad_par... with %d threads  dim_vox= %d\n",nthreads,dim_vox);
//#ifdef USE_TBB // Enables Intel's Threading Building Blocks parallel routines
//				*pdummy = FOPS::convoluteK_nopad_tbb(*pvol,*pdummy,kernel,dim_vox,nthreads);
//#endif
//#ifdef USE_PTHREAD // Enables PThread parallel routines
				*pdummy = FOPS::convoluteK_nopad_par(*pvol,*pdummy,kernel,dim_vox,nthreads,threads_data);
//#endif
			}
			else
			{
				// fprintf(stderr,"before convoluteK_nopad\n");
				//				fprintf(stderr,"Ejecutando: convoluteK_nopad... %d threads  dim_vox= %d\n", nthreads,dim_vox);
				*pdummy = FOPS::convoluteK_nopad(*pvol,*pdummy,kernel,dim_vox);
				//fprintf(stderr,"after convoluteK_nopad\n");
				//FOPS::writeFile(*pdummy,"pdummy.sit");
			}
#else
			*pdummy = FOPS::convoluteK_nopad(*pvol,*pdummy,kernel,dim_vox);
#endif

			temp = *pvol;
			*pvol = *pdummy;
			*pdummy = temp;
		break;

		default:
			printf("Msg(project_pdb): Please select a valid PDB projection method!\n");
			exit(1);
		break;
	}
//FOPS::writeFile(*pvol,"pvol.sit");
	if(debug)
	{
		rt_timer2.stopTimer();
		timer2 = (float) rt_timer2.getElapsedTime();
		fprintf(stderr,"projection= %f s   convolution= %f\n",timer1,timer2);
	}
//	if(debug)
//		printf("(%s)\n",ht_timer.print_time());
}

// Accepts or rejects a candidate based on Simulated Annealing
// "next" should be higher than "actual"
bool anneal( double actual, double next, double temp )
{
	double dice,prob;
	// Kirkpatrick et al.
	// if ep < e --> P = 1
	// else P = exp( (e-ep)/T )

	if(temp > 0.0)
	{
		if( next > actual )
			prob = exp( (actual - next) / temp );
		else
			prob = 1;
	}
	// This fixes some minor "up-hill" peaks in the score profile (22/9/2009)
	else // if temp <= 0.0
	{
		if( next >= actual )
			return false;
		else
			return true;
	}

//	dice = (double) rand()/RAND_MAX; // Playing dice! [0,1]
	dice = rg->Random(); // Playing dice! [0:1) (Mersenne)
//	printf("Msg(anneal): actual=%6.4f next=%6.4f T=%4e dice=%5.3f prob=%4e\n",actual,next,temp,dice,prob);
	if( prob >= dice )
	{
//		printf("Msg(anneal): ACCEPTED!!!\n");
		return true;
	}
//	printf("Msg(anneal): REJECTED!!!\n");
	return false;
}

// Projects a Map over a Macromolecule, and outputs a profile with the densities
// within a cube from each residue atom. (The profile has N-residues elements)
// radius = sphere radius in Amstrongs (defines neighborhood)
// size = kernel size (to save time)
// Profile memory is allocated if *p_profile=NULL.
void project_map2pdb(vlVolume *map, Macromolecule *mol,double **p_profile,float radius,int size)
{
    int i,j,z,vi,vj,vz,di,dj,dz;
    int cont=0;
    float voxel;
    int stepx,stepy,stepz,offset;
    int border= (size-1)/2; // respect to the kernell size
    Tcoor at_pos; // ,vox_pos; // positions
    int index=0;
	pdbIter *iter_atom,*iter_res;
	Residue *res;
	double *profile;

	if(*p_profile == NULL)
	{
		profile = (double *) malloc(sizeof(double) * mol->get_num_fragments());
		check_pointer(profile,"project_map2pdb --> allocating profile");
		*p_profile = profile;
	}
	else
		profile = *p_profile;

	vlStep step=map->stepping();
    stepx=step.x();
    stepy=step.y();
    stepz=step.z();
    vlDim dim=map->dim();
    float *container= (float*)map->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    int rang_vox = (int) ceil( radius / map->units().x() ); // // maximum distance considered (# voxels)

	vlPoint3f pos;
	map->getPosition(&pos); // map corner position
	int num_res = mol->get_num_fragments();
//	FOPS::writeFile(map,"map2pdb.sit");

	// SCREENING MACROMOLECULE FIRST
	iter_res = new pdbIter(mol);

	for(iter_res->pos_fragment=0; !iter_res->gend_fragment(); iter_res->next_fragment())
	{
		profile[index] = 0.0; // reset residue value
		cont = 0; // reset residue voxel counter
		res = (Residue *) iter_res->get_fragment();
		iter_atom = new pdbIter( res );
		for( iter_atom->pos_atom=0; !iter_atom->gend_atom(); iter_atom->next_atom() )
		{
			// gets atom position (PDB)
			(iter_atom->get_atom())->getPosition(at_pos);
			// interpolates the atom position in the map (gets voxel value)

			// Atom's voxel indices
			vi = (int) roundf( ( at_pos[0] - pos.x() ) / map->units().x() );
			vj = (int) roundf( ( at_pos[1] - pos.y() ) / map->units().y() );
			vz = (int) roundf( ( at_pos[2] - pos.z() ) / map->units().z() );

		    	// Screens atom's voxel bounding box
		    for(i=vi-rang_vox;i<=vi+rang_vox;i++)
		      for(j=vj-rang_vox;j<=vj+rang_vox;j++)
		        for(z=vz-rang_vox;z<=vz+rang_vox;z++)
		        {
		        	di = i-vi;
		        	dj = j-vj;
		        	dz = z-vz;
		        	// checks whether the atom is inside the voxel bounding sphere (integer distance; less accurate, but faster!!!)
		        	if( (int)sqrtf( di*di + dj*dj + dz*dz ) <= rang_vox ) // if inside bounding sphere
			        if( i>=border && j>=border && z>=border && i<dim.x()-border && j<dim.y()-border && z<dim.z()-border )
			        { // if inside map
				    	  offset = (i*stepx)+(j*stepy)+(z*stepz);
				    	  voxel= *(container+offset); // voxel density
					  	  profile[index] += voxel;
		        		  cont++; // counts residue voxel number
		        	}
		        }

		}
		delete iter_atom;
		if(cont!=0) // avoids zero division!!!
			profile[index] /= cont; // Computes Average
		else
			profile[index] = 0.0;
		index ++;
	}
	delete iter_res;
}

// Smooths a profile, allocating memory for the new profile.
// method = 0 --> Windowed average
// method = 1 --> Gaussian Kernel
double *smooth_profile(double *profile, int size, int kernel, int method)
{
	double avg;
	int cont;
	double *profile_new = (double *) malloc( sizeof(double) * size);

	switch(method)
	{
		case 0: // Windowed average
			for( int i = 0; i < size; i++ )
			{
				avg = 0.0;
				cont = 0;
				for( int j = i - kernel; j < i + kernel; j++ )
					if( j >= 0 && j < size )
					{
						avg += profile[j];
						cont++;
					}
				avg /= cont;
				profile_new[i] = avg;
			}
		break;

		case 1:
		break;
	}
	return profile_new;
}

// Normalizes a profile, (overwriting it!)
// New values will range from "min" to "max" (lineal transformation)
void norm_profile(double *profile, int size, double min_new, double max_new)
{
	double min=1e10, max=-1e10, factor, factor2;

	// Max and Min values...
	for( int i = 0; i < size; i++ )
	{
		if(profile[i] > max)
			max = profile[i];
		if(profile[i] < min)
			min = profile[i];
	}

	// Normalization
	factor = max - min;
	factor2 = max_new - min_new;
//	if(abs(factor2) > 1e-6 && abs(factor) > 1e-6) // Avoids division by zero
	if(abs(factor) > 1e-6) // Avoids division by zero
		for( int i = 0; i < size; i++ )
			profile[i] = ( (profile[i] - min) / factor ) * factor2 + min_new;
}
void norm_profile(float *profile, int size, float min_new, float max_new)
{
	float min=1e10, max=-1e10;
	double factor, factor2;

	// Max and Min values...
	for( int i = 0; i < size; i++ )
	{
		if(profile[i] > max)
			max = profile[i];
		if(profile[i] < min)
			min = profile[i];
	}

	// Normalization
	factor = max - min;
	factor2 = max_new - min_new;
//	if(abs(factor2) > 1e-6 && abs(factor) > 1e-6) // Avoids division by zero
	if(abs(factor) > 1e-6) // Avoids division by zero
		for( int i = 0; i < size; i++ )
			profile[i] = ( (profile[i] - min) / factor ) * factor2 + min_new;
}

// Computes normalized cross-correlation (from two vectors)
double corr_vector( double *a, double *b, int size )
{
	double corr = 0, avg_a = 0, avg_b = 0, sig_a = 0, sig_b = 0;
	// Computing averages
	for ( int i = 0; i < size ; i++ )
	{
		avg_a += *( a + i );
		avg_b += *( b + i );
	}
	avg_a /= size;
	avg_b /= size;

	// Computing sigmas
	for ( int i = 0; i < size ; i++ )
	{
		sig_a += pow( *( a + i ) - avg_a, 2);
		sig_b += pow( *( b + i ) - avg_b, 2);
	}
	sig_a /= size;
	sig_b /= size;
	sig_a = sqrt( sig_a );
	sig_b = sqrt( sig_b );

	// Computing normalized cross-correlation
	for(int i=0; i<size; i++)	corr += ( *(a + i) - avg_a ) * ( *(b + i) - avg_b );
	corr /= size * sig_a * sig_b;
	return( corr );
}

// B-factor profile computation function (Using ADP library)
// Cartesian eigenvectors needed.
// It does not allocate "profile" memory. "nm"-->[0,#modes-1]
// ("CA" atom only!)
void bfact_profile(Macromolecule *mol, double *evec, double *profile, int nm)
{
	int num_atoms = mol->get_num_atoms();
	int nm_index=3*num_atoms*nm; // <-- mode number index to compute B-fact
	int atom_index=0,res_index=0;
	pdbIter *iter_res,*iter_atom;
	Tcoor pos;
	Residue *res;

	iter_res = new pdbIter(mol);
	for(iter_res->pos_fragment=0; !iter_res->gend_fragment(); iter_res->next_fragment())
	{
		res = (Residue *) iter_res->get_fragment();
		iter_atom = new pdbIter( res );

		for( iter_atom->pos_atom=0; !iter_atom->gend_atom(); iter_atom->next_atom() )
		{
	//		printf("atom name: -%s-\n",(iter_atom->get_atom())->getName() );
			if( strcmp( (iter_atom->get_atom())->getName(), " CA  ") == 0 )
			{
	//			printf("\tIt's a CA!\n");
				profile[res_index] = sqrt( pow(evec[nm_index+atom_index],2) + pow(evec[nm_index+atom_index+1],2) + pow(evec[nm_index+atom_index+2],2) );
			}
			atom_index+=3;
		}
		delete iter_atom;
		res_index++;
	}
	delete iter_res;
}

// Computes "sigma" from a vector (array)
float sig_vector(float *v, int size)
{
	double avg=0,sig=0;
	for( int i = 0; i < size; i++ )
		avg += v[i];
	avg /= size;
	for( int i = 0; i < size; i++ )
		sig += pow( v[i] - avg, 2 );
	sig /= size;
	return((float)sig);
}


// Creates, whether necessary, an array with "num_atoms" elemets
// copied from an array with "num_res) elements (the same values will be
// set for the same residue atoms)
// "p_out" = NULL --> will lead to memory allocation
void array_res2atom(Macromolecule *mol, double *array, double **p_out)
{
   int natoms,res_index=0,atom_index=0;
   Residue * res;
   pdbIter *iter_res,*iter_atom;
   double *out;

   natoms = mol->get_num_atoms();

   if(*p_out == NULL)
   {
   		if( !(out=(double *)malloc(sizeof(double)*natoms)) )
   		{
   			fprintf(stderr,"Msg(array_res2atom): Unable to allocate memory!\n");
   			exit(1);
   		}
   }
   else
   		out = *p_out; // avoids segmentation fault

   iter_res = new pdbIter( mol );
   for( iter_res->pos_fragment = 0; !iter_res->gend_fragment();	iter_res->next_fragment() )
   {
   		// get residues
   		res = (Residue *) iter_res->get_fragment();
   		iter_atom = new pdbIter( res );
	    for( iter_atom->pos_atom = 0; !iter_atom->gend_atom(); iter_atom->next_atom() )
	    {
	    	out[atom_index] = array[res_index]; // loads "atomic" array with the same order
//	    	printf("%f ",out[atom_index]);
	    	atom_index++;
	    }
	    delete iter_atom;
	    res_index++;
   }
//   printf("\n");
   delete iter_res;
   *p_out = out; // outputs array
}

// Normalizes lineally in order to match sigma to "max"
void norm_vol(vlVolume *vol,float max,float sig)
{
	float sigma;
	if(sig<=0)
		sigma = FOPS::sigma(vol); // computes sigma
	else
		sigma=sig; // normalizes with previously computed sigma
	FOPS::mul(vol,max/sigma);
}

// Note: To see "DFFE score functions" made since 24/7/2008, please
// take a look to "ADP_Workspace_11_01_2009.tgz"

// Mon made (6/10/2008)
// Difference map Energy
// map = model map (already allocated)
// target = target map (already allocated)
// offsets = array with the voxel offsets with any chance to be >0 (voxels inside target map)
// n_max = number of elements in "offsets" array (# of voxels that could be >0)
// ("diff" dffe_method)
float dffe_diff2(vlVolume *map, vlVolume *target, int *offsets, int n_max)
{
    float voxelD;
    double energy=0.0;
//	int n_vox;
//  vlDim dim = diffmap->dim();
//	n_vox = dim.x() * dim.y() * dim.z(); // total voxels (normalization)
    float *container_map = (float*)map->getVoxelVoidPtr(vlPoint3ui(0,0,0));
    float *container_target = (float*)target->getVoxelVoidPtr(vlPoint3ui(0,0,0));

	// Adding positive mass
	for( int i = 0; i < n_max; i++ )
	{
	    voxelD = *( container_target + offsets[i] ); // target map voxel
	    voxelD -= *( container_map + offsets[i] ); // model map voxel
        if( voxelD > 0.0 )
	       	energy += voxelD;
	}

//    return energy/n_vox;
    return (float) energy;
}

// Mon made (25/11/2008)
// UNDER CONSTRUCTION
// It checks both maps maximum dimensions, and makes a new map to where both maps fit.
vlVolume *fitsize(vlVolume *map, vlVolume *map2)
{
	int dx,dy,dz,dx2,dy2,dz2;
	vlPoint3f orig,orig2;
	vlPoint3f oldPosition, newPosition;
	vlVolume *out=NULL;

	// Get sizes and origins
	dx = map->dim().x();
	dy = map->dim().y();
	dz = map->dim().z();
	dx2 = map2->dim().x();
	dy2 = map2->dim().y();
	dz2 = map2->dim().z();
    map->getPosition(&orig);
    map2->getPosition(&orig2);

//    newPosition.x(oldPosition.x()-(float)(margin.x()* old->units().x() ));
//    newPosition.y(oldPosition.y()-(float)(margin.y()* old->units().y() ));
//    newPosition.z(oldPosition.z()-(float)(margin.z()* old->units().z() ));
//    padded->setPosition(newPosition);

  return out;
}


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
  Macromolecule *(*fptr)(Macromolecule *), // (OPTIONAL)
	// Then, "mol" will be used for moveing and "mol2" & "fptr(mol)" for scoring.
	// fptr = Macromolecule *(*fptr)(Macromolecule *)  (this is mandatory in: CA-model)
  Macromolecule **p_mol2 // (OPTIONAL) p_mol2 --> Initial Scoring Macromolecule
	// (Note that "mol2" will be a reference to an atoms sub-set in "mol")
#ifdef USE_PTHREAD // Enables PThread parallel routines
  , int nthreads // Number of threads used in parallelization
#endif
  )
{
	int verbose=0;
	int sel_chain=-1;
	int refine_iter = 4; // Number of iterations per Rot.&Tras.-lational refinement
	double deltas[6],result[6]; //
	double corr,corr_new,corr_iter; // actual correlation
	double fx,fa,fb;
	Macromolecule *mol_iter,*mol_iter2,*mol_delta,*mol_delta2,*buff; // new dummys
	float shake;
	bool apply_pose;

	corr_iter = corr = corr0;
	if( verbose > 0 )
		printf("pose_ref> Correlation before any refinement: %f\n",corr);

	// COPY INPUT MACROMOLECULES (because motion may not be successful)
	if(fptr!=NULL) // CA-model
	{
		mol_iter2 = new Macromolecule(*p_mol); // loading PDB copy (from origin)
		mol_iter = fptr(mol_iter2); // just scoring // NEEDS TO BE INSIDE
		mol_delta2 = new Macromolecule(*p_mol); // loading PDB copy (from origin)
		mol_delta = fptr(mol_delta2); // just scoring // NEEDS TO BE INSIDE
	}
	else
	{
		mol_iter = new Macromolecule(*p_mol); // loading PDB copy (from origin)
		mol_delta = new Macromolecule(*p_mol); // loading PDB copy (from origin)
	}

	int index_chain=0;
	// Multiple chains loop
	do {
		for(int n_ref=0; n_ref<refine_iter; n_ref++)
		{
			// Chain choice
			sel_chain = chain_order[index_chain]; // "rotrans" will act only over the seleted chain

			if( verbose > 1 )
				printf("pose_ref> Correlation before refinement iter %d: %f (chain= %d)\n",n_ref+1,corr_iter,sel_chain);

//			if(poseRand_switch)
//			{
//				for(int i=0; i<3; i++)
//					result[i] = ( ( (double) rg->Random() * 2 ) - 1 ) * shake_mov; // (-1:1) Playing dice with Mersenne!
//				for(int i=3; i<6; i++)
//					result[i] = ( ( (double) rg->Random() * 2 ) - 1 ) * shake_rot; // (-1:1) Playing dice with Mersenne!
//				n_ref = refine_iter; // Forces loop exit
//			}
//			else
				// Numerical Recipes, pag 402.
				// 10.2 Parabolic Interpolation (1D x 6)
				// INVERSE PARABOLIC INTERPOLATION
				fx = 1 - corr_iter;
				for(int j=0; j<6; j++)
					result[j] = 0.0; // reset

				// Iters each dimension (6D)
				for(int i=0; i<6; i++)
				{
					if(fptr!=NULL) // CA-model
					{
//						mol_delta2 = new Macromolecule(mol_iter2); // loading PDB copy (from origin)
//						mol_delta = fptr(mol_delta2); // just scoring // NEEDS TO BE INSIDE
						mol_delta2->copy_coordinates(mol_iter2); // copy current coords.
					}
					else
					{
//						mol_delta = new Macromolecule(mol_iter); // loading PDB copy (from origin)
						mol_delta->copy_coordinates(mol_iter); // copy current coords.
					}

					for(int j=0; j<6; j++)
						deltas[j] = 0.0; // reset

					// MOVEMENT "+"
					if( i < 3 ) // TRASLATION
					{
						// "shake" decreases as "refine_iter" increases (to force convergence)
						shake = shake_mov/powf(2.0,(float)n_ref);
						deltas[i] = shake;
					}
					else // ROTATION
					{
						shake = shake_rot/powf(2.0,(float)n_ref);
						deltas[i] = shake;
					}

					// Rotates and traslates a Macromolecule, checking whether its necessary!
					if(fptr!=NULL) // CA-model
						rotrans(mol_delta2, (float) deltas[0], (float) deltas[1], (float) deltas[2],
								(float) deltas[3], (float) deltas[4], (float) deltas[5], sel_chain );
					else
						rotrans(mol_delta, (float) deltas[0], (float) deltas[1], (float) deltas[2],
								(float) deltas[3], (float) deltas[4], (float) deltas[5], sel_chain );

#ifdef USE_PTHREAD // Enables PThread parallel routines
					// PROJECTION AND FILTRATION
					project_pdb(mol_delta,p_volF,p_vol_dummy,kernel,dim_vox,filter,resolution,false,nthreads);
#else
					project_pdb(mol_delta,p_volF,p_vol_dummy,kernel,dim_vox,filter,resolution,false);
#endif

					// If there were some negative values, cross-corr may be negative...
					FOPS::threshold( *p_volF, model_thr );

					/**
					* Mon made (3/6/2008)
					* Cross.Corr. inside a cubic-frame (nozero=false) or inside the final volume (nozero=true, default)
					* size = kernel width
					* vol = Target (Final) volume
					* If we have avg1 and sig1 (from: vol), we should provide them!
					*/
					// float correlation_frame(vlVolume *vol2, vlVolume *vol, vlDim size, bool nozero, bool volcalc, float volavg, float volsig);

					// Evaluating Cross.corr.
					fb = 1 - FOPS::correlation_frame(*p_volF,*p_volR,pad,corr_nozero,false,avgR,sigR); // >0 movement

					// MOVEMENT "-"
					deltas[i] *= -2; // because we must undo the +++ movement !!!

					// Rotates and traslates a Macromolecule, checking whether its necessary!
					if(fptr!=NULL) // CA-model
						rotrans(mol_delta2, (float) deltas[0], (float) deltas[1], (float) deltas[2],
								(float) deltas[3], (float) deltas[4], (float) deltas[5], sel_chain );
					else
						rotrans(mol_delta, (float) deltas[0], (float) deltas[1], (float) deltas[2],
								(float) deltas[3], (float) deltas[4], (float) deltas[5], sel_chain );

#ifdef USE_PTHREAD // Enables PThread parallel routines
					// PROJECTION AND FILTRATION
					project_pdb(mol_delta,p_volF,p_vol_dummy,kernel,dim_vox,filter,resolution,false,nthreads);
#else
					project_pdb(mol_delta,p_volF,p_vol_dummy,kernel,dim_vox,filter,resolution,false);
#endif

					// If there were some negative values, cross-corr may be negative...
					FOPS::threshold( *p_volF, model_thr );

					// Evaluating Cross.corr.
					fa = 1 - FOPS::correlation_frame(*p_volF,*p_volR,pad,corr_nozero,false,avgR,sigR); // <0 movement

					// INVERSE PARABOLIC INTERPOLATION
					result[i] = parab_iterpol(fx, fa, fb, shake);

					if( verbose > 1 )
						printf("pose_ref> (D%d) (fx=%5.4f fa=%5.4f fb=%5.4f) (shake=%9.6f) = %9.6f\n",
								i+1,fx,fa,fb,shake,result[i]);

//					if(fptr!=NULL) // CA-model
//					{
//						mol_delta->erase_level(pdb_atom); // remove CA selection
//						delete mol_delta2;
//					}
//					delete mol_delta;
				}

			if( verbose > 1 )
				printf( "pose_ref> REF. (iter= %d): %5.3f %5.3f %5.3f (Rot,º) %5.3f %5.3f %5.3f (Trans,A)\n"
						, n_ref+1, result[3]*(180/M_PI), result[4]*(180/M_PI), result[5]*(180/M_PI)
						, result[0], result[1], result[2] );

			// Whether all movement componets are smaller than tolerance (refine_tol)
			if( abs(result[0]) < refine_tol
					&& abs(result[1]) < refine_tol
					&& abs(result[2]) < refine_tol
					&& abs(result[3])*(180/M_PI) < refine_tol
					&& abs(result[4])*(180/M_PI) < refine_tol
					&& abs(result[5])*(180/M_PI) < refine_tol )
			{
				if( verbose > 1 )
					printf("pose_ref> SMALL MOVEMENT! (iter= %d)--> DON'T MOVING!\n",n_ref);
				n_ref= refine_iter; // forces loop exit!
			}
			else // Significant movement! Accepting it!
			{
				// Computing new correlation
				if(fptr!=NULL) // CA-model
					rotrans(mol_iter2, (float) result[0], (float) result[1], (float) result[2],
							(float) result[3], (float) result[4], (float) result[5], sel_chain );
				else
					rotrans(mol_iter, (float) result[0], (float) result[1], (float) result[2],
							(float) result[3], (float) result[4], (float) result[5], sel_chain );

#ifdef USE_PTHREAD // Enables PThread parallel routines
				project_pdb(mol_iter,p_volF,p_vol_dummy,kernel,dim_vox,filter,resolution,false,nthreads);
#else
				project_pdb(mol_iter,p_volF,p_vol_dummy,kernel,dim_vox,filter,resolution,false);
#endif

				// If there were some negative values, cross-corr may be negative...
				FOPS::threshold( *p_volF, model_thr );

				corr_iter = FOPS::correlation_frame(*p_volF,*p_volR,pad,corr_nozero,false,avgR,sigR);
				if( verbose > 1 )
					printf("pose_ref> Correlation after refinement iter %d: %f\n",n_ref+1,corr_iter);

//				sprintf(dummy_string,"%s_refpose%d%d.pdb",name,index_chain,n_ref);
//				mol_iter2->writePDB( dummy_string );

			}
		}
		index_chain++;
	} while( index_chain < n_chain && sel_chain >= 0 ); // iters chains

	if(fptr!=NULL) // CA-model
	{
		mol_delta->erase_level(pdb_atom); // remove CA selection
		delete mol_delta2;
	}
	delete mol_delta;

	for(int i=0;i<n_chain;i++)
		chain_order[i] = -1; // initialization


	if(corr_iter > corr) // Allways accepted
	{
		apply_pose = true;
		if(verbose > 0)
			printf( "pose_ref> POSE REFINEMENT: IMPROVED! old=%f new=%f UPDATED!\n",corr,corr_iter );
	}
	else // Monte Carlo
	{
		apply_pose = false;
		if(verbose > 0)
			printf( "pose_ref> POSE REFINEMENT: NOT-IMPROVED! old=%f new=%f AND REJECTED!\n",corr,corr_iter);
	}

	if(apply_pose) // Standard Cross-correlation check
	{
		// MOTION ACCEPTED!
		if(fptr!=NULL) // CA-model
		{
			(*p_mol2)->erase_level(pdb_atom); // remove CA selection
			delete *p_mol2;
			delete *p_mol;
			*p_mol = mol_iter2; // moving macromolecule
			*p_mol2 = mol_iter; // scoring macromolecule
		}
		else
		{
			delete *p_mol;
			*p_mol = mol_iter; // Macromolecule moved and updated
		}

		return corr_iter; // correlation improved
	}
	else
	{
		if(fptr!=NULL) // CA-model
		{
			mol_iter->erase_level(pdb_atom);
			delete mol_iter2;
		}
		delete mol_iter;
		if( verbose > 0 )
			printf( "pose_ref> POSE REFINEMENT (chain %d): Corr. NOT IMPROVED! old=%f new=%f DISCARDED\n",sel_chain,corr0,corr_iter );

		return -1.0; // correlation didn't improved
	}
}

//// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
//// allocating profile memory if (*p_profile==NULL).
//// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
//// Damma & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
//// Comparison for Flexible Proteins and Predicted Protein Structures".
//// Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
//void gaussian_weight(Macromolecule *mol, Macromolecule *mol2, double **p_profile, double c)
//{
//	int num_atoms;
//	num_atoms = mol->get_num_atoms();
//
//	// some initial checking
//	if(num_atoms != mol2->get_num_atoms())
//	{
//		fprintf(stderr,"Msg(gaussian_weight): Sorry, different number of atoms %d and %d\n",num_atoms,mol2->get_num_atoms());
//		exit(1);
//	}
//
//	double *profile;
//	if(*p_profile == NULL)
//	{
//		profile = (double *) malloc(sizeof(double) * num_atoms);
//		check_pointer(profile,"gaussian_weight weights profile allocation");
//		*p_profile = profile;
//	}
//	else
//		profile = *p_profile;
//
//	pdbIter *iter,*iter2;
//	Tcoor pos,pos2;
//	double d2;
//	iter = new pdbIter(mol);
//	iter2 = new pdbIter(mol2);
////	fprintf(stderr,"profile:\n");
//	for( iter->pos_atom=0,iter2->pos_atom=0; !iter->gend_atom(); iter->next_atom(),iter2->next_atom() )
//	{
//		(iter->get_atom())->getPosition(pos);
//		(iter2->get_atom())->getPosition(pos2);
//		d2 = pow(pos[0]-pos2[0],2) + pow(pos[1]-pos2[1],2) + pow(pos[2]-pos2[2],2);
//		profile[iter->pos_atom] = exp( -d2/c );
////		fprintf(stderr,"%f ",profile[iter->pos_atom]);
//	}
////	fprintf(stderr,"\n");
//	delete iter;
//	delete iter2;
////	fprintf(stderr,"Msg(gaussian_weight): Success!\n");
//}
//
//// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
//// allocating profile memory if (*p_profile==NULL).
//// (17/7/2012) --> Taking into account sequence alignment (at atomic level) for different size molecules.
////                 *The Weights profile its relative to the first molecule (mol)
//// Weights are computed according to eq.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
//// Damma & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
////                   Comparison for Flexible Proteins and Predicted Protein Structures".
////                   Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
//void gaussian_weight(Macromolecule *mol, Macromolecule *mol2, double **p_profile, bool *mask1, bool *mask2, double c)
//{
//	bool debug = false;
//	int num_atoms;
//	num_atoms = mol->get_num_atoms();
//
////	// some initial checking
////	if(num_atoms != mol2->get_num_atoms())
////	{
////		fprintf(stderr,"Msg(gaussian_weight): Sorry, different number of atoms %d and %d\n",num_atoms,mol2->get_num_atoms());
////		exit(1);
////	}
//
//	double *profile;
//	if(*p_profile == NULL)
//	{
//		profile = (double *) malloc(sizeof(double) * num_atoms);
//		check_pointer(profile,"gaussian_weight weights profile allocation");
//		*p_profile = profile;
//	}
//	else
//		profile = *p_profile;
//
//	// Initialization (mandatory if masks are provided!!!)
//	for(int i=0; i<num_atoms; i++)
//		profile[i] = 0.0; // some negative number for marking atoms that should not be considered
//
//	pdbIter *iter1,*iter2;
//	Atom *at1,*at2;
//	Tcoor pos1,pos2;
//	double d2;
//	iter1 = new pdbIter(mol);
//	iter2 = new pdbIter(mol2);
//	int ncommon = 0;
////	fprintf(stderr,"profile:\n");
//
//
////	pdbIter *iter1,*iter2;
////	iter1 = new pdbIter(this,false);
////	iter2 = new pdbIter(mol2,false);
//
//	iter2->pos_atom=0;
//	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
//	{
//		if(mask1[iter1->pos_atom]) // Mol 1 matching atom
//		{
//			while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding atom
//				iter2->next_atom();
//
//			if(!iter2->gend_atom())
//			{
//				// get atoms
//				at1 = iter1->get_atom();
//				at2 = iter2->get_atom();
//
//				// get positions
//				at1->getPosition(pos1);
//				at2->getPosition(pos2);
//
//				if(debug)
//					printf("Getting atoms %d (%s) and %d (%s)\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName());
//
//				// Computing Weighted RMSD profile
//				d2 = pow(pos1[0]-pos2[0],2) + pow(pos1[1]-pos2[1],2) + pow(pos1[2]-pos2[2],2);
//				profile[iter1->pos_atom] = exp( -d2/c );
////				profile[ncommon] = exp( -d2/c );
//
//				if(debug)
//					printf("profile[%d]= %f\n",iter1->pos_atom,profile[iter1->pos_atom]);
////					printf("profile[%d]= %f\n",iter1->pos_atom,profile[ncommon]);
//
//				ncommon++;
//				iter2->next_atom();
//			}
//		}
//	}
//
////	// Some checking
////	if(num_atoms != iter1->pos_atom)
//
//
//	delete iter1;
//	delete iter2;
////	fprintf(stderr,"Msg(gaussian_weight): Success!\n");
//}
