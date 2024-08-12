/** ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
 main.cpp  -  description ------------------- begin                : Wed Sep 22 18:08:48 CEST 2004
 copyright            : (C) 2004 by Jose Ignacio Garzon email                :
 ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff */

/** ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

This program is free software; you can redistribute it and/or modify  *
 it under the terms of the GNU General Public License as published by  *
 the Free Software Foundation; either version 2 of the License, or     *
 (at your option) any later version.                                   *

fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff */

#include <stdio.h> // needed by some compilers
#include "cmdl/CmdLine.h"
//#include "libnma/nma.h"
#include "libnma/include/libnma_io.h"
#include "libnma/include/libnma_diag.h"
#include "libnma/include/libnma_def.h"
//#include <iostream>
//#include <libpdb/world.h>
//#include <libtools/timer.h>
//#include <libtools/Matrix4x4.h>
//#include "libpdb/pdbIter.h"
//#include "libtools/timer.h"
//#include "libpdb/ResIni.h"
//#include "nma.h"
//#include "def.h"
//#include <string>
//#include <stdlib.h>

// input variables
using namespace TCLAP;
void parseOptions(int argc, char** argv);
float power;
char file_pdb[28]; // pdb input
int method;
float cutoff;
bool deform=true;
bool no_mass=false;
bool saveno=true;
bool savea=false;

float cte=40;

void parseOptions(int argc, char** argv)
{
        std::string temp;
//        CmdLine cmd("nmac","   ", "0.99" );
        CmdLine cmd("defprot","   ", "1.00" );

        try {

        //
        // Define required arguments no labeled
        //
        UnlabeledValueArg<std::string> pdb("pdb","pdb input file","default","pdb");
        cmd.add( pdb );
        //
        // Define labeled arguments
        //


        SwitchArg def("d","deform", "Turn off deformability calculations", true);
        cmd.add( def );

        SwitchArg nomass("n","nomass", "no mass weighting (mass=1)", true);
        cmd.add( nomass );

        SwitchArg Hin("m", "hinsen","Hinsen distance criterion",  true);
        cmd.add( Hin );

        SwitchArg Savea("", "ascii","Save modes ascii",  false);
        cmd.add( Savea );

        SwitchArg Saveno("", "nosave","Skip saving the modes",  true);
        cmd.add( Saveno );

        ValueArg<float> Cte("k", "cons","Strenght constant (default 40.0)",  false,40.0,"float");
        cmd.add( Cte );

        ValueArg<float>  Power("p","power", "conectivity (1/d)", false,6.0,"float");
        cmd.add( Power);

        ValueArg<float> Cutoff("c","cutoff"," distance threshold", false, 8,"float");
        cmd.add( Cutoff );


        // Parse the command line.
        cmd.parse(argc,argv);

        strcpy(file_pdb,((temp=pdb.getValue()).c_str()));

        // esto  no debia estar aqui...
        FILE *f;
        if ( (f=fopen(file_pdb, "r"))==NULL) {
        fprintf(stderr, "\n  Error->Cannot open file '%s'\n\n", file_pdb);
        exit(1);
        } fclose(f);


        cutoff = Cutoff.getValue();
        power = Power.getValue();
        cte = Cte.getValue();
        method=3;
        if (Power.isSet())  method=1;
        if (Hin.isSet())  method=2;
        if (Cutoff.isSet())  method=0;

        if (def.isSet())  deform = def.getValue();

        if (nomass.isSet())  no_mass = true;

        if (Savea.isSet())  savea=true;
        if (Saveno.isSet())  saveno=false;


        } catch ( ArgException& e )
        { std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl; }

        //cmd.~CmdLine();

}




int main( int argc, char * argv[] )
{
	// Initialize aminoacids and nucleotids
	init_aminoacids();

  parseOptions(argc,argv);

  // reading pdb
  Macromolecule * molr = new Macromolecule( "pdb" );
  molr->readPDB( file_pdb );
  int num_atoms=0;
  num_atoms = molr->get_num_atoms();
  fprintf(stderr, "nmac>\nnmac> Atoms read -> %d\n", num_atoms );
  if (num_atoms==0)  {
  fprintf(stderr, "nmac>\nnmac> Error %s seems to be a not valid PDB\n", file_pdb );
  exit(0);
  }
  //molr->info();

  // CA  Conditions
  Condition * calpha = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
  calpha->add( " CA " );
  Conditions * calpha2 = new Conditions();
  calpha2->add( calpha );

  // if CA selection
  Macromolecule * mol = molr->select( calpha2 );
  num_atoms = mol->get_num_atoms();
  fprintf(stderr, "nmac> CA atoms detected -> %d\n", num_atoms );
  if (num_atoms>5000)  {
   fprintf(stderr, "nmac>\nnmac> Error many CA atoms (5000 is the limit) \n");
    exit(0);
  }

  mol->writePDB( "ca.pdb");
  fprintf(stderr,"nmac> saved CA only PDB file->ca.pdb\n");

  // compute pair distance matrix
  double * dist_matrix;
  mol->distanceMatrix( &dist_matrix );

  // check min and max pair distances
  int i, j, index;
  double maxdist, mindist, currdist, currmindist;


  mindist = 1e20;
  maxdist = -1;

  /* check min and max pairwise distances */
  mindist = 1e20;
  maxdist = -1;
  for ( i = 0; i < num_atoms; i++ )
  {
    currmindist = 1e20;
    for ( j = 0; j < num_atoms; j++ )
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

  // input parameters

/* compute background contact matrix */
float * cont_matrix;
cont_matrix = ( float * ) malloc( num_atoms * ( num_atoms - 1 ) / 2 * sizeof( float ) );
for ( i = 0; i <  num_atoms * ( num_atoms - 1 ) / 2 ; i++ ) cont_matrix[i] = 0.0;



switch(method){

 case 0:
//
// DISTANCE CUTOFF METHOD
//

   fprintf(stderr, "nmac>\nnmac> Cutoff Distance method\nnmac>\n");
   fprintf(stderr, "nmac> Range of closest neighbor distances: %5.2f - %5.2f.\n", mindist, maxdist );
   fprintf(stderr, "nmac> Distance cutoff %f\n", cutoff );
   if ( cutoff < maxdist )
   {
     fprintf(stderr, "nmac> Cutoff must be at least %5.2f to ensure that all atoms have at least one connected neighbor\n\n", maxdist );
     exit( 1 );
   }

   for ( i = 0; i < num_atoms * ( num_atoms - 1 ) / 2; i++ )
    if ( dist_matrix[i] <= cutoff)  cont_matrix[i] = cte;
   break;


 case 1:
//
// POWER DISTANCE METHOD
//

/* input power of distance for contact matrix */
fprintf(stderr, "nmac>\nnmac> Power distance method\nnmac>\n");
fprintf(stderr, "nmac> Distance power: %f\n", power );

for ( i = 0; i < num_atoms * ( num_atoms - 1 ) / 2; i++ )
    if (dist_matrix[i]<=4.0)
     cont_matrix[i] = cte / ( 1.0 + pow( dist_matrix[i] / 3.8, (double)power ) );
      else
       cont_matrix[i] = cte / ( 1.0 + pow( dist_matrix[i] / 3.8, (double)power ) );
  break;


case 2:
//
// HINSEN DISTANCE METHOD
//

/* input power of distance for contact matrix */
fprintf(stderr, "nmac> Following Hinsen's distance criterion\n");

for ( i = 0; i < num_atoms * ( num_atoms - 1 ) / 2; i++ )
   if (dist_matrix[i]<=4.0)
    cont_matrix[i] = 86000*dist_matrix[i]-23900;
    else
    cont_matrix[i] = 128/ (pow( dist_matrix[i]/10, 6 ) );
break;

default:

/* input power of distance for contact matrix */
power=6.0;
fprintf(stderr, "nmac>\nnmac> Power distance method\nnmac>\n");
fprintf(stderr, "nmac> Distance power: %f\n", power );

for ( i = 0; i < num_atoms * ( num_atoms - 1 ) / 2; i++ )
    if (dist_matrix[i]<=4.0)
     cont_matrix[i] = cte / ( 1.0 + pow( dist_matrix[i] / 3.8, (double)power ) );
      else
       cont_matrix[i] = cte / ( 1.0 + pow( dist_matrix[i] / 3.8, (double)power ) );
  break;

 }



  /* getting coordinates single row */
  float * coord;
  mol->coordMatrix( & coord );

  fprintf(stderr, "nmac>\nnmac> Computing Hessian Matrix\n");
  /* compute Hessian matrix */
  float * hess_matrix;
  int size = num_atoms * 3;
  hess_matrix = ( float * ) malloc( size * size * sizeof( float ) );
  for ( i = 0; i < size*size; i++ ) hess_matrix[i] = 0.0;


  double r[3], rsqu;
  int m, n;

/* off-diagonal */
for ( i = 0; i < num_atoms; ++i )
  for ( j = i+1; j < num_atoms; ++j )
  {
    r[0] = coord[i * 3] - coord[j * 3];
    r[1] = coord[i * 3 + 1] - coord[j * 3 + 1];
    r[2] = coord[i * 3 + 2] - coord[j * 3 + 2];
    index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2;
    rsqu = dist_matrix[index] * dist_matrix[index];
    for ( m = 0; m < 3; ++m )
      for ( n = 0; n < 3; ++n )
      {
        hess_matrix[( 3 * i + m ) + size * ( 3 * j + n )] = -cont_matrix[index] * r[m] * r[n] / rsqu;
        hess_matrix[( 3 * j + m ) + size * ( 3 * i + n )] = -cont_matrix[index] * r[m] * r[n] / rsqu;
      }
  }


/* on-diagonal */
for ( i = 0; i < num_atoms; ++i )
  for ( j = i+1; j < num_atoms; ++j )
  {
    r[0] = coord[i * 3] - coord[j * 3];
    r[1] = coord[i * 3 + 1] - coord[j * 3 + 1];
    r[2] = coord[i * 3 + 2] - coord[j * 3 + 2];

    index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2;
    rsqu = dist_matrix[index] * dist_matrix[index];

    for ( m = 0; m < 3; ++m )
      for ( n = 0; n < 3; ++n )
      {
        hess_matrix[( 3 * i + m ) + size * ( 3 * i + n )] += cont_matrix[index] * r[m] * r[n] / rsqu;
        hess_matrix[( 3 * j + m ) + size * ( 3 * j + n )] += cont_matrix[index] * r[m] * r[n] / rsqu;
      }
  }



  // free( dist_matrix );
  free( cont_matrix );

  float * mass;
  mass = ( float * ) malloc( num_atoms * sizeof( float ) );
  pdbIter * iter2 = new pdbIter( mol );
  Residue * res;
   double dump=0;

  if (no_mass)
    for ( i = 0; i < num_atoms; i++ )
     mass[i] = 1;
    else
    {
//  init_aminoacids();


  i = 0;
  for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
  {
    res = ( Residue * ) iter2->get_fragment();
    j = resnum_from_resname( res->getName() );
    mass[i++] = AA[j].mass;
    dump+=mass[i-1];
//   printf("%d %s %f %f\n", i-1, AA[j].aa_name3, mass[i-1],dump);
  }


  for ( i = 0; i < num_atoms; i++ )
    mass[i] /= dump;

    }

   /* multiply hess_matrix by mass matrix, i.e., form F = M^(-1/2).H.M^(-1/2) */
  for ( i = 0; i < num_atoms; ++i )
    for ( j = 0; j < num_atoms; ++j )
      for ( m = 0; m < 3; ++m )
        for ( n = 0; n < 3; ++n ) {
          hess_matrix[( 3 * i + m ) + size * ( 3 * j + n )] /= sqrt( mass[i] * mass[j] );
        //  hess_matrix[( 3 * j + m ) + size * ( 3 * i + n )] /= sqrt( mass[i] * mass[j] );
        }


/* input */
int neigval, null_modes=0;
neigval = size;
float * eigval, * eigvect;


// hessian diagonalization
fprintf(stderr, "nmac>\nnmac> Computing Hessian Diagonalization\n");
null_modes = nma (size, neigval, hess_matrix, &eigval, &eigvect);
fprintf(stderr,"nmac> 1st not NULL eigenvector %d scalefactor %f\n",
null_modes+1, sqrt(eigval[null_modes]));

free(hess_matrix);

/* switch from mass-weighted coords to standard cartesian coords. */
 for(i=0; i<neigval; i++)
   for(j=0; j<num_atoms; j++)
     for(m=0; m<3; m++)
       eigvect[i*size+3*j+m] /= sqrt(mass[j]);
free( mass );

// Saving modes

if (saveno) {

if (savea) {
  save_ascii_modes (size, null_modes, neigval, eigval, eigvect);
  fprintf(stderr, "nmac>\nnmac> Saving  modes in ascii files\n");
}

save_ptraj_modes (size, null_modes, neigval, eigval, eigvect);
fprintf(stderr, "nmac>\nnmac> Saving  modes in ptraj file nmac_ptraj.evec\n");

}

if (deform) {
// computing deformability
double *deform, *mob;
double *bf;
char *list;
bf = (double *) malloc(num_atoms * sizeof(double));

mol-> get_Pdbfact(bf);


compute_def (num_atoms,  coord, dist_matrix, eigval,  eigvect, neigval, null_modes, &deform, &mob);
fprintf(stderr,"nmac>\nnmac> Saving deformability/mobility in defmob.tab\n");

// saving deformability/mobility
list=mol->getResInfo();
save_ascii_defmob (num_atoms, deform, mob, bf, list);

//saving mobility in PDB
norm_defmob (num_atoms, deform);
mol-> exchange_Pdbfact (deform);
printf("nmac> Saving def.pdb (CA only)\n");
mol->writePDB("def.pdb");

molr-> exchange_Pdbfact(deform);
printf("nmac> Saving def_all.pdb (full atoms)\n");
molr->writePDB("def_all.pdb");


//saving mobilitiy in PDB
printf("nmac> Saving mob.pdb (CA only)\n");
norm_defmob (num_atoms, mob);
mol-> exchange_Pdbfact (mob);
mol->writePDB("mob.pdb");

free(eigval); free(eigvect);
free( dist_matrix );
free(deform);
free(mob);


} else {

  free(eigval); free(eigvect);
  free( dist_matrix );

}


fprintf(stderr,"nmac>\nnmac> The End\nnmac>\nnmac>\n");


}
