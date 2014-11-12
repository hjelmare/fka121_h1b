/*
 MD_main.c
 
 Created by AL on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#include "MD_functions.h"

/* Main program */
int main()
{
	// REMEMBER: \m/ METAL UNITS \m/
	// simulation settings
	double totalTime = 1;
	double timestep = 0.1;
		
	// physical parameters
	int dim = 3;		// there's a few functions that are hardcoded for dim=3, not sure how to fix...
	int nCells = 4;
	int nParticles = 4*pow(nCells,dim);
	double mass = 0.00279636;  // 26.9815 u
	double latticeParameter = 4.05;
	double maxDeviation = 0.05;

	// storage of physical quantities
	double pos[nParticles][dim];
	double vel[nParticles][dim];
	double force[nParticles][dim];
	double energy, potentialEnergy, kineticEnergy;

	// derived quantities
	int nSteps = (int) totalTime/timestep;
	double supercellLength = nCells * latticeParameter;

	// other stuff
	int i,j,k;   

	srand(time(NULL));
	double random_value;

	init_fcc(pos, 4, latticeParameter);

	for(i = 0; i<nParticles; i++){
		for (j = 0; j<3; j++){
			random_value = (double) rand() / (double) RAND_MAX;
			pos[i][j] = pos[i][j] + random_value * maxDeviation * latticeParameter;
 		}
	}

	// set velocities to zero
	for (i = 0; i<nParticles; i++){
		for (j = 0; j<3; j++){
			vel[i][j] = 0;
		}
	}
	
	potentialEnergy = get_energy_AL(pos, supercellLength, nParticles);
	kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
	energy = potentialEnergy + kineticEnergy;

	printf("energy = %e \n", potentialEnergy);

	// so, what needs to be done? (for task 1, to begin with...)
	// a main loop which
	// -gets forces (there's a function for that)
	// -moves things (verlet? or something else? did i miss something?)
	// -calculates pe, ke, and E, and saves them for matlab plotting

	//Saving initial data:
	FILE *kineticEnergyFile;
	kineticEnergyFile = fopen("kineticEnergy.data","w");
	fprintf(kineticEnergyFile, "%e \t %e \n", 0.0, kineticEnergy );

	FILE *potentialEnergyFile;
	potentialEnergyFile = fopen("potentialEnergy.data","w");
	fprintf(potentialEnergyFile, "%e \t %e \n", 0.0, potentialEnergy);
	
	FILE *totEnergyFile;
	totEnergyFile = fopen("totEnergy.data","w");
	fprintf(totEnergyFile, "%e \t %e \n", 0.0, energy);


	
	for (i=1;i<nSteps;i++) {
		// Update velocities and positions
		for (j=0; j<nParticles; j++) {
			for (k = 0; k<dim; k++) {
				vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
				pos[j][k] = pos[j][k] + timestep*vel[j][k];
			}
		}
		
		get_forces_AL(force, pos, supercellLength, nParticles);

		// Update velocities again
		for (j=0; j<nParticles; j++) {
			for (k=0; k<dim; k++) {
				vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
			}
		}

		potentialEnergy = get_energy_AL(pos, supercellLength, nParticles);
		kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
		energy = potentialEnergy + kineticEnergy;
		printf("e %e, pe %e, ke %e, v11 %e \n", energy, potentialEnergy, kineticEnergy, vel[0][0]);

		fprintf(kineticEnergyFile, "%e \t %e \n", i*timestep, kineticEnergy );
		fprintf(potentialEnergyFile, "%e \t %e \n", i*timestep, potentialEnergy);
		fprintf(totEnergyFile, "%e \t %e \n", i*timestep, energy);
	}

    /*
     Descriptions of the different functions in the files initfcc.c and alpotential.c are listed below.
     */
    
    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of primitive cells in each direction and a0 is the lattice parameter. The positions of all the atoms are stored in pos which should be a matrix of the size N x 3, where N is the number of atoms. The first, second and third column correspond to the x,y and z coordinate respectively.
     */
    /*
     init_fcc(pos, Nc, a0);
    */
    
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be a matrix containing the positions of all the atoms, L is the length of the supercell and N is the number of atoms.
     */
    /*
     double energy;
     energy = get_energy_AL(pos, L, N);
     */
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix containing the positions of all the atoms, L is the length of the supercell and N is the number of atoms.
     */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
     */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the forces are stored in f which should be a matrix of size N x 3, where N is the number of atoms and column 1,2 and 3 correspond to the x,y and z component of the force resepctively . pos should be a matrix containing the positions of all the atoms, L is the length of the supercell and N is the number of atoms.
     */
    /*
     get_forces_AL(f,pos, L, N);
     */
}
