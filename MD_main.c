// MD_main.c

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#include "MD_functions.h"

#define PI 3.14159265368979

/* Main program */
int main()
{
  // REMEMBER: \m/ METAL UNITS \m/
  // simulation settings
  double equilibrationTime = 25;
  double productionTime = 25;
  double timestep = 0.01;
  //timestep = 0.1;
  //timestep = 0.001;

  double msdAverageSteps = 50;
    
  int nSpectrumPoints = 1000;
  double spectrumInterval = PI;
  
  double maxCorrelationTime = 1.0;
  
  double timeConstantT = 0.02;
  double timeConstantP = 0.05;
  timeConstantT = 1;
  timeConstantP = timeConstantT;

  // physical parameters
  int dim = 3;    // do not change
  int nCells = 4;
  int nParticles = 4*pow(nCells,dim);
  double wantedTemperature = 500+273;   // Target temperature
  wantedTemperature = 700+273;
//  wantedTemperature = 900+273;
  double wantedTemperature2 = 900+273;
  double wantedPressure =  6.32420934* 0.0000001; // 1 atm in metal units
  double mass = 0.00279636;  // 26.9815 u * 1.0364 * 0.0001
  double latticeParameter = 4.05;
  double maxDeviation = 0.05;

  // derived quantities
  double nEquilibrationSteps = equilibrationTime/timestep;
  double nProductionSteps = productionTime/timestep;
  int maxCorrelationSteps = maxCorrelationTime/timestep;


  // storage of physical quantities
  double pos[nParticles][dim];
  double vel[nParticles][dim];
  double force[nParticles][dim];
  double savedValuesT[maxCorrelationSteps];
  double meanValuesT[maxCorrelationSteps];
  double savedValuesP[maxCorrelationSteps];
  double meanValuesP[maxCorrelationSteps];
  double energy, potentialEnergy, kineticEnergy;
  double alphaT = 1;
  double alphaP = 1;
  double currentTemperature;
  double currentPressure;
  double sqrtAlphaT;
  double curtAlphaP;
  double virial;
  double volume;
  double nAverageSteps;
  double meanTemperature, meanSquareTemperature;
  double meanTemperature_ik[maxCorrelationSteps-1];
  double meanPressure, meanSquarePressure;
  double meanPressure_ik[maxCorrelationSteps-1];
  double phiTemperature[maxCorrelationSteps];
  double phiPressure[maxCorrelationSteps];
  double meanTemperatureSquare, meanPressureSquare;
  double varTemperature, varPressure;
  double sTemperature, sPressure;
  double savedPos[maxCorrelationSteps][nParticles][dim];
  double savedVelocities[maxCorrelationSteps][nParticles][dim];
  double meanVelocityScalar[maxCorrelationSteps][nParticles];
  double meanVelocityAverage[maxCorrelationSteps];
  double spectrum[nSpectrumPoints];
  double diffusionCoefficient;
  double meanDistanceScalar[maxCorrelationSteps][nParticles];
  double meanDistanceAverage[maxCorrelationSteps];
  double msdDiffCoeff = 0;
  int MTemperature, MPressure;

  // other stuff
  int i,j,k,m;   

  srand(time(NULL));
  double random_value;

// End of variable declarations, start of actual code

  printf("Initializing");
  init_fcc(pos, nCells, latticeParameter);

  // Applying random perturbations from strict fcc positions
  for(i = 0; i<nParticles; i++){
    for (j = 0; j<3; j++){
      random_value = (double) rand() / (double) RAND_MAX;
      pos[i][j] = pos[i][j] + (random_value-0.5) * maxDeviation * \
      latticeParameter;
    }
  }

  // Initialize velocities and forces to zero
  for (i = 0; i<nParticles; i++){
    for (j = 0; j<3; j++){
      vel[i][j] = 0;
      force[i][j] = 0;
    }
  }
  
  // Calculating initial energies
  potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, nParticles);
  kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
  energy = potentialEnergy + kineticEnergy;

  // Saving initial energies on the first line in the file
  FILE *energyFile;
  energyFile = fopen("energy.data","w");
  fprintf(energyFile, "%e \t %e \t %e \t %e \n", 0.0, energy, \
  potentialEnergy, kineticEnergy);

// Initialization is done, moving on to equilibration

  FILE *ptFile;
  ptFile = fopen("pt700.data","w");

  printf("\t"); // Progress indicator stuff
  for (i = 0; i < ( nEquilibrationSteps > nProductionSteps  ? \
  nEquilibrationSteps:nProductionSteps )/100; i++) {
      printf(".");
  }
  printf("\tDone!\nEquilibration\t");
///*
  for (i=0;i<nEquilibrationSteps;i++) {   // Start of equilibration loop--------------------------------------------------------------------------------------------------------------
    // Update velocities and positions
    for (j=0; j<nParticles; j++) {
      for (k = 0; k<dim; k++) {
        vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
        pos[j][k] = pos[j][k] + timestep*vel[j][k];
      }
    }

    // Calculate forces so that we can take the next half step
    get_forces_AL(force, pos, nCells*latticeParameter, nParticles);

    // Update velocities again
    for (j=0; j<nParticles; j++) {
      for (k=0; k<dim; k++) {
        vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
      }
    }
    
    // Calculating energies
    potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, \
    nParticles);
    kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
    energy = potentialEnergy + kineticEnergy;
    fprintf(energyFile,"%e \t %e \t %e \t %e \n", (i+1)*timestep, energy, \
    potentialEnergy, kineticEnergy);    

    // Get temperature and pressure
    currentTemperature = GetInstantTemperature(vel, nParticles, mass, dim);

    volume = pow(nCells*latticeParameter, 3);
    virial = get_virial_AL(pos, nCells*latticeParameter, nParticles);
    currentPressure = GetPressure(currentTemperature, volume, virial, \
    nParticles);

    fprintf(ptFile, "%e \t %e \t %e \t %e \n", (i+1)*timestep, currentPressure, \
    currentTemperature, latticeParameter);

    // Calculate scaling parameters for equilibration
    alphaT = GetAlphaT(wantedTemperature2, currentTemperature, timestep, \
    timeConstantT);
    alphaP = GetAlphaP(wantedPressure, currentPressure, timestep, \
    timeConstantP);
    sqrtAlphaT = sqrt(alphaT);
    curtAlphaP = pow(alphaP, 1.0 / 3);

    // Rescale velocity and position
    latticeParameter = latticeParameter * curtAlphaP;
    for (j=0; j<nParticles; j++) {
      for (k=0; k<dim; k++) {
        vel[j][k] = vel[j][k] * sqrtAlphaT;
        pos[j][k] = pos[j][k] * curtAlphaP;
      }
    }
//-Extra get_forces here... not sure we should be doing it twice per step, but it works...    
    // Calculate forces so that we can take the next half step
    get_forces_AL(force, pos, nCells*latticeParameter, nParticles);

    if( i % 100 == 0) {   // Progress indicator
      printf("I");
      fflush(stdout);
    }
  } // End of equilibration loop
//*/
// End of equilibration, start of production--------------------------------------------------------------------------------------------------------------------------------------------
   // Start of equilibration loop 2--------------------------------------------------------------------------------------------------------------
  for (i=0;i<nEquilibrationSteps;i++) {    
    // Update velocities and positions
    for (j=0; j<nParticles; j++) {
      for (k = 0; k<dim; k++) {
        vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
        pos[j][k] = pos[j][k] + timestep*vel[j][k];
      }
    }

    // Calculate forces so that we can take the next half step
    get_forces_AL(force, pos, nCells*latticeParameter, nParticles);

    // Update velocities again
    for (j=0; j<nParticles; j++) {
      for (k=0; k<dim; k++) {
        vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
      }
    }
    
    // Calculating energies
    potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, \
    nParticles);
    kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
    energy = potentialEnergy + kineticEnergy;
    fprintf(energyFile,"%e \t %e \t %e \t %e \n", (i+1)*timestep, energy, \
    potentialEnergy, kineticEnergy);    

    // Get temperature and pressure
    currentTemperature = GetInstantTemperature(vel, nParticles, mass, dim);

    volume = pow(nCells*latticeParameter, 3);
    virial = get_virial_AL(pos, nCells*latticeParameter, nParticles);
    currentPressure = GetPressure(currentTemperature, volume, virial, \
    nParticles);

    fprintf(ptFile, "%e \t %e \t %e \t %e \n", (i+1)*timestep, currentPressure, \
    currentTemperature, latticeParameter);

    // Calculate scaling parameters for equilibration
    alphaT = GetAlphaT(wantedTemperature, currentTemperature, timestep, \
    timeConstantT);
    alphaP = GetAlphaP(wantedPressure, currentPressure, timestep, \
    timeConstantP);
    sqrtAlphaT = sqrt(alphaT);
    curtAlphaP = pow(alphaP, 1.0 / 3);

    // Rescale velocity and position
    latticeParameter = latticeParameter * curtAlphaP;
    for (j=0; j<nParticles; j++) {
      for (k=0; k<dim; k++) {
        vel[j][k] = vel[j][k] * sqrtAlphaT;
        pos[j][k] = pos[j][k] * curtAlphaP;
      }
    }
//-Extra get_forces here... not sure we should be doing it twice per step, but it works...    
    // Calculate forces so that we can take the next half step
    get_forces_AL(force, pos, nCells*latticeParameter, nParticles);

    if( i % 100 == 0) {   // Progress indicator
      printf("I");
      fflush(stdout);
    }
  } // End of equilibration loop

// End of equilibration, start of production--------------------------------------------------------------------------------------------------------------------------------------------

  FILE *positionFile;
  positionFile = fopen("position.data","w");
  printf("\tDone!\nProduction\t");
  

  for (i=0;i<nProductionSteps;i++) {
    // Save positions of one particle to check whether the
    // aluminium is solid or liquid
    fprintf(positionFile, "%d \t %e \t %e \t %e \n", i, \
    pos[0][0], pos[0][1], pos[0][2]);

    // Update velocities and positions
    for (j=0; j<nParticles; j++) {
      for (k = 0; k<dim; k++) {
        vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
        pos[j][k] = pos[j][k] + timestep*vel[j][k];
      }
    }
    
    // Calculate forces so that we can take the next half step
    get_forces_AL(force, pos, nCells*latticeParameter, nParticles);

    // Update velocities again
    for (j=0; j<nParticles; j++) {
      for (k=0; k<dim; k++) {
        vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
      }
    }

    potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, \
    nParticles);
    kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
    energy = potentialEnergy + kineticEnergy;

    currentTemperature = GetInstantTemperature(vel, nParticles, mass, dim);
    volume = pow(nCells*latticeParameter, 3);
    virial = get_virial_AL(pos, nCells*latticeParameter, nParticles);
    currentPressure = GetPressure(currentTemperature, volume, virial, \
    nParticles);

    //Saves temp and pressure values in order to calculate s.
    if (i < maxCorrelationSteps) {
      savedValuesT[i] = currentTemperature;
      savedValuesP[i] = currentPressure;    
    } else {  
      // updates the saved values.
      for (j = 0; j < maxCorrelationSteps - 1; j++){  
        savedValuesT[j] = savedValuesT[j+1];
        savedValuesP[j] = savedValuesP[j+1];  
      }
      savedValuesT[maxCorrelationSteps - 1] = currentTemperature;
      savedValuesP[maxCorrelationSteps - 1] = currentPressure;
      
      // calculates and saves f_i*f_k
      for (j = 0; j<maxCorrelationSteps; j++){
        meanTemperature_ik[j] += savedValuesT[maxCorrelationSteps - 1] \
        * savedValuesT[maxCorrelationSteps - 1 - j];
        meanPressure_ik[j] += savedValuesP[maxCorrelationSteps - 1] \
        * savedValuesP[maxCorrelationSteps - 1 - j];
      }
      meanTemperature += currentTemperature;
      meanPressure += currentPressure;
      meanSquareTemperature += currentTemperature * currentTemperature;
      meanSquarePressure += currentPressure * currentPressure;
    }

    // Saving data for the velocity correlation function
    if (i < maxCorrelationSteps) {
      for (m = 0; m < nParticles ; m++) {
        savedVelocities[i][m][0] = vel[m][0];
        savedVelocities[i][m][1] = vel[m][1];
        savedVelocities[i][m][2] = vel[m][2];
      }
    } else {
      for (j = 0; j< maxCorrelationSteps- 1; j++) {
        for (m = 0; m<nParticles;  m++) {
          savedVelocities[j][m][0] = savedVelocities[j+1][m][0];
          savedVelocities[j][m][1] = savedVelocities[j+1][m][1];
          savedVelocities[j][m][2] = savedVelocities[j+1][m][2];
        }
      }
      for (m = 0; m<nParticles; m++) {
        savedVelocities[maxCorrelationSteps - 1][m][0] = vel[m][0];
        savedVelocities[maxCorrelationSteps - 1][m][1] = vel[m][1];
        savedVelocities[maxCorrelationSteps - 1][m][2] = vel[m][2];
      }

      for (j = 0; j<maxCorrelationSteps; j++) {
        for (m = 0; m < nParticles; m++) {
          meanVelocityScalar[j][m] += ScalarProduct( \
          savedVelocities[0][m], savedVelocities[j][m]);
        }
      }
    }

    if(i <  maxCorrelationSteps) {
      for(j = 0; j<nParticles; j++) {
        for (k = 0; k<dim; k++) {
          savedPos[i][j][k] = pos[j][k];
        } 
      }
    }else {
      for (j = 0; j<maxCorrelationSteps-1; j++) {
        for (k = 0; k<nParticles; k++) {
          for (m = 0; m<dim; m++) {
            savedPos[j][k][m] = savedPos[j+1][k][m];
          }
        }
      }
      for (j = 0; j<nParticles; j++) {
        for (k=0; k<dim; k++) {
          savedPos[maxCorrelationSteps-1][j][k] = pos[j][k];
        }
      }
      for(j = 0; j<maxCorrelationSteps; j++){
        for(k = 0; k<nParticles; k++){
          meanDistanceScalar[j][k] += getDistanceSquared( \
          savedPos[0][k], savedPos[j][k]); 
        }
      }
    } 


    if(i%100 == 0) {  // Progress indicator
      printf("I");
      fflush(stdout);
    }

  } // End of production loop

  printf("\tDone!\nCleaning up\n");

// The production part of the simulation is now done.
// The code below forms some averages and other stuff that
// is best done with all data on hand

  FILE *valuesFile;   // For saving various values (non-array data)
  valuesFile = fopen("values.data","w");
  
  nAverageSteps = nProductionSteps-maxCorrelationSteps;
      // take averages over this many steps
  
  FILE *msdFile;
  msdFile = fopen("msd.data","w");

  for( i = 0; i < maxCorrelationSteps; i++) {
    meanDistanceAverage[i] = 0;
    for (j = 0; j<nParticles; j++) {
      meanDistanceScalar[i][j] = meanDistanceScalar[i][j]/nAverageSteps;
      meanDistanceAverage[i] += meanDistanceScalar[i][j];
    }
    meanDistanceAverage[i] = meanDistanceAverage[i]/nParticles;
    fprintf(msdFile, "%e\t%e\n", i*timestep, meanDistanceAverage[i]);
  }

  // Diffusion coeff from MSD
  for (i = maxCorrelationSteps - msdAverageSteps; \
      i < maxCorrelationSteps; i++) {
    msdDiffCoeff += 1.0/(6*i*timestep)*meanDistanceAverage[i];
  }
  msdDiffCoeff = msdDiffCoeff / msdAverageSteps;
  fprintf(valuesFile,"%e\tDiffusion coeff (MSD)\n", msdDiffCoeff);

  // Final processing and saving of velocity correlation
  // function and diffusion coeff from time integral
  FILE *velcorFile;
  velcorFile = fopen("velcor.data","w");

  diffusionCoefficient = 0;

  for( i = 0; i < maxCorrelationSteps; i++) {
    meanVelocityAverage[i] = 0;
    for (j = 0; j<nParticles; j++) {
      meanVelocityScalar[i][j] = meanVelocityScalar[i][j]/nAverageSteps;
      meanVelocityAverage[i] += meanVelocityScalar[i][j];
    }
    meanVelocityAverage[i] = meanVelocityAverage[i]/nParticles;
    fprintf(velcorFile, "%e\t%e\n", i*timestep, meanVelocityAverage[i]);
    
    diffusionCoefficient += meanVelocityAverage[i];
  }

  diffusionCoefficient = (1.0/3.0) * \
  diffusionCoefficient / maxCorrelationSteps;
  fprintf(valuesFile,"%e\tDiffusion coeff (time integral)\n", \
  diffusionCoefficient);

  // Calculating and saving spectrum integral thingy
  FILE *spectrumFile;
  spectrumFile = fopen("spectrum.data","w");

  for( i = 0; i < nSpectrumPoints; i++) {
    spectrum[i] = 0;
    for( j = 0; j < maxCorrelationSteps; j++) {
      spectrum[i] += meanVelocityAverage[j] * cos(spectrumInterval * \
      i * j / nSpectrumPoints);
    }
    spectrum[i] = 2 * spectrum[i]/maxCorrelationSteps;
    fprintf(spectrumFile,"%e \t %e \n",i*spectrumInterval/nSpectrumPoints, \
    spectrum[i]);
  }

  // Setting up for calculating correlation data for temperature and pressure
  meanTemperature = meanTemperature/nAverageSteps;
  meanSquareTemperature = meanSquareTemperature/nAverageSteps;
  meanPressure = meanPressure/nAverageSteps;
  meanSquarePressure = meanSquarePressure/nAverageSteps;

  meanTemperatureSquare = meanTemperature*meanTemperature;
  meanPressureSquare = meanPressure*meanPressure;

  for(i = 0; i < maxCorrelationSteps+1; i++){
    meanTemperature_ik[i] = meanTemperature_ik[i]/nAverageSteps;
    meanPressure_ik[i] = meanPressure_ik[i]/nAverageSteps;
  }

  fprintf(valuesFile, "%e\t Mean Temperature\n", meanTemperature);
  fprintf(valuesFile,"%e\t Mean Pressure\n", meanPressure);

  // Correlation data for temperature and pressure
  FILE *phiTFile;
  phiTFile = fopen("phiTemperature.data","w");

  FILE *phiPFile;
  phiPFile = fopen("phiPressure.data","w");
  
  for(i = 0; i<maxCorrelationSteps; i++){
    phiTemperature[i] = (meanTemperature_ik[i] - meanTemperatureSquare)/ \
    (meanSquareTemperature - meanTemperatureSquare);
    phiPressure[i] = (meanPressure_ik[i] - meanPressureSquare)/ \
    (meanSquarePressure - meanPressureSquare);
    fprintf(phiTFile, "%e \t %e \n", (i+1)*timestep, phiTemperature[i]);
    fprintf(phiPFile, "%e \t %e \n", (i+1)*timestep, phiPressure[i]);
  }

  // Finding M (the longest time for which correlation is above some limit)
  double limit = exp(-2.0);
  i = 0;
  while( (phiTemperature[i] > limit) && (i < maxCorrelationSteps) ) {
    i++;
  }
  MTemperature = i;

  i = 0;
  while( (phiPressure[i] > limit) && (i < maxCorrelationSteps) ) {
    i++;
  }
  MPressure = i;

  //Sum over phi to get statistical inefficiency for temperature
  for(i = 0; i<MTemperature; i++){
    sTemperature += phiTemperature[i];
  }
  sTemperature = sTemperature*2;
  // phi is symmetric and we sum from 0 to M instead of from -M to M

  //Sum over phi to get statistical inefficiency for pressure
  for(i = 0; i<MPressure; i++){
    sPressure += phiPressure[i];
  }
  sPressure = sPressure*2;
  // phi is symmetric and we sum from 0 to M instead of from -M to M

  // Calculating variances
  varTemperature = meanSquareTemperature - meanTemperatureSquare;
  varPressure = meanSquarePressure - meanPressureSquare;

  varTemperature = varTemperature*sTemperature/nProductionSteps;
  varPressure = varPressure*sPressure/nProductionSteps;
  fprintf(valuesFile,"%e\tsTemperature\n", sTemperature);
  fprintf(valuesFile,"%e\tsPressure\n", sPressure);
  fprintf(valuesFile,"%e\tvarTemperature\n", varTemperature);
  fprintf(valuesFile,"%e\tvarPressure\n",varPressure);
  
  FILE *donefile;   // Just an empty file to indicate completion
  donefile = fopen("done.data", "w"); // useful when using ssh+screen
  fprintf(donefile, "done");

  close(energyFile);
  close(ptFile);
  close(positionFile);
  close(valuesFile);
  close(msdFile);
  close(velcorFile);
  close(spectrumFile);
  close(phiTFile);
  close(phiPFile);
  close(donefile);

  printf("Done!\n");
  return 0;
}
