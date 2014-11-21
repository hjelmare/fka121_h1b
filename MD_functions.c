#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define BOLTZMANN 8.6173324*0.00001 // in eV/K


double getDistanceSquared(double vec1[3], double vec2[3])
{
  double sum = 0;
  double temporary;
  int i;
  

  for(i = 0; i<3; i++){
    temporary = vec1[i] - vec2[i];
    sum += pow(temporary,2);
  }
  return sum;
}

double ScalarProduct(double vec1[3], double vec2[3])
{
  double product = 0;
  int i;

  for(i = 0 ; i < 3 ; i++) {
    product += vec1[i]*vec2[i];
  }

  return product;
}


double GetPressure(double temperature, double volume, double virial, \
int nParticles)
{
  double pressure = 0;

  pressure = (nParticles * BOLTZMANN * temperature + virial)/volume;

  return pressure;
}

double GetInstantTemperature(double vel [][3], int nParticles, \
double mass, int dim)


{
  double temperature = 0;
  int i,j;

  double factor = 2/(3*nParticles*BOLTZMANN);

  for (i=0;i<nParticles;i++) {
    for (j=0;j<dim;j++) {
      temperature += factor * pow(vel[i][j],2)*mass/2;
    }
  } 

  return temperature;
}

// Calculates the kinetic energy
double GetKineticEnergy(double vel[][3], double mass, int nParticles)
{
  int i,j;
  double sum;
  double energy = 0;

  for(i = 0; i<nParticles; i++) {
    sum = 0;
    for(j = 0; j<3; j++) {
      sum +=vel[i][j]*vel[i][j];
    }
    energy += sum*mass/2;
  }
  return energy;
}

// Calculate alphaT (the modifier to get the right temperature)
double GetAlphaT(double wantedTemp, double currentTemp, \
                        double timestep, double timeConstant)
{
  double alpha;
  alpha = 1 + timestep/timeConstant * (wantedTemp-currentTemp)/currentTemp;

  return alpha;
}



// Calculate alphaP (the modifier to get the right pressure)
double GetAlphaP(double wantedPressure, double currentPressure, \
                          double timestep, double timeConstant)
{
  double kappa = 2.21901454; // angstrom^3/eV (at 300K)
  double alpha;

  alpha = 1 - kappa*timestep/timeConstant* \
        (wantedPressure - currentPressure);
  
  return alpha;
}
