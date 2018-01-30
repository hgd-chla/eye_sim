#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h> //required to use 'M_PI'
#include <stdlib.h> //required to use 'rand()'
#include <time.h> //required to use 'srand(time(NULL))'

/* ************************* */
/* 
Routine: gradient_finder
Purpose: find the velcoity gradient via opposite square weight at any given location within the defined sphere
Coordinates: Spherical (radius {0, R}, theta {0, 2*pi}, phi {0, pi})
Unit: SI 
To-do: -  
*/
/* ************************* */
int main () {

  //Initializing rand()
  srand(time(NULL));
  int rad_rand_max = 6;
  int theta_rand_max = 3;
  int phi_rand_max = 2;

  //Declaring mass and velocity arrays
  int particle_num = 10;
  float radius_max = 10;
  float theta_max = 2*M_PI;
  float phi_max = M_PI;
  
  float masses[particle_num];
  
  float pos_radius[particle_num];
  float pos_theta[particle_num];
  float pos_phi[particle_num];
  
  float vel_radius[particle_num];
  float vel_theta[particle_num];
  float vel_phi[particle_num];

  for(int i=0; i < particle_num; i++)
  { 
    //Set initial radial positions for all particles
    pos_radius[i] = rand()%rad_rand_max;

    //Set initial angular positions for all particles
    //Normalizing random numbers to theta and phi maxes
    float temp = rand()%theta_rand_max;
    temp = temp/theta_rand_max;
    pos_theta[i] = 2*M_PI*temp; 
    temp = rand()%phi_rand_max;
    temp = temp/phi_rand_max;
    pos_phi[i] = M_PI*temp;

    //Set initial radial and angular velocities for all particles
    //Keeps the velocities small enough that one time step won't take the particle out of the sphere
    temp = rand()%11;
    vel_radius[i] = 0.001*temp;
    temp = rand()%11;
    vel_theta[i] = 0.001*temp;
    temp = rand()%11;
    vel_phi[i] = 0.001*temp;      
  }

  //Temporal evolution
  float t_max = 10;
  float t_step;

  for(int j=0; j < 10*t_max; j++)
  {
    t_step = 0.1*j;
    int i = j/t_max;
  
    pos_radius[i] = pos_radius[i] + vel_radius[i]*t_step;
    if (fabsf(pos_radius[i]) >= radius_max)
    {
	pos_radius[i] = pos_radius[i] - radius_max;
	vel_radius[i] = -1*vel_radius[i];
    }
      
    pos_theta[i] = pos_theta[i] + vel_theta[i]*t_step; 
    //Maintain periodicity 
    if (pos_theta[i] > theta_max)
    {
	pos_theta[i] = pos_theta[i] - theta_max;
    }
      
    pos_phi[i] = pos_phi[i] + vel_phi[i]*t_step;
    //Maintain periodicity 
    if (pos_phi[i] > phi_max)
    {
	pos_phi[i] = pos_phi[i] - phi_max;
    }

    for(int j=0; j < i; j++)
    {
      if((pos_radius[i] == pos_radius[j]) && (pos_theta[i] == pos_theta[j]) && (pos_phi[i] == pos_phi[j]))
      {
	//Case where each particle is moving in a different direction
	if(((vel_radius[i] < 0) && (vel_radius[j] > 0)) || ((vel_radius[i] > 0) && (vel_radius[j] < 0)))
	{
	  vel_radius[i] = -1*vel_radius[i];
	  vel_radius[j] = -1*vel_radius[j];
	  
	  vel_theta[i] = -1*vel_theta[i];
	  vel_theta[j] = -1*vel_theta[j];
	  
	  vel_phi[i] = -1*vel_phi[i];
	  vel_phi[j] = -1*vel_phi[j];
	}
	//Case where both particles are moving in the same direction, but one has a higher velocity
      }
    }
  }
  
  return(0);
}
