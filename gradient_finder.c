/* ************************* */
/* 
Routine: gradient_finder
Purpose: find the velcoity gradient via opposite square weight at any given location within the defined sphere
Coordinates: Spherical (radius {0, R}, theta {0, 2*pi}, phi {0, pi})
Unit: SI 
Commenting: Comment above a block indicates the block's purpose; comment next to line provides clarification 
To-do: 
      -Replace arrays structs where appropriate
      -Incorporate Harvey's speed up suggestions
      -Pull redundant lines out into subroutines
      -Implement a more efficient means of comparing particle positions  
      -Set max velocities
*/
/* ************************* */

//include and define statements
#define _USE_MATH_DEFINES //required to use 'M_PI'

#define mass_rand_max 76; //Units: 0.1 mg
#define rad_rand_max 13 //Units: mm
#define theta_rand_max 3 //Units: pi radians
#define phi_rand_max 2 //Units: pi radians
#define vel_r_rand_max 10 //Units: 0.001 mm/s
#define vel_t_rand_max 10 //Units: 0.001 pi radians/s
#define vel_p_rand_max 10 //Units: 0.001 pi radians
#define particle_num 10
#define b_drag 0.005 //Drag coefficient (F_drag = -b*v)

#include <stdio.h>
#include <math.h> //required to use 'M_PI'
#include <stdlib.h> //required to use 'rand()'
#include <time.h> //required to use 'srand(time(NULL))'

//Main routine
int main () {

  //Initializing rand() and setting maximums for rand()
  srand(time(NULL));

  //Set max values
  float mass_max = 7.5; // mg; one-thousandth of average adult human eye mass (7.5 g)
  float radius_max = 12; // average adult human eye diameter 12 mm 
  float theta_max = 2*M_PI; // 2*pi radians
  float phi_max = M_PI; // pi radians
  
  //Declaring mass and velocity arrays
  float masses[particle_num]; 
  
  float pos_radius[particle_num];
  float pos_theta[particle_num];
  float pos_phi[particle_num];
  
  float vel_radius[particle_num];
  float vel_theta[particle_num];
  float vel_phi[particle_num];

  //Setting initial values for all properties
  for(int i=0; i < particle_num; i++)
  {
    //Set masses for all particles
    float mass_temp = rand()%mass_rand_max; //rand() only works with ints
    mass_temp = mass_temp/mass_rand_max; //Gets non-integer values
    masses[i] = 0.0001*mass_temp; //Scales correctly (mg)
    
    //Set initial radial positions for all particles
    pos_radius[i] = rand()%rad_rand_max; //(mm)

    //Set initial angular positions for all particles
    float ang_temp = rand()%theta_rand_max; //rand() only works with ints
    ang_temp = ang_temp/theta_rand_max; //Gets non-integer values
    pos_theta[i] = 2*M_PI*ang_temp; //Scales correctly (radians)
    
    ang_temp = rand()%phi_rand_max; //rand() only works with ints
    ang_temp = ang_temp/phi_rand_max; //Gets non-integer values
    pos_phi[i] = M_PI*ang_temp; //Scales correctly (radians)

    //Set initial radial and angular velocities for all particles; multiplying by 0.001 keeps steps small
    float vel_temp = rand()%vel_r_rand_max; //rand() only works with ints
    vel_temp = vel_temp/vel_r_rand_max; //Gets non-integer values
    vel_radius[i] = 0.001*vel_temp; //Units: mm/s

    vel_temp = rand()%vel_t_rand_max; //rand() only works with ints
    vel_temp = vel_temp/vel_t_rand_max; //Gets non-integer values
    vel_theta[i] = 0.001*M_PI*vel_temp; //Scales correctly (radians/s)

    vel_temp = rand()%vel_p_rand_max; //rand() only works with ints
    vel_temp = vel_temp/vel_p_rand_max; //Gets non-integer values
    vel_phi[i] = 0.001*M_PI*vel_temp; //Scales correctly (radians/s) 
  }

  //Temporal evolution of particles' positions and motion
  float t_max = 10; //Units: seconds
  float t_step;
  
  float temp = particle_num;
  float delta_t_step = 1/temp;

  for(int j=0; j < particle_num*t_max; j++) 
  {
    //Time step position and velocity 
    //Sets index for arrays 
    int i = j/t_max;
    
    //Sets time step
    t_step = delta_t_step*j;
    
    //Evolves particle positions
    float beta = b_drag/masses[i];
    
    pos_radius[i] = pos_radius[i] + vel_radius[i]*t_step - 0.5*beta*vel_radius[i]*t_step*t_step; //Approx. a_drag = -(b/m)*v 
    //Keeps particle within sphere
    if (fabsf(pos_radius[i]) >= radius_max) 
    {
	pos_radius[i] = pos_radius[i] - radius_max;
	vel_radius[i] = -1*vel_radius[i];
    }
      
    pos_theta[i] = pos_theta[i] + vel_theta[i]*t_step - 0.5*beta*vel_theta[i]*t_step*t_step; //Approx. a_drag = -(b/m)*v 
    //Maintain 2*pi periodicity 
    if (pos_theta[i] > theta_max)
    {
	pos_theta[i] = pos_theta[i] - theta_max;
    }
      
    pos_phi[i] = pos_phi[i] + vel_phi[i]*t_step - 0.5*beta*vel_phi[i]*t_step*t_step; //Approx. a_drag = -(b/m)*v 
    //Maintain pi periodicity 
    if (pos_phi[i] > phi_max)
    {
      pos_phi[i] = -1*(2*M_PI - pos_phi[i]);
      pos_theta[i] = pos_theta[i] + M_PI;
    }
    
    //At end of each time step, check if any particles moved to same position
    int a, b;
    for(a=0; a < particle_num; a++)
    {
      for(b=a; b < particle_num; b++)
      {
	if (a != b)
	{
	  if ( (pos_radius[a] == pos_radius[b]) && (pos_theta[a] == pos_theta[b] ) && (pos_phi[a] == pos_phi[b]) )
	  {
	    /*printf("i: %d j: %d delta_t_step: %f t_step: %f\n", i, j, delta_t_step, t_step);
	    printf("a: %d b: %d testing collisions\n", a, b);
	    printf("Before collision:\n");
	    printf("Particle a - v_rad: %f, v_theta: %f, v_phi: %f \n", vel_radius[a], vel_theta[a], vel_phi[a]);
	    printf("Particle b - v_rad: %f, v_theta: %f, v_phi: %f \n", vel_radius[b], vel_theta[b], vel_phi[b]);*/
	    
	    //Transform to center of momentum frame, so p = 0
	    float V_rad_cm = (masses[a]*vel_radius[a] + masses[b]*vel_radius[b])/(masses[a] + masses[b]);
	    float vel_rad_a_cm = vel_radius[a] - V_rad_cm;
	    float vel_rad_b_cm = vel_radius[b] - V_rad_cm;

	    float V_theta_cm = (masses[a]*vel_theta[a] + masses[b]*vel_theta[b])/(masses[a] + masses[b]); //As the R in the ang. mom. is the same for both particles, preemptively diving it out
	    float vel_theta_a_cm = vel_theta[a] - V_theta_cm; 
	    float vel_theta_b_cm = vel_theta[b] - V_theta_cm; 

	    float V_phi_cm = (masses[a]*vel_phi[a] + masses[b]*vel_phi[b])/(masses[a] + masses[b]); //As the R in the ang. mom. is the same for both particles, preemptively diving it out
	    float vel_phi_a_cm = vel_phi[a] - V_phi_cm; 
	    float vel_phi_b_cm = vel_phi[b] - V_phi_cm;

	    //Set initial kinetic energies
	    float KE_radius = 0.5*masses[a]*vel_rad_a_cm*vel_rad_a_cm + 0.5*masses[b]*vel_rad_b_cm*vel_rad_b_cm;
	    float KE_theta = 0.5*masses[a]*vel_theta_a_cm*vel_theta_a_cm + 0.5*masses[b]*vel_theta_b_cm*vel_theta_b_cm; //As the R in the ang. mom. is the same for both particles, and before and after collisions, preemptively diving it out
	    float KE_phi = 0.5*masses[a]*vel_phi_a_cm*vel_phi_a_cm + 0.5*masses[b]*vel_phi_b_cm*vel_phi_b_cm; //As the R in the ang. mom. is the same for both particles, and before and after collisions, preemptively diving it out

	    //Perfectly elastic collision: solve kinetic energy conservation and momentum conservation to get velocities of particles
	    float v_temp = masses[a]*(1 + masses[a]/masses[b]);
	    v_temp = (2*KE_radius)/v_temp;
	    vel_radius[a] = sqrt(v_temp);
	    vel_radius[a] = vel_radius[a] + V_rad_cm; //Transform back to lab frame
	    vel_radius[b] = -1*(masses[a]/masses[b])*vel_radius[a];
	    vel_radius[b] = vel_radius[b] + V_rad_cm; //Transform back to lab frame

	    v_temp = masses[a]*(1 + masses[a]/masses[b]);
	    v_temp = (2*KE_theta)/v_temp;
	    vel_theta[a] = sqrt(v_temp);
	    vel_theta[a] = vel_theta[a] + V_theta_cm; //Transform back to lab frame
	    vel_theta[b] = -1*(masses[a]/masses[b])*vel_theta[a];
	    vel_theta[b] = vel_theta[b] + V_theta_cm; //Transform back to lab frame

	    v_temp = masses[a]*(1 + masses[a]/masses[b]);
	    v_temp = (2*KE_phi)/v_temp;
	    vel_phi[a] = sqrt(v_temp);
	    vel_phi[a] = vel_phi[a] + V_phi_cm; //Transform back to lab frame
	    vel_phi[b] = -1*(masses[a]/masses[b])*vel_phi[a];
	    vel_phi[b] = vel_phi[b] + V_phi_cm; //Transform back to lab frame

	    /*printf("After collision:\n");
	    printf("Particle a - v_rad: %f, v_theta: %f, v_phi: %f \n", vel_radius[a], vel_theta[a], vel_phi[a]);
	    printf("Particle b - v_rad: %f, v_theta: %f, v_phi: %f \n", vel_radius[b], vel_theta[b], vel_phi[b]);*/
	  }
	}
      }
    }
  }
  
  return(0);
}
/* ************************* */
