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
#define t_max 10 //Units: seconds
#define delta_t 0.1 //Units: seconds
#define n_t 10 //Units: inverse seconds (granularity of temporal evolution)

#include <stdio.h>
#include <math.h> //required to use 'M_PI'
#include <stdlib.h> //required to use 'rand()'
#include <time.h> //required to use 'srand(time(NULL))'

//Subroutines

//Center of Mass frame
float vel_CM(float mass_a, float mass_b, float vel_a, float vel_b)
{
  float vel_CM = (mass_a*vel_a + mass_b*vel_b)/(mass_a + mass_b);
  return vel_CM;
}

//Kinetic Energy
float kinetic_energy(float mass_1, float mass_2, float vel_a_cm, float vel_b_cm)
{
  float KE = 0.5*mass_1*vel_a_cm*vel_a_cm + 0.5*mass_2*vel_b_cm*vel_b_cm;
  return KE;
}
/* ************************* */
//Elastic Collsion Tools
//Final Velocity of Particle A
float elastic_collsion_vel_1(float mass_1, float mass_2, float vel_a, float vel_cm, float KE)
{
  float v_temp = mass_1*(1 + mass_1/mass_2);
  v_temp = (2*KE)/v_temp;
  float vel_1 = sqrt(v_temp);
  vel_1 = vel_a + vel_cm; //Transform back to lab frame
  return vel_1;
}

//Final Velocity of Particle B
float elastic_collsion_vel_2(float mass_1, float mass_2, float vel_a, float vel_b, float vel_cm)
{
  float vel_2 = -1*(mass_1/mass_2)*vel_a;
  vel_2 = vel_b + vel_cm; //Transform back to lab frame
  return vel_2;
}
/* ************************* */
//Gradient Finder Tools
//Derivative
float derivative(float pos_0, float pos_1, float pos_2, float vel_0, float vel_1, float vel_2)
{
  float temp_1 = (vel_1 - vel_0)/(pos_1 - pos_0) ;
  float temp_2 = (vel_2 - vel_1)/(pos_2 - pos_1);
  float derivative = 0.5*(temp_1 + temp_2);
  return derivative; 
}

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
  
  float pos_radius_t[particle_num][t_max*n_t];
  float pos_theta_t[particle_num][t_max*n_t];
  float pos_phi_t[particle_num][t_max*n_t];
  

  float vel_radius_t[particle_num][t_max*n_t];
  float vel_theta_t[particle_num][t_max*n_t];
  float vel_phi_t[particle_num][t_max*n_t];

  //Setting initial values for all properties
  for(int i=0; i < particle_num; i++)
  {
    //Set masses for all particles
    float mass_temp = rand()%mass_rand_max; //rand() only works with ints
    mass_temp = mass_temp/mass_rand_max; //Gets non-integer values
    masses[i] = 0.0001*mass_temp; //Scales correctly (mg)
    
    //Set initial radial positions for all particles
    pos_radius_t[i][0] = rand()%rad_rand_max; //(mm)

    //Set initial angular positions for all particles
    float ang_temp = rand()%theta_rand_max; //rand() only works with ints
    ang_temp = ang_temp/theta_rand_max; //Gets non-integer values
    pos_theta_t[i][0] = 2*M_PI*ang_temp; //Scales correctly (radians)
    
    ang_temp = rand()%phi_rand_max; //rand() only works with ints
    ang_temp = ang_temp/phi_rand_max; //Gets non-integer values
    pos_phi_t[i][0] = M_PI*ang_temp; //Scales correctly (radians)

    //Set initial radial and angular velocities for all particles; multiplying by 0.001 keeps steps small
    float vel_temp = rand()%vel_r_rand_max; //rand() only works with ints
    vel_temp = vel_temp/vel_r_rand_max; //Gets non-integer values
    vel_radius_t[i][0] = 0.001*vel_temp; //Units: mm/s

    vel_temp = rand()%vel_t_rand_max; //rand() only works with ints
    vel_temp = vel_temp/vel_t_rand_max; //Gets non-integer values
    vel_theta_t[i][0] = 0.001*M_PI*vel_temp; //Scales correctly (radians/s)

    vel_temp = rand()%vel_p_rand_max; //rand() only works with ints
    vel_temp = vel_temp/vel_p_rand_max; //Gets non-integer values
    vel_phi_t[i][0] = 0.001*M_PI*vel_temp; //Scales correctly (radians/s) 
  }

  //Let's find the velocity gradient as a function of time
  for(int t=1; t < t_max*n_t; t++)
  {
    for(int i=0; i < particle_num; i++)
    {
      //*************************
      //Evolves particle positions
      //and velocities
      //************************* 
      float beta = b_drag/masses[i];
      float t_step = delta_t*t;
    
      pos_radius_t[i][t] = pos_radius_t[i][t-1] + vel_radius_t[i][t-1]*t_step - 0.5*beta*vel_radius_t[i][t-1]*t_step*t_step; //Approx. a_drag = -(b/m)*v
      vel_radius_t[i][t] = vel_radius_t[i][t-1] - beta*vel_radius_t[i][t-1]*t_step;
      //Keeps particle within sphere
      if (fabsf(pos_radius_t[i][t]) >= radius_max) { pos_radius_t[i][t] = pos_radius_t[i][t] - radius_max; vel_radius_t[i][t] = -1*vel_radius_t[i][t]; }
      
      pos_theta_t[i][t] = pos_theta_t[i][t-1] + vel_theta_t[i][t-1]*t_step - 0.5*beta*vel_theta_t[i][t-1]*t_step*t_step; //Approx. a_drag = -(b/m)*v
      vel_theta_t[i][t] = vel_theta_t[i][t-1] - beta*vel_theta_t[i][t-1]*t_step;
      //Maintain 2*pi periodicity 
      if (pos_theta_t[i][t] > theta_max) { pos_theta_t[i][t] = pos_theta_t[i][t] - theta_max; }
      
      pos_phi_t[i][t] = pos_phi_t[i][t-1] + vel_phi_t[i][t-1]*t_step - 0.5*beta*vel_phi_t[i][t-1]*t_step*t_step; //Approx. a_drag = -(b/m)*v
      vel_phi_t[i][t] = vel_phi_t[i][t-1] - beta*vel_phi_t[i][t-1]*t_step;
      //Maintain pi periodicity 
      if (pos_phi_t[i][t] > phi_max) { pos_phi_t[i][t] = -1*(2*M_PI - pos_phi_t[i][t]); pos_theta_t[i][t] = pos_theta_t[i][t] + M_PI; }
    }
    
    //*************************
    //Check for any collisions
    //*************************
    int a, b;
    for(a=0; a < particle_num; a++)
    {
      for(b=a; b < particle_num; b++)
      {
	if (a != b)
	{
	  if ( (pos_radius_t[a][t] == pos_radius_t[b][t]) && (pos_theta_t[a][t] == pos_theta_t[b][t]) && (pos_phi_t[a][t] == pos_phi_t[b][t]) )
	  {
	    //Transform to center of momentum frame, so p = 0
	    float V_rad_cm = vel_CM(masses[a], masses[b], vel_radius_t[a][t], vel_radius_t[b][t]);
	    float vel_rad_a_cm = vel_radius_t[a][t] - V_rad_cm; float vel_rad_b_cm = vel_radius_t[b][t] - V_rad_cm;
	    
	    float V_theta_cm = vel_CM(masses[a], masses[b], vel_theta_t[a][t], vel_theta_t[b][t]);
	    float vel_theta_a_cm = vel_theta_t[a][t] - V_theta_cm; float vel_theta_b_cm = vel_theta_t[b][t] - V_theta_cm; 
	    
	    float V_phi_cm = vel_CM(masses[a], masses[b], vel_phi_t[a][t], vel_phi_t[b][t]);
	    float vel_phi_a_cm = vel_phi_t[a][t] - V_phi_cm; float vel_phi_b_cm = vel_phi_t[b][t] - V_phi_cm;
	    
	    //Set initial kinetic energies	   
	    float KE_radius = kinetic_energy(masses[a], masses[b], vel_rad_a_cm, vel_rad_b_cm);
	    float KE_theta = kinetic_energy(masses[a], masses[b], vel_theta_a_cm, vel_theta_b_cm);
	    float KE_phi = kinetic_energy(masses[a], masses[b], vel_phi_a_cm, vel_phi_b_cm);
	    
	    //Perfectly elastic collision: solve kinetic energy conservation and momentum conservation to get velocities of particles
	    vel_radius_t[a][t] = elastic_collsion_vel_1(masses[a], masses[b], vel_radius_t[a][t], V_rad_cm, KE_radius);
	    vel_radius_t[b][t] = elastic_collsion_vel_2(masses[a], masses[b], vel_radius_t[a][t], vel_radius_t[b][t], V_rad_cm);
	    
	    vel_theta_t[a][t] = elastic_collsion_vel_1(masses[a], masses[b], vel_theta_t[a][t], V_theta_cm, KE_theta);
	    vel_theta_t[b][t] = elastic_collsion_vel_2(masses[a], masses[b], vel_theta_t[a][t], vel_theta_t[b][t], V_theta_cm);
	    
	    vel_phi_t[a][t] = elastic_collsion_vel_1(masses[a], masses[b], vel_phi_t[a][t], V_phi_cm, KE_phi);
	    vel_phi_t[b][t] = elastic_collsion_vel_2(masses[a], masses[b], vel_phi_t[a][t], vel_phi_t[b][t], V_phi_cm);
	  }
	}
      }
    }
  }

  //*************************
  //Set velocity gradient at
  //a given particle's position
  //*************************

  float grad_v_radius_t[particle_num][t_max*n_t]; //The [i][t] element gives the gradient of the velocity 
  float grad_v_theta_t[particle_num][t_max*n_t];  //at the position of particle i at time t*delta_t seconds
  float grad_v_phi_t[particle_num][t_max*n_t];
  
  for(int t=1; t < (t_max*n_t)-1; t++)
  {
    for(int i=0; i < particle_num; i++)
    {
      grad_v_radius_t[i][t] = derivative(pos_radius_t[i][t-1], pos_radius_t[i][t], pos_radius_t[i][t+1], vel_radius_t[i][t-1], vel_radius_t[i][t], vel_radius_t[i][t+1]);
      grad_v_theta_t[i][t] = derivative(pos_theta_t[i][t-1], pos_theta_t[i][t], pos_theta_t[i][t+1], vel_theta_t[i][t-1], vel_theta_t[i][t], vel_theta_t[i][t+1])/(pos_radius_t[i][t]*sin(pos_phi_t[i][t]));
      grad_v_phi_t[i][t] = derivative(pos_phi_t[i][t-1], pos_phi_t[i][t], pos_phi_t[i][t+1], vel_phi_t[i][t-1], vel_phi_t[i][t], vel_phi_t[i][t+1])/pos_radius_t[i][t];
    }
  }

  for(int i=0; i < particle_num; i++)
  {
    float temp = (vel_radius_t[i][1] - vel_radius_t[i][0])/(pos_radius_t[i][1] - pos_radius_t[i][0]); 
    grad_v_radius_t[i][0] = temp;
    
    temp = (vel_theta_t[i][1] - vel_theta_t[i][0])/(pos_theta_t[i][1] - pos_theta_t[i][0]);
    grad_v_theta_t[i][0] = temp/(pos_radius_t[i][0]*sin(pos_phi_t[i][0]));
    
    temp = (vel_phi_t[i][1] - vel_phi_t[i][0])/(pos_phi_t[i][1] - pos_phi_t[i][0]); 
    grad_v_phi_t[i][0] = temp/pos_radius_t[i][0];

    temp = (vel_radius_t[i][(t_max*n_t)-1] - vel_radius_t[i][(t_max*n_t)-2])/(pos_radius_t[i][(t_max*n_t)-1] - pos_radius_t[i][(t_max*n_t)-2]); 
    grad_v_radius_t[i][(t_max*n_t)-1] = temp;
    
    temp = (vel_theta_t[i][(t_max*n_t)-1] - vel_theta_t[i][(t_max*n_t)-2])/(pos_theta_t[i][(t_max*n_t)-1] - pos_theta_t[i][(t_max*n_t)-2]); 
    grad_v_theta_t[i][(t_max*n_t)-1] = temp/(pos_radius_t[i][(t_max*n_t)-1]*sin(pos_phi_t[i][(t_max*n_t)-1]));
    
    temp = (vel_phi_t[i][(t_max*n_t)-1] - vel_phi_t[i][(t_max*n_t)-2])/(pos_phi_t[i][(t_max*n_t)-1] - pos_phi_t[i][(t_max*n_t)-2]);
    grad_v_phi_t[i][(t_max*n_t)-1] = temp/pos_radius_t[i][(t_max*n_t)-1];    
  }

  /*
  //print testing
  for(int t=1; t < (t_max*n_t)-1; t++)
  {
    printf("t = %d\n", t);
    for(int i=0; i < particle_num; i++)
    {
      printf("i = %d\n", i);
      printf("%f %f %f \n", vel_theta_t[i][t], vel_theta_t[i][t], vel_phi_t[i][t]);
    }
  }
  */
  
  
  return(0);
}
/* ************************* */
