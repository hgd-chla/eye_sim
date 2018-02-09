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
Bugs:
      -theta and phi runaway 
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
#include <string.h>

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
//Derivative: dv/vx = (1/v)*dv/dt
float derivative(float t_0, float t_1, float t_2, float vel_0, float vel_1, float vel_2)
{
  float del_v_1 = vel_1 - vel_0;
  float del_t_1 = t_1 - t_0;
  float temp_1 = del_v_1/del_t_1;
  
  float del_v_2 = vel_2 - vel_1;
  float del_t_2 = t_2 - t_1;
  float temp_2 = del_v_2/del_t_2;
  
  float derivative = 0.5*(temp_1 + temp_2);
  if( vel_1 == 0.0000 ){ derivative = derivative/0.0000001; } else { derivative = derivative/vel_1; }
  return derivative; 
}

//Main routine
int main () {
  
  //*************************
  //Initializing rand()
  //*************************
  srand(time(NULL));

  //*************************
  //Set max values for rand()
  //*************************
  float mass_max = 7.5; // mg; one-thousandth of average adult human eye mass (7.5 g)
  float radius_max = 12; // average adult human eye diameter 12 mm 
  float theta_max = 2*M_PI; // 2*pi radians
  float phi_max = M_PI; // pi radians
  
  //*************************
  //Declare arrays 
  //*************************
  float masses[particle_num]; 
  
  float pos_radius_t[particle_num][t_max*n_t];
  float pos_theta_t[particle_num][t_max*n_t];
  float pos_phi_t[particle_num][t_max*n_t];
  
  float vel_radius_t[particle_num][t_max*n_t];
  float vel_theta_t[particle_num][t_max*n_t];
  float vel_phi_t[particle_num][t_max*n_t];

  //*************************
  //Randomly fill velocities and positions
  //*************************
  for(int t=0; t < t_max*n_t; t++)
  {
    for(int i=0; i < particle_num; i++)
    {
      //Set initial radial positions for all particles
      pos_radius_t[i][t] = rand()%rad_rand_max; //(mm)

      //Set initial angular positions for all particles
      float ang_temp = rand()%theta_rand_max; //rand() only works with ints
      ang_temp = ang_temp/theta_rand_max; //Gets non-integer values
      pos_theta_t[i][t] = 2*M_PI*ang_temp; //Scales correctly (radians)
    
      ang_temp = rand()%phi_rand_max; //rand() only works with ints
      ang_temp = ang_temp/phi_rand_max; //Gets non-integer values
      pos_phi_t[i][t] = M_PI*ang_temp; //Scales correctly (radians)

      //Set initial radial and angular velocities for all particles; multiplying by 0.001 keeps steps small
      float vel_temp = rand()%vel_r_rand_max; //rand() only works with ints
      vel_temp = vel_temp/vel_r_rand_max; //Gets non-integer values
      vel_radius_t[i][t] = 0.001*vel_temp; //Units: mm/s

      vel_temp = rand()%vel_t_rand_max; //rand() only works with ints
      vel_temp = vel_temp/vel_t_rand_max; //Gets non-integer values
      vel_theta_t[i][t] = 0.001*M_PI*vel_temp; //Scales correctly (radians/s)

      vel_temp = rand()%vel_p_rand_max; //rand() only works with ints
      vel_temp = vel_temp/vel_p_rand_max; //Gets non-integer values
      vel_phi_t[i][t] = 0.001*M_PI*vel_temp; //Scales correctly (radians/s) 
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
      float temp = derivative(delta_t*(t-1), delta_t*t, delta_t*(t+1), vel_radius_t[i][t-1], vel_radius_t[i][t], vel_radius_t[i][t+1]);
      grad_v_radius_t[i][t] = temp;
      
      temp = derivative(delta_t*(t-1), delta_t*t, delta_t*(t+1), vel_theta_t[i][t-1], vel_theta_t[i][t], vel_theta_t[i][t+1]);
      if ( (pos_radius_t[i][t]*sin(pos_phi_t[i][t])) == 0.0000 ) { grad_v_theta_t[i][t] = temp/0.0000001;} else { grad_v_theta_t[i][t] = temp/(pos_radius_t[i][t]*sin(pos_phi_t[i][t])) ;}

      temp = derivative(delta_t*(t-1), delta_t*t, delta_t*(t+1), vel_phi_t[i][t-1], vel_phi_t[i][t], vel_phi_t[i][t+1]);
      if (pos_radius_t[i][t] == 0.0000) { grad_v_phi_t[i][t] = temp/0.0000001; } else { grad_v_phi_t[i][t] = temp/pos_radius_t[i][t]; }
    }
  }

  //Set velocity gradient at initial and final time steps
  for(int i=0; i < particle_num; i++)
  {
    float temp = (vel_radius_t[i][1] - vel_radius_t[i][0])/delta_t;
    if(vel_radius_t[i][0] == 0.0000){ temp = temp/0.0000001; } else { temp = temp/vel_radius_t[i][0]; }
    grad_v_radius_t[i][0] = temp;
    
    temp = (vel_theta_t[i][1] - vel_theta_t[i][0])/delta_t;
    if(vel_theta_t[i][0] == 0.0000){ temp = temp/0.0000001; } else { temp = temp/vel_theta_t[i][0]; }
    if( (pos_radius_t[i][0]*sin(pos_phi_t[i][0])) == 0.0000) { grad_v_theta_t[i][0] = temp/0.0000001; } else { grad_v_theta_t[i][0] = temp/(pos_radius_t[i][0]*sin(pos_phi_t[i][0])); }
    
    temp = (vel_phi_t[i][1] - vel_phi_t[i][0])/delta_t;
    if(vel_phi_t[i][0] == 0.0000){ temp = temp/0.0000001; } else { temp = temp/vel_phi_t[i][0]; } 
    if(pos_radius_t[i][0] == 0.0000){ grad_v_phi_t[i][0] = temp/0.0000001; } else { grad_v_phi_t[i][0] = temp/pos_radius_t[i][0]; }

    temp = (vel_radius_t[i][(t_max*n_t)-1] - vel_radius_t[i][(t_max*n_t)-2])/( delta_t*( ((t_max*n_t)-2) - ((t_max*n_t)-1) ) );
    if(vel_radius_t[i][(t_max*n_t)-1] == 0.0000){ temp = temp/0.0000001; } else { temp = temp/vel_radius_t[i][(t_max*n_t)-1]; }
    grad_v_radius_t[i][(t_max*n_t)-1] = temp;
    
    temp = (vel_theta_t[i][(t_max*n_t)-1] - vel_theta_t[i][(t_max*n_t)-2])/( delta_t*( ((t_max*n_t)-2) - ((t_max*n_t)-1) ) );
    if(vel_theta_t[i][(t_max*n_t)-1] == 0.0000){ temp = temp/0.0000001; } else { temp = temp/vel_theta_t[i][(t_max*n_t)-1]; }
    if( (pos_radius_t[i][(t_max*n_t)-1]*sin(pos_phi_t[i][(t_max*n_t)-1])) == 0.0000 ){ grad_v_theta_t[i][(t_max*n_t)-1] = temp/0.0000001; } else { grad_v_theta_t[i][(t_max*n_t)-1] = temp/(pos_radius_t[i][(t_max*n_t)-1]*sin(pos_phi_t[i][(t_max*n_t)-1])); }
    
    temp = (vel_phi_t[i][(t_max*n_t)-1] - vel_phi_t[i][(t_max*n_t)-2])/( delta_t*( ((t_max*n_t)-2) - ((t_max*n_t)-1) ) );
    if(vel_phi_t[i][(t_max*n_t)-1] == 0.0000){ temp = temp/0.0000001; } else { temp = temp/vel_phi_t[i][(t_max*n_t)-1]; }
    if(pos_radius_t[i][(t_max*n_t)-1] == 0.0000 ){grad_v_phi_t[i][(t_max*n_t)-1] = temp/0.0000001;} else {grad_v_phi_t[i][(t_max*n_t)-1] = temp/pos_radius_t[i][(t_max*n_t)-1];}
       
  }

  //*************************
  //Write the velocities,
  //positions, and velocity gradient
  //for each particle to separate files
  //*************************
  for(int i=0; i < particle_num; i++)
  {
    
    char pos_file_name [17];
    strcpy(pos_file_name,"pos_output_");

    int length = snprintf(NULL, 0, "%d", i);
    char* str = malloc( length + 1 );
    snprintf(str, length + 1, "%d", i);

    strcat(pos_file_name,str);
    strcat(pos_file_name,".txt");
   
    FILE *f_pos = fopen(pos_file_name, "w");
    if (f_pos == NULL)
    {
      printf("Error opening file!\n");
      exit(1);
    }

    char vel_file_name [17];
    strcpy(vel_file_name,"vel_output_");

    int length2 = snprintf(NULL, 0, "%d", i);
    char* str2 = malloc( length2 + 1 );
    snprintf(str2, length2 + 1, "%d", i);

    strcat(vel_file_name,str2);
    strcat(vel_file_name,".txt");
   
    FILE *f_vel = fopen(vel_file_name, "w");
    if (f_vel == NULL)
    {
      printf("Error opening file!\n");
      exit(1);
    }

    char grad_file_name [18];
    strcpy(grad_file_name,"grad_output_");

    int length3 = snprintf(NULL, 0, "%d", i);
    char* str3 = malloc( length3 + 1 );
    snprintf(str3, length3 + 1, "%d", i);

    strcat(grad_file_name,str3);
    strcat(grad_file_name,".txt");
   
    FILE *f_grad = fopen(grad_file_name, "w");
    if (f_grad == NULL)
    {
      printf("Error opening file!\n");
      exit(1);
    }
    
    for(int t=0; t < (t_max*n_t)-1; t++)
    {
      fprintf(f_pos,"%f %.5f %.5f %.5f \n", delta_t*t, pos_radius_t[i][t], pos_theta_t[i][t], pos_phi_t[i][t]);
      fprintf(f_vel,"%f %.5f %.5f %.5f \n", delta_t*t, vel_radius_t[i][t], vel_theta_t[i][t], vel_phi_t[i][t]);
      fprintf(f_grad,"%f %.5f %.5f %.5f \n", delta_t*t, grad_v_radius_t[i][t], grad_v_theta_t[i][t], grad_v_phi_t[i][t]);
    }
    fclose(f_pos);
    fclose(f_vel);
    fclose(f_grad);
  }
  
  return(0);
}
/* ************************* */
