//include and define statements
#define _USE_MATH_DEFINES //required to use 'M_PI'

#define tracer_num 10
#define mass_rand_max 76 //Units: 0.1 mg
#define rad_rand_max 13 //Units: mm
#define theta_rand_max 3 //Units: pi radians
#define phi_rand_max 2 //Units: pi radians
#define vel_r_rand_max 10 //Units: 0.001 mm/s
#define vel_t_rand_max 10 //Units: 0.001 pi radians/s
#define vel_p_rand_max 10 //Units: 0.001 pi radians
#define b_drag 0.005 //Drag coefficient (F_drag = -b*v)
#define t_max 10 //Units: seconds
#define delta_t 0.1 //Units: seconds
#define n_t 10 //Units: inverse seconds (granularity of temporal evolution)

#include <stdio.h>
#include <math.h> //required to use 'M_PI'
#include <stdlib.h> //required to use 'rand()'
#include <time.h> //required to use 'srand(time(NULL))'
#include <string.h>

//Data Structure Definitions
struct tracer
{ 
  float position[t_max*n_t][3]; //[time step][0: rho (mm), 1: theta (radians), 2: phi (radians)]
  float velocity[t_max*n_t][3]; //[time step][0: rho (mm), 1: theta (radians), 2: phi (radians)]
  float gradient[t_max*n_t][3]; //[time step][0: radial gradient, 1: azimuthal gradient, 2: polar gradient]

  //We'll use these later for actual evolution, but we're just testing gradient finder right now
  //float mesh[t_max*n_t][3]; //[time step][0: delta rho (mm), 1: delta theta (radians), 2: delta phi (radians)] 
  //float density; //fluid density within mesh (g*mm^-3)
};

//Routine Defintions

//An easier to read (pseudo)random number generator
float random_gen(int lower, int upper)
{
  int num = (rand() + lower) % (upper + 1);
  return num;    
}

//Assigns a random value to argument val 
void random_fill(float *val, int max, float scale)
{ 
  float temp = random_gen(0, max); //random only works with ints
  temp = temp/max; //Gets non-integer values
  *val = scale*temp; //Scales correctly
}

//*************************
//TO DO: Pass pointer to
//array, not array itself
//*************************
//Calculates spatial derivative using chain rule: dv/vx = (1/v)*dv/dt
float derivative(float velocity[t_max*n_t][3], int t, int i)
{
  float der = 0;
  if (t == 0) {
    
    float del_v_1 = velocity[t+1][i] - velocity[t][i];
    der = del_v_1/delta_t;
    
  } else if (t == (t_max*n_t)-1) {

    float del_v_1 = velocity[t][i] - velocity[t-1][i];
    der = del_v_1/delta_t;
 
  } else {
    
    float del_v_1 = velocity[t][i] - velocity[t-1][i];
    float der_1 = del_v_1/delta_t; 
  
    float del_v_2 = velocity[t+1][i] - velocity[t][i];
    float der_2 = del_v_1/delta_t;

    der = 0.5*(der_1 + der_2);

  }

  if( velocity[t][i] == 0.0000 ){ der = der/0.0000001; } else { der = der/velocity[t][i]; }
    
  return der; 
}

//The titular routine
void gradient(struct tracer *tr, int t)
{
  float der[3];
  for(int i=0; i<3; i++)
  {
    der[i] = derivative(tr->velocity, t, i);
  }

  float scale = 1;
  tr->gradient[t][0] = scale*der[0]; //radial gradient at the position of the tracer at time t*delta_t

  if ( ( (tr->position[t][0])*(sin(tr->position[t][2])) ) == 0.0000 ) { scale = 1/0.0000001;} else { scale = 1/( (tr->position[t][0])*(sin(tr->position[t][2])) ) ;} 
  tr->gradient[t][1] = scale*der[1]; //azimuthal gradient at the position of the tracer at time t*delta_t

  if ( (tr->position[t][0]) == 0.0000) { scale = 1/0.0000001; } else {scale = 1/(tr->position[t][0]); }
  tr->gradient[t][2] = scale*der[2]; //polar gradient at the position of the tracer at time t*delta_t  
}


//Main
int main()
{
  srand(time(NULL)); //initializing rand()

  //Set max values for rand()
  float mass_max = 7.5; // mg; one-thousandth of average adult human eye mass (7.5 g)
  float radius_max = 12; // average adult human eye diameter 12 mm 
  float theta_max = 2*M_PI; // 2*pi radians
  float phi_max = M_PI; // pi radians
  
  struct tracer compendium[tracer_num];

  //Randomly fill tracer positions and velocities
  for(int i=0; i<10; i++)
  {
    for(int t=0; t<t_max*n_t; t++)
    {
      random_fill(&compendium[i].position[t][0], rad_rand_max, 1); // (mm)
      random_fill(&compendium[i].position[t][1], theta_rand_max, 2*M_PI); // (radians)
      random_fill(&compendium[i].position[t][2], phi_rand_max, M_PI); // (radians)

      random_fill(&compendium[i].velocity[t][0], vel_r_rand_max, 0.001); // (mm/s)
      random_fill(&compendium[i].velocity[t][1], vel_t_rand_max, 0.001*2*M_PI); // (radians/s)
      random_fill(&compendium[i].velocity[t][2], vel_p_rand_max, 0.001*M_PI); // (radians/s)
    }
  }

  //Find velocity gradient at each tracer's position
  for(int t=0; t < t_max*n_t; t++)
  {
    for(int i=0; i < tracer_num; i++)
    {
      gradient(&compendium[i], t);
    }
  }


  //*************************
  //TO DO: Rewrite as a routine
  //*************************
  //Write the velocities, positions, and velocity gradients for each tracer to separate files
  for(int i=0; i < tracer_num; i++)
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
      fprintf(f_pos,"%f %.5f %.5f %.5f \n", delta_t*t, compendium[i].position[t][0], compendium[i].position[t][1],compendium[i].position[t][2]);
      fprintf(f_vel,"%f %.5f %.5f %.5f \n", delta_t*t, compendium[i].velocity[t][0], compendium[i].velocity[t][1],compendium[i].velocity[t][2]);
      fprintf(f_grad,"%f %.5f %.5f %.5f \n", delta_t*t, compendium[i].gradient[t][0], compendium[i].gradient[t][1],compendium[i].gradient[t][2]);
    }
    fclose(f_pos);
    fclose(f_vel);
    fclose(f_grad);
}

  
  return 0;
}
