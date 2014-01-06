#ifndef SIM_GRAV_H
#define SIM_GRAV_H

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  SIMS  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Globals //
   float glb_alpha_min, glb_alpha_max, c_par_min, c_par_max;
   float pdet_min, pdet_max;

void sim_inits();

#endif;


// Read in sim inits //
void sim_inits()
{
   string tmp;
   ifstream init_file("inits_sims.ini");

   init_file >> tmp;
   init_file >> glb_alpha_min;
   init_file >> glb_alpha_max;
   init_file >> tmp;
   init_file >> c_par_min;
   init_file >> c_par_max;
   init_file >> tmp;
   init_file >> pdet_min;
   init_file >> pdet_max;
   init_file >> tmp;
   init_file >> fixed_d;
   init_file >> tmp;
   init_file >> sim;
   init_file >> tmp;
   init_file >> fit_pdet;
   init_file >> tmp;
   init_file >> to_year;

   init_file.close();
   if(fixed_d != 0)
      d_par = fixed_d;

}


// Sim lakes and sources //
// <<< Instead just use 2EB

// Initiallize (seed) the largest lake //
// <<< Instead just use actual seed lakes
// <<< But we will need to clear the current records of
//     last_abs, last_obs, last_obs_uninv, discovered


// Simulate the invasion //


// Simulate the detection process //


// Fit the presence only model //


