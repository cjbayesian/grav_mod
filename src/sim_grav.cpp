//////////////////////////////////////////////////////////
// Simulates and Tests a Bayesian implimentation of 
// gravity mediated spread of invasive species
// using a presence only observation model
// 
// 	Corey Chivers, 2014
// 	McGill University
//	Department of Biology
//	corey.chivers@mail.mcgill.ca
//
// See Leung et al. 2006 for underlying 
// gravity and likelihood equation defs
// 
//////////////////////////////////////////////////////////

#include"grav.h"
#include"sim_grav.h"

// Read in sim inits //
int main(int argc,char *argv[])
{
   // decode arguments
	args(argc,argv);

   // Sim lakes and sources //
   // <<< Instead just use 2EB
   inits();
   sim_inits();

   // Read in Data // 
   cout<< "Reading Spatial Data\n";
   read_data();

   // Set initial Params //
   cout<< "Initiallizing...\n";
   // Initiallize (seed) the largest lake //
   // <<< Instead just use actual seed lakes
   // <<< But we will need to clear the current records of
   //     last_abs, last_obs, last_obs_uninv, discovered

   //for faster sims:
   n_lakes = 500;
   n_sources = 100;
   clear_traf_mat();


   n_sampled = 0;
   // Drop situations where very few lakes are observed (< %10)
   //while(n_sampled < 30)
   //{
      for(int i=1;i<=n_lakes;i++)
      {
         //clear all but the seed lake(s)
         if(lakes(i).discovered != from_year)
         {
            lakes(i).discovered = 0;
            lakes(i).invaded = 0;
            lakes(i).last_abs = 0;
            lakes(i).status_2010 = 0;
         }
      }

      init_state();
      init_t(); 
      calc_state();
      generate_params();
      gamma_par = 0;
      // Simulate the invasion //
      calc_traf();
      calc_Uj();
      sim_spread(); cout << "\n";   
      // Simulate the detection process //
      detect();
      which_sampled_or_valid();
   //}

   // Fit the presence only model //
   pdet = 0;
   //glb_alpha = glb_alpha*10;
   for(int i=1;i<=39;i++)
   {
      pdet = pdet + 0.025; 
      for(int ss=1;ss<=5;ss++)
      {
         sim_spread();
         cout << pdet << "\t" << l_hood_detp() << "\n";
      }
   }

   return 0;
}


