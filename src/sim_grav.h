#ifndef SIM_GRAV_H
#define SIM_GRAV_H

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  SIMS  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Globals //
   float glb_alpha_min, glb_alpha_max, c_par_min, c_par_max;
   float d_par_min, d_par_max, e_par_min, e_par_max, pdet_min, pdet_max;

void sim_inits();
void generate_params();
#endif;


// Read in sim inits //
void sim_inits()
{
   string tmp;
   ifstream init_file("inits_sim.ini");

   init_file >> tmp;
   init_file >> glb_alpha_min;
   init_file >> glb_alpha_max;
   init_file >> tmp;
   init_file >> c_par_min;
   init_file >> c_par_max;
   init_file >> tmp;
   init_file >> d_par_min;
   init_file >> d_par_max;
   init_file >> tmp;
   init_file >> e_par_min;
   init_file >> e_par_max;
   init_file >> tmp;
   init_file >> pdet_min;
   init_file >> pdet_max;

   init_file.close();
}

void generate_params()
{
   glb_alpha = runif(glb_alpha_min,glb_alpha_max);
   c_par = runif(c_par_min,c_par_max);
   d_par = runif(d_par_min,d_par_max);
   e_par = runif(e_par_min,e_par_max);
   pdet = runif(pdet_min,pdet_max);

   cout << glb_alpha << "\t" << c_par << "\t" << d_par << "\t" << e_par << "\t" << pdet << "\n";
}

void detect()
{
   // Which sites to monitor
   int n_monitored = n_lakes;
   _vbc_vec<int> lake_inds(1,n_lakes);
   for(int i=1;i<=n_lakes;i++) {lake_inds(i)=i;}
   lake_inds = sample_wo_replace(lake_inds,n_monitored);

   // Monitor until detected or to_year
   int lake_index;
   for(int i=1;i<=n_monitored;i++)
   {
      lake_index = lake_inds(i);
      for(int t=t_vec(lake_index);t<=to_year;t++)
      {
         if(runif(0,1) < pdet && lakes(lake_index).discovered != from_year)
         {
            //cout << t_vec(lake_index) << "\t" << t << "\n"; 
            lakes(lake_index).discovered = t;
            break;
         }
      }

   }
   //cout << "from_year: " << from_year << "\n";
}




