// To Compile use:
// g++ -w -I/home/cchivers/c++/library/ -I/usr/share/R/include -o ../gb *.cpp /home/cchivers/c++/library/libCJ -lRmath -lm -fopenmp -O3 -ftree-vectorize -msse2
// g++ -w -I/home/cchivers/SchoolBackUp/c++/library/ -I/usr/share/R/include -o ../gb *.cpp /home/cchivers/SchoolBackUp/c++/library/libCJ -lRmath -lm -fopenmp -O3 -ftree-vectorize -msse2

// Simulates and Tests a Bayesian implimentation of 
// gravity mediated spread of invasive species
// 
// 	Corey Chivers, 2012
// 	McGill University
//	Department of Biology
//	corey.chivers@mail.mcgill.ca
//
// See Leung et al. 2006 for underlying 
// gravity and likelihood equation defs
// 


// * CREATE A MAP OF SUITABILITY * 
// 

#include"grav.h"

int main(int argc,char *argv[])
{
 	// decode arguments
	args(argc,argv);
   inits();
   if(run_type==6)
      to_year=2010;

    // Read in Data //
    cout<< "Reading Data\n";
    read_data();

    // Set initial Params //
    cout<< "Initiallizing...\n";
    init_state();
    init_t(); 
    calc_state();
    calc_traf();
    calc_traf_mat();
    calc_pp(); //!!


   ofstream ll_("output/ll_test.dat");
   if(FALSE)
   {
      int tmp_t;
       for(int lake=1;lake<=n_lakes;lake++)
       {
           if(lakes(lake).invaded == 0 && lakes(lake).last_abs == 0 )
           {
               for(int t=from_year;t<=to_year+2;t++)
               {
                   // sample_t();
                   tmp_t=t_vec(lake);
                   t_vec(lake)=t;
                   calc_state();
                   calc_pp();
                   cout << lake << "\t" << t <<endl;
                   ll_ << lake << "\t" << t <<"\t"<<l_hood()<<endl;
               }
               t_vec(lake)=tmp_t;
           }
       }
   }

   if(FALSE) // likelihood profile of d and e //if init_t is random, MLE d=e=0 (no effect of distance or size: CHECK)
   {
       d_par=-1;
       e_par=0;
       for(int i=0;i<=60;i++)
       {
           d_par=d_par+0.1;
           e_par=0;
           for(int j=0;j<=10;j++)
           {
           e_par=e_par + 0.1;
           calc_traf();
           calc_traf_mat();
           calc_pp();
           ll_ << d_par << "\t" << e_par << "\t" << l_hood()<< endl;
           cout <<  d_par <<"\t"<< e_par << "\t" << l_hood()<< endl;
           }
       }
   }
   if(FALSE) //likelihood profile of alpha
   { 
      chem_pars(1)=0;
       for(int i=0;i<=1000;i++)
       {
           chem_pars(1)+=0.0001;
           ll_ << chem_pars(1) << "\t" << l_hood()<< endl;
       }
   }

   //Just sim under the true alpha to see the t_vec distribution
   // to compare with sim_dat.R
   if(FALSE)
   {
       for(int lake_index=1;lake_index<=n_lakes;lake_index++)
       {
           lakes(lake_index).discovered=0;
           lakes(lake_index).last_abs=0;
       }
       for(int i=1;i<=1000;i++)
       {
           sim_spread();
           write_t();
       }
   }



   /// SEEDS ///
   chem_pars(1)=-8;
   chem_pars(2)=0;
   chem_pars(3)=0; //-2:1 MLE ~0
   chem_pars(4)=0; //-1.5:1.5 MLE ~0
   chem_pars(5)=0; //-0.4:0.2 MLE ~0
   chem_pars(6)=0.05; //0:0.1 MLE 0.05    ***
   chem_pars(7)=1.3; //-0.4:0.2 MLE ~1.3  **** 
   chem_pars(8)=-0.31; //-0.7:0.1 MLE ~-0.31 ****
   chem_pars(9)=-0.01; //-0.4:0.2 MLE ~-0.01 **
   chem_pars(10)=0.06; //-0.1:0.1 MLE ~0.06 *
   chem_pars(11)=0; //-0.16:0.1 MLE ~0
   chem_pars(12)=0; //-0.2:0.2 MLE ~0
   chem_pars(13)=0; //-0.04:0.01 MLE ~0
   chem_pars(14)=0.1; //-0.3:0.3 MLE ~0.1

   int test_ch=14;
   d_par=0.175;
   //float bb = l_hood();

   if(ll)
   {
       float tmplhood;
       for(int i=1;i<=20;i++)
       {
          cerr << i << "\n";
           //chem_pars(test_ch)=chem_pars(test_ch)+0.0013;
           d_par=d_par+0.05;
           calc_traf();
          cerr << "A" << "\n";
           calc_traf_mat();
          cerr << "B" << "\n";

           sim_spread();
          cerr << "C" << "\n";
           //write_t();
           tmplhood = l_hood();
          cerr << "D" << "\n";
           ll_ << chem_pars(test_ch) << "\t" << tmplhood << "\n";
           cout << d_par << "\t" << tmplhood << "\n";

       }
   }

   ll_.close();

cout << "# sampled\t" << n_sampled << "\n";

if(run_type==1)
{
   // FIT ON TRAF_PARS & SPREAD PARS ONLY (NO ENV) //
   // need a likelihood function wrapper to call l_hood() multiple times and average the result
   // to smooth out stochastic surface.
   // BOOTSTRAP RESAMPLING OF DATA (SAMPLED LAKES) TO GENERATE CI //
   float garbage=l_hood();
   int n_reps = 100;
   int n_pars=4; //no e (area exponent)
   _vbc_vec<float>params1(1,n_pars);
   params1(1)=1.27;
   params1(2)=1.48;
   params1(3)=0.0000489;
   params1(4)=0.00105;  
   _vbc_vec<float> dat1(1,n_pars);
   _vbc_vec<float> MLE_params(1,n_pars);

   _vbc_vec<int> tmp_index_sampled;
   tmp_index_sampled = sampled_index;
   ofstream par_file;
   par_file.open("output/pred_pars.tab");

   if(boot)
   {
      ofstream boot_file;
      boot_file.open("output/boot_lakes.tab");
      for(int i=1;i<=n_reps;i++)
      {
         //Bootstrap resample //
         sampled_index = sample_w_replace(tmp_index_sampled);
         for(int j=1;j<=n_sampled;j++)
            boot_file << sampled_index(j) << "\t";
         boot_file << "\n";
         boot_file.flush();
         // --- //

         simplex::clsSimplex<float> gertzen_rep;
         //gertzen_rep.set_param_small(1e-3);
         gertzen_rep.start(&dat1,&params1, &MLE_l_hood,n_pars, 1e-1);
         gertzen_rep.getParams(&MLE_params);

         cout << "\n\nMLE "<< i << " of " << n_reps << "\n\n";
         for(int p=1;p<=n_pars;p++)
            par_file << MLE_params(p) <<"\t";
         par_file << "\n";
         par_file.flush();
      }
      boot_file.close();
   }else{
      simplex::clsSimplex<float> gertzen_rep;
      //gertzen_rep.set_param_small(1e-3);
      gertzen_rep.start(&dat1,&params1, &MLE_l_hood,n_pars, 1e-1);
      gertzen_rep.getParams(&MLE_params);

      cout << "\n\nMLE\n";
      for(int p=1;p<=n_pars;p++)
         par_file << MLE_params(p) <<"\t";
      par_file << "\n";
      par_file.flush();
   }
   par_file.close();

}



if(run_type==2)
{
   //MCMC lib
	string mcmc_file("output/lib.mcmc");
   if(!no_env)
   {
      _vbc_vec<float> params(1,4+n_chem_var);
      _vbc_vec<float> prop_width(1,4+n_chem_var,1,4+n_chem_var);
      prop_width(1)=0.001;
      prop_width(2)=0.05;
      prop_width(3)=0.01;
      prop_width(4)=0.000001;


      params(1)=1;
      params(2)=1;
      params(3)=1;   
      for(int i=1;i<=n_chem_var+1;i++)
      {
         prop_width(i+3)=0.000001;
         params(i+3)=chem_pars(i);
      }
      prop_width(4)=0.1;

      _vbc_vec<float> prop_sigma;
      prop_sigma = diag(prop_width);

      // Print out prop_sigma
      for(int i=1;i<=n_chem_var+4;i++)
      {
         for(int j=1;j<=n_chem_var+4;j++)   
            cout << prop_sigma(i,j) << " | ";
         cout << "\n";
      }
     
      mcmcMD::run_mcmc(params, 
         prop_sigma, 
         &likelihood_wrapperMCMC_MD,
         &prior_MD, 
         &restrict_MCMC_MD, 
         50000, 
         1, 
         1, 
         mcmc_file.c_str(),
         true,
         true,
         true,
         500,
         4);
   }else
   {
   /// No env.
      _vbc_vec<float> params(1,4);
      _vbc_vec<float> prop_width(1,4,1,4);
      prop_width(1)=0.05;
      prop_width(2)=0.05;
      prop_width(3)=0.05;
      prop_width(4)=0.000001;


      params(1)=1;
      params(2)=1;
      params(3)=2;   
      params(4)=0.0001;         

      _vbc_vec<float> prop_sigma;
      prop_sigma = diag(prop_width);

      // Print out prop_sigma
      for(int i=1;i<=4;i++)
      {
         for(int j=1;j<=4;j++)   
            cout << prop_sigma(i,j) << " | ";
         cout << "\n";
      }
     
      mcmcMD::run_mcmc(params, 
         prop_sigma, 
         &likelihood_wrapperMCMC_MD,
         &prior_MD, 
         &restrict_MCMC_MD, 
         50000, 
         1, 
         1, 
         mcmc_file.c_str(),
         true,
         true,
         true,
         500,
         4);
   }


/*	mcmcMD::run_mcmc(pms, 
      props, 
      &like, 
      &prior, 
      &restrictions, 
      500000, 
      1, 
      100, 
      file_name.c_str(),
      TRUE,
      TRUE,
      1000);
*/
}


/// Sim from posterior ///
if(run_type==3)
   sim_spread_posterior();

/// Traf tests ////
if(run_type==4)
{
   ofstream traf_ll_file("output/traf_ll.dat");

   d_par=1;
   e_par=0;

   calc_traf();
   calc_traf_mat();
   calc_pp();
   sim_spread();
   cout<< l_hood() <<"\n";

   for(int i=1;i<=10;i++)
   {
      e_par=e_par+0.2;
      c_par=0;
         calc_traf();
         calc_traf_mat();
         calc_pp();
      for(int j = 1;j<=10;j++)
      {
         c_par=c_par+0.2;
         glb_alpha=0; //from MLE
         for(int k=1;k<=500;k++)
         {
            glb_alpha=glb_alpha+0.00001;
            sim_spread();
            sim_spread();
            sim_spread();
            cout << e_par << "\t" << c_par << "\t" << glb_alpha << "\t" << l_hood()  << "\n"; 
            traf_ll_file << e_par << "\t" << c_par << "\t" << glb_alpha << "\t" << l_hood()  << "\n"; 
         }
      }
   } 

   traf_ll_file.close();
   calc_pp();
   write_pp(); 
   write_inv_stat();
}

/// Just whatever tests ////
if(run_type==5)
{  
   ofstream ll_file("output/ll.dat");
   int n_pars=5;
   _vbc_vec<float>params1(1,n_pars);
   params1(1)=1.27; //according to gertzen 2011
   params1(2)=1;
   params1(3)=2;
   params1(4)=0;  
   params1(5)=0;

   calc_traf();
   calc_traf_mat();
   write_traf_mat();
   /*
   _vbc_vec<float> dat1(1,n_pars);

   for(int i=1;i<=30;i++)
   {
      params1(5) +=0.005;
      for(int p=1;p<=n_pars;p++)
         ll_file << params1(p) << "\t";
      ll_file << MLE_l_hood(&params1, &dat1) << "\n";
   }

   ll_file.close();

   /*
   ofstream n_inv_file("output/n_inv_file.dat");
   e_par=1;
   d_par=1;
   c_par=0;
   glb_alpha=0.005;

   calc_traf();
   calc_traf_mat();
   calc_pp();


   for(int j=1;j<=80;j++)
   {
      e_par=e_par+0.05;

   calc_traf();
   calc_traf_mat();
   calc_pp();

      for(int i=1;i<=30;i++)
      {
         sim_spread();

         //Write cummulative_total invaded for that sim.
         for(int t=from_year+1;t<=to_year;t++)
            n_inv_file << state(t).n_inv << "\t";
         n_inv_file << "\n";

         ll_file<< e_par << "\t" << l_hood() << "\n";
         cout << e_par << "\t"<< l_hood() << "\n";
      }
   }

   n_inv_file.close();
   ll_file.close();

    // MLE at d=e=1,c=1.5 : glb_alpha=0.00081
   
   e_par=0.25;
   d_par=1;
   c_par=1.5;
   glb_alpha=0.005;

   /*
   calc_traf();
   calc_traf_mat();
   calc_pp();
      sim_spread();
      sim_spread();
      sim_spread();
   cout << "Simming...\n";
   for(int i=1;i<=10000;i++)
   {
      sim_spread();
      write_t();
   }
     
   t_file.close();
   */

}
/// Predictions from Bootstapped MLE ////
if(run_type==6)
{  
   //m reps of the bootstrap/posterior
   int m_pars = wc_l("output/pred_pars.tab");
   int n_pars = 4;
   cout << "# Generating a " << n_val_lakes << " by " << m_pars << " prediction matrix\n";

   _vbc_vec<float>params1(1,m_pars,1,n_pars);

   // Read parameters values from file //
   ifstream pred_pars;
   pred_pars.open("output/pred_pars.tab");

   for(int i=1;i<=m_pars;i++)
   {
      for(int j=1;j<=n_pars;j++)
      {
         pred_pars >> params1(i,j);
      }
   }
   pred_pars.close();
   // -- //

   //for(int i=1;i<=val_lakes_index.UBound();i++) //
   //   cout << i <<" "<< val_lakes_index(i) << "\n";

   _vbc_vec<float> preds;
   preds = predict_p(params1,val_lakes_index,m_pars); //pars,indicies of validation set

   // Write predictions to file //
   ofstream pred_p_file;
   pred_p_file.open("output/pred_p.tab");
   cout << "---------------- " << "\n";
   for(int i=1;i<=n_val_lakes;i++)
      pred_p_file << val_lakes_index(i) << "\t";
   pred_p_file << "\n";
   for(int m=1;m<=m_pars;m++)
   {
      for(int i=1;i<=n_val_lakes;i++)
         pred_p_file << preds(m,i) << "\t";
      pred_p_file << "\n";
   }
   pred_p_file.close();
   // -- //

   // Write Validation lakes outcomes to file //
   ofstream val_d_file;
   val_d_file.open("output/val_lakes.dat");
   for(int i=1;i<=n_val_lakes;i++)
      val_d_file << lakes(val_lakes_index(i)).val_invaded << "\n";
   val_d_file.close();
   // -- //


   //TEST SAMPLING FUNCTION sample_w_replace()//
   /*
   _vbc_vec<int> test_vec(1,4);
   for(int i=1;i<=4;i++)
      test_vec(i) = i;

   _vbc_vec<int> rnd_vec;
   
   for(int i=1;i<=100;i++)
   {
      rnd_vec = sample_w_replace(test_vec);
      for(int j=1;j<=4;j++)
         cout << rnd_vec(j) << "\t";
      cout << "\n";
   }
   */
}

   return 0;
}
