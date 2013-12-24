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
    clear_traf_mat();

   ofstream ll_("output/ll_test.dat");
   if(FALSE)
   {
      calc_traf();
      calc_traf_mat();
      calc_pp(); //!!
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
//i       d        e           c     B_o       NAUT        KKUT      MGUT
//8394 1.06513 0.531416 0.000125572 -13.713 0.00405104 -0.00475339 0.0038181
//      CAUT     PPUT1  SIO3UR      DOC       COLTR     ALKTI         ALKT
//0.00403173 0.0571088 1.28896 -0.31447 -0.00654376 0.0621327 -0.000480646
//        PH     COND25 SECCHI.DEPTH       NA
//0.00852002 0.00335544    0.0995729 -380.192

//0.407766	1.44137	0.417034	-10.8476	0.708086 -0.0150081	-0.404532	-0.512538	0.689779	0.43177-0.584563	-0.224491	0.786624	1.0269	0.868935	-0.0982729	0.37267

   chem_pars(1)=-10;
   chem_pars(2)=0.7;
   chem_pars(3)=-0.01; //-2:1 MLE ~0
   chem_pars(4)=-0.38; //-1.5:1.5 MLE ~0
   chem_pars(5)=-0.5; //-0.4:0.2 MLE ~0
   chem_pars(6)=0.68; //0:0.1 MLE 0.05    ***
   chem_pars(7)=0.43; //-0.4:0.2 MLE ~1.3  **** 
   chem_pars(8)=-0.58; //-0.7:0.1 MLE ~-0.31 ****
   chem_pars(9)=-0.22; //-0.4:0.2 MLE ~-0.01 **
   chem_pars(10)=0.78; //-0.1:0.1 MLE ~0.06 *
   chem_pars(11)=1.02; //-0.16:0.1 MLE ~0
   chem_pars(12)=0.85; //-0.2:0.2 MLE ~0
   chem_pars(13)=-0.09; //-0.04:0.01 MLE ~0
   chem_pars(14)=0.37; //-0.3:0.3 MLE ~0.1

   int test_ch=14;
   d_par=1.54;
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
   // need a likelihood function wrapper to call l_hood() 
   // multiple times and average the result
   // to smooth out stochastic surface.
   // BOOTSTRAP RESAMPLING OF DATA (SAMPLED LAKES) TO GENERATE CI //

   float garbage=l_hood();
   int n_reps = 1000;
   ofstream par_file;

   int n_pars; 
   _vbc_vec<float> params1;
   _vbc_vec<float> dat1;
   _vbc_vec<float> MLE_params;

    if(sim)
        par_file.open("sims/gb_output/pred_pars.tab");
    else
        par_file.open("output/pred_pars.tab");    

    parse_par_seeds(parseeds_file,&params1);
    cout << env << " <<==\n";

    n_pars = params1.UBound();
    
    dat1.redim(1,n_pars);
    MLE_params = params1;

    _vbc_vec<int> tmp_index_sampled;
    tmp_index_sampled = sampled_index;

    if(boot)
    {
        ofstream boot_file;
        if(sim)
            boot_file.open("sims/gb_output/boot_lakes.tab");
        else
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
            gertzen_rep.start(&dat1,&params1, &MLE_l_hood,n_pars, 1e-2);
            gertzen_rep.getParams(&MLE_params);

            cout << "\n\nMLE "<< i << " of " << n_reps << "\n\n";
            for(int p=1;p<=n_pars;p++)
                par_file << MLE_params(p) <<"\t";
            par_file << "\n";
            par_file.flush();
            params1 = MLE_params; //seed the next run with the MLE from this one.
        }
        boot_file.close();
    }else{
        simplex::clsSimplex<float> gertzen_rep;
        //gertzen_rep.set_param_small(1e-3);
        gertzen_rep.start(&dat1,&params1, &MLE_l_hood,n_pars, 1e-2);
        gertzen_rep.getParams(&MLE_params);

        cout << "\n\nMLE\n";
        for(int p=1;p<=n_pars;p++)
         par_file << MLE_params(p) <<"\t";
        par_file << "\n";
        par_file.flush();

        // Print out distribution of alpha values at MLE
        ofstream alphas_file;
        alphas_file.open("output/alphas.tab",std::fstream::app);
        for(int i=1;i<=n_lakes;i++)
        {
            alphas_file << calc_alpha(i) << "\n";
        }
        alphas_file.close();
    }
    par_file.close();

}



if(run_type==2)
{
    //MCMC lib
   n_sim_for_smooth = 1; 
	string mcmc_file("output/lib.mcmc");

   _vbc_vec<float> params;
   _vbc_vec<float> prop_width;

   parse_par_seeds(parseeds_file,&params);
   cout << env << " <<==\n";

   int n_pars = params.UBound();
   prop_width.redim(1,n_pars);

   prop_width(1)=0.05; //d
   prop_width(2)=0.05; //e
   prop_width(3)=0.05; //c

    for(int i=4;i<=n_pars;i++)
        prop_width(i)=0.0001;

    _vbc_vec<float> prop_sigma;
    prop_sigma = diag(prop_width);

      mcmcMD::run_mcmc(params, 
         prop_sigma, 
         &likelihood_wrapperMCMC_MD,
         &prior_MD, 
         &restrict_MCMC_MD, 
         50000, 
         50, 
         1, 
         mcmc_file.c_str(),
         true,
         true,
         true,
         500,
         4);
}


/// Sim from posterior ///
if(run_type==3)
   sim_spread_posterior();

/// Traf tests ////
if(run_type==4)
{
   ofstream traf_ll_file("output/traf_ll.dat");

   d_par=1.54;
   e_par=2;
   c_par=0.8;
   calc_traf();
   calc_traf_mat();
   calc_pp();
   sim_spread();
   cout << l_hood() <<"\n";
   
    /*
    for(int i=1;i<=70;i++)
    {
        e_par=e_par+0.05;
        cout << e_par << "\t";
        calc_traf();
        calc_traf_mat();
        calc_pp();
        sim_spread();
        traf_ll_file << e_par  << "\t" << l_hood()  << "\n"; 
        cout << l_hood() <<"\t";
        sim_spread();
        traf_ll_file << e_par  << "\t" << l_hood()  << "\n"; 
        cout << l_hood() <<"\t";
        sim_spread();
        traf_ll_file << e_par  << "\t" << l_hood()  << "\n"; 
        cout << l_hood() <<"\n";
        traf_ll_file << e_par  << "\t" << l_hood()  << "\n"; 
    }


   for(int i=1;i<=10;i++)
   {
      e_par=e_par+0.2;
      c_par=0.8;
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
    */
   traf_ll_file.close();
   calc_pp();
   write_pp(); 
   write_inv_stat();
}

/// Holdout sets for internal AUC ////
if(run_type==5)
{  
   _vbc_vec<float>params1;
   parse_par_seeds(parseeds_file,&params1);
   cout << env << " <<==\n";

   int n_pars = params1.UBound();
   //...
}

/// Predictions from Bootstapped MLE ////
if(run_type==6)
{  

   int m_pars, n_pars;
   _vbc_vec<float> params1;
   ifstream pred_pars;
   if(!env)
   {
      //m reps of the bootstrap/posterior
      if(sim)
        m_pars = wc_l("sims/gb_output/pred_pars.tab");
      else
        m_pars = wc_l("output/pred_pars.tab");    
      n_pars = 4;
      cout << "# Generating a " << n_val_lakes << " by " << m_pars << " prediction matrix\n";

      params1.redim(1,m_pars,1,n_pars);

      // Read parameters values from file //
      if(sim)
          pred_pars.open("sims/gb_output/pred_pars.tab");
      else
          pred_pars.open("output/pred_pars.tab");

   }else{
      //m reps of the bootstrap/posterior
      if(sim)
          m_pars = wc_l("sims/gb_output/pred_parsENV.tab");
      else
          m_pars = wc_l("output/pred_parsENV.tab");
      n_pars=18; //13 env + intercept + d,e,c,gamma
      cout << "# Generating a " << n_val_lakes << " by " << m_pars << " prediction matrix\n";

      params1.redim(1,m_pars,1,n_pars);

      // Read parameters values from file //
      if(sim)
          pred_pars.open("sims/gb_output/pred_parsENV.tab");
      else
          pred_pars.open("output/pred_parsENV.tab");
   }
    // read in predicted parameters from bootstrap or Bayesian
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
   write_traf_mat();
   // Write predictions to file //
   ofstream pred_p_file;
   if(sim)
       pred_p_file.open("sims/gb_output/pred_p.tab");
   else
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
   if(sim)
       val_d_file.open("sims/gb_output/val_lakes.dat");
   else
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


/// Validation site Predictions from Bayesian Posterior ////
if(run_type==7)
{  

    int m_pars, n_pars;
    _vbc_vec<float> params1;
    ifstream pred_pars;

    m_pars = wc_l("output/mcmct.lib");    
    //m reps of the bootstrap/posterior
    n_pars=17; //13 env + intercept + d,e,c
    cout << "# Generating a " << n_val_lakes << " by " << m_pars << " prediction matrix\n";

    params1.redim(1,m_pars,1,n_pars);
    pred_pars.open("output/lib.mcmc500000t");
    // read in predicted parameters from bootstrap or Bayesian
    float jnk;
    for(int i=1;i<=m_pars;i++)
    {
        pred_pars >> jnk; //index
        for(int j=1;j<=n_pars;j++)
        {
            pred_pars >> params1(i,j);
        }
        pred_pars >> jnk; //log-likelihood
    }
    pred_pars.close();
    // -- //

    //for(int i=1;i<=val_lakes_index.UBound();i++) //
    //   cout << i <<" "<< val_lakes_index(i) << "\n";

    _vbc_vec<float> preds;
    preds = predict_p(params1,val_lakes_index,m_pars); //pars,indicies of validation set
    //write_traf_mat();
    // Write predictions to file //
    ofstream pred_p_file;
    
    pred_p_file.open("output/pred_p_bayes.tab");
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
   if(sim)
       val_d_file.open("sims/gb_output/val_lakes.dat");
   else
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
