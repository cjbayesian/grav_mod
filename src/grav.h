#include<iostream>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<cmath>
#include<string>
#include<sstream>
#include<modRank.h>
#include<vbc_vector.h>
#include<Simplex.h>
#include<metro_hastingsMD.h>
#include<iomanip>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#ifndef GRAV_H
#define GRAV_H

using namespace std;
using namespace modRank;
using namespace vbc_lib;
using namespace simplex;
using namespace mcmcMD;

// GLOBALS //
    int run_type=2; //1:MLE, 2:MCMC, 3:sim spread from posterior
    int post_length=0;
    bool env=FALSE; //set to FALSE for reproducing Gertzen.
    bool fit_pdet=FALSE;
    int n_sim_for_smooth=100; //number of sims for smoothing the likelihood surface
    bool ll=FALSE;
    int est_env=13;
    bool sim=FALSE;
    bool gridded=TRUE;
    bool boot=TRUE;
    float pdet=0.2; //detection probability to be fit
    float d_par=2; //0.288344;   //to be fit
    float e_par=1;        //to be fit
    float c_par=1;          //to be fit
    float gamma_par=0;      //to be fit
    float glb_alpha=0.00001;          //to be fit (as a function of chem_vars)
    int n_lakes=1646;
    int n_sources=213;      //4496; with non grided Oi
    int n_val_lakes;    //Should be 100 -- need to review data parsing code.
    int n_chem_var=13;
    int from_year=1989;
    int to_year=2009;
    int n_sampled=0;
    float fixed_d=0;
    bool new_tr_par=TRUE;
    string parseeds_file;
    _vbc_vec<float> d_matrix;
    _vbc_vec<float> chem_pars(1,n_chem_var+1); //to be fit
    _vbc_vec<float> traf_mat;
    _vbc_vec<int> t_vec;
    _vbc_vec<int> sampled_index;
    _vbc_vec<int> val_lakes_index(1,n_val_lakes);
    _vbc_vec<int> val_lakes_inv_index(1,n_lakes);
    _vbc_vec<string> parnames;
//debugging output //
    ofstream track_lake("output/tl.dat");
    ofstream t_file("output/t_mcmc.dat");
    ofstream l_file("output/l_mcmc.dat");
    ofstream par_file("output/par_mcmc.dat");
    ofstream alpha_mcmc_file("output/alpha_mcmc.dat");
// CLASSES //
class clsLake
	{
	public:
		float Wj,x,y,alpha,Uj;
		bool invaded;
      bool val_invaded;
      int discovered;
      int last_abs;
      int status_2010;
      _vbc_vec<float> chem;
		_vbc_vec<float> Uij;
		_vbc_vec<float> Qjt;
      _vbc_vec<float> pp;
	};
class clsSource
	{
	public:
		float Ai;
		int Oi;
        _vbc_vec<float> Xit; //prop_invaded
        _vbc_vec<float> Dij;
		_vbc_vec<float> Gij;
	};

class clsState
	{
	public:
        int year;
        int n_inv;
        int n_new_inv;
        int n_u_inv;
        int n_calc_pp;
        _vbc_vec<int> inv;
        _vbc_vec<int> new_inv;
		_vbc_vec<int> u_inv;
		_vbc_vec<int> calc_pp_index;    
        
	};

_vbc_vec<clsLake> lakes(1,n_lakes);
_vbc_vec<clsSource> sources(1,n_sources);
_vbc_vec<clsState> state(from_year,to_year);


void args(int ,char *);
void get_gen_pars();
void sample_t_l(int);

float calc_alpha(int);
void calc_pp();
void calc_pp_validation(_vbc_vec<int>);
void update_sim_pp(int);
void update_pp_l_hood(int, int);
int which_min(int,int);
int which_max(int,int);
void write_t();
void write_par();
void write_traf_mat();
void calc_traf_mat_part(_vbc_vec<int>);
void clear_traf_mat();
void write_inv_stat();
void write_pp();
float l_hood();
float l_hood_detp();

float average(_vbc_vec<float> );
int wc_l(string);
float MLE_l_hood(_vbc_vec<float> *,_vbc_vec<float> *);

void likelihood_wrapperMCMC(_vbc_vec<float> *, float *,int );
void likelihood_wrapperMCMC_MD(_vbc_vec<float> *, float *,int );
float prior(float,int);
bool restrict_MCMC(float,int);

float prior_MD(_vbc_vec<float>,int);
bool restrict_MCMC_MD(_vbc_vec<float>);
_vbc_vec<float> predict_p(_vbc_vec<float>,_vbc_vec<int>,int); //pars,indicies
_vbc_vec<int> sample_w_replace(_vbc_vec<int>);
_vbc_vec<int> sample_wo_replace(_vbc_vec<int>,int);

void link_pars(_vbc_vec<string> ,_vbc_vec<float> );
void parse_par_seeds(string,_vbc_vec<float> *);
#endif;

////////////////////////////////////////////////////////////////////
void args(int argc,char *argv[])
{
	// decode arguments //
	int optind=1;
	while ((optind < argc) && (argv[optind][0]=='-')) 
	{
      string sw = argv[optind]; //-s flag for simulated data
      if (sw=="-seeds") 
		{    
			optind++;
            parseeds_file = argv[optind];
	   }
      else if (sw=="-runtype")
      {
         optind++;
         run_type = atoi(argv[optind]);
      }
      else if (sw=="-ll") //-ll flag to call l_hood() with given pars
      {
            ll=TRUE; 
            sim = true;  
            optind++;         
            n_sources = atoi(argv[optind]);
            optind++;
            n_lakes = atoi(argv[optind]);
            optind++;
            d_par = atof(argv[optind]);
            optind++;
            e_par = atof(argv[optind]);
            get_gen_pars();
        }
        else if (sw=="-lld") //-ll flag to call l_hood() with given pars
        {    
            ll=TRUE; 
            sim = false;  
            optind++;         
            d_par = atof(argv[optind]);
            optind++;
            e_par = atof(argv[optind]);
        }
        else
		{
         cerr << "Unknown switch: " << argv[optind] << endl 
		   	<< "Options:\n\t-seeds <seed params file>  \n\t-runtype: <integer run type>\n" 
            << "Usage: ./gb -s <n_sources> <n_lakes>\n";
            exit(1);
		}
        optind++;
	} // end of arg-decoding
}

void inits()
{
   string tmp;
   ifstream init_file("inits.ini");

   init_file >> tmp;
   init_file >> post_length;
   init_file >> tmp;
   init_file >> n_sim_for_smooth;
   init_file >> tmp;
   init_file >> gridded;
   init_file >> tmp;
   init_file >> boot;
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



void get_gen_pars()
{
    ifstream par_file;
    par_file.open("sims/ch_params.csv");
    for(int ch=1;ch<=n_chem_var+1;ch++)
        par_file >> chem_pars(ch);
    //par_file >> d_par;
    //par_file >> e_par;
    par_file.close();
}
void read_data()
{
    // Lakes and chemistry //
    ifstream l_file;
    if(sim){
        l_file.open("sims/simmed_lakes.csv");
        n_lakes = wc_l("sims/simmed_lakes.csv");
    }
    else{    
        //l_file.open("../2010_bytho_data/lakes_processed.csv"); //lakes_processed_erin_fix.csv");
        //n_lakes = wc_l("../2010_bytho_data/lakes_processed.csv");
        l_file.open("../2010_bytho_data/lakes_processed_normalized.csv"); 
        n_lakes = wc_l("../2010_bytho_data/lakes_processed_normalized.csv");
    }
    /* File structure:
    "Hectares",
    "UTM_X",
    "UTM_Y",
    "invaded",
    "B_DISC_YEAR",
    "LAST_ABS",
    "status_2010",
    "NAUT",
    "KKUT",
    "MGUT",
    "CAUT",
    "PPUT1",
    "SIO3UR",
    "DOC",
    "COLTR",
    "ALKTI",
    "ALKT",
    "PH",
    "COND25",
    "SECCHI.DEPTH"
    */
    traf_mat.redim(1,n_lakes,1,n_lakes);
    t_vec.redim(1,n_lakes);
    _vbc_vec<int> tmp_index_sampled(1,n_lakes);
    _vbc_vec<int> tmp_val_lakes_index(1,n_lakes);    
    _vbc_vec<int> tmp_val_inv_lakes_index(1,n_lakes); 
    int n_validation = 0;
    int n_invaded_in_validation_set = 0;
    for(int i = 1;i<=n_lakes;i++)
	{
        lakes(i).chem.redim(1,n_chem_var);
        lakes(i).pp.redim(from_year,to_year);
        lakes(i).pp(from_year)=0; //in li_hood() lakes invaded in 1989 don't contribute to likelihood.
        l_file >> lakes(i).Wj;
        l_file >> lakes(i).x;
        l_file >> lakes(i).y;
        l_file >> lakes(i).invaded;
        l_file >> lakes(i).discovered;
        l_file >> lakes(i).last_abs;
        l_file >> lakes(i).status_2010;

        // Catch for to_year < 2010.
        // modify lakes(i) to reflect missing information.
        /*
         lakes(i).val_invaded=0;
         if(lakes(i).discovered > 2009)
         {
            lakes(i).invaded=0;
            lakes(i).val_invaded=1;
            lakes(i).discovered=0;
            n_validation++;
            tmp_val_lakes_index(n_validation)=i;
            n_invaded_in_validation_set++;
            tmp_val_inv_lakes_index(n_invaded_in_validation_set)=i;
         }
         if(lakes(i).last_abs > 2009)
         {
            lakes(i).last_abs = 0;
            n_validation++;
            tmp_val_lakes_index(n_validation)=i;
         }

         */
        // Count validation lakes //
        if(lakes(i).status_2010 != 0)
        { 
           n_validation++;
           tmp_val_lakes_index(n_validation)=i;
           if(lakes(i).status_2010 == 2)
           {
              lakes(i).val_invaded = 1;
              n_invaded_in_validation_set++;
              tmp_val_inv_lakes_index(n_invaded_in_validation_set)=i;
           }
        }
        if( (lakes(i).discovered != 0 || lakes(i).last_abs != 0) && lakes(i).discovered != from_year) //sampled lakes (excluding the seed(s) and validation lakes)
        {
            if( (fit_pdet && lakes(i).discovered != 0) ||  !fit_pdet )
            {
                n_sampled++;
                tmp_index_sampled(n_sampled)=i;
            }
        }
        for(int j = 1;j<=n_chem_var;j++)
        {
            l_file >> lakes(i).chem(j);
        }
    }
    cout << "Number of lakes in validation set = "<< n_validation << "\n";

    sampled_index.redim(1,n_sampled);
    for(int i = 1;i<=n_sampled;i++)
	    sampled_index(i)=tmp_index_sampled(i);

    n_val_lakes = n_validation;
    val_lakes_index.redim(1,n_validation);
    for(int i = 1;i<=n_validation;i++)
      val_lakes_index(i)=tmp_val_lakes_index(i);

    val_lakes_inv_index.redim(1,n_invaded_in_validation_set);
    for(int i = 1;i<=n_invaded_in_validation_set;i++)
      val_lakes_inv_index(i)=tmp_val_inv_lakes_index(i);

    l_file.close();

    // Distance matrix //
    ifstream d_file;
    if(sim){
        d_file.open("sims/distance_matrix.csv"); //distance_matrix.csv
        n_sources = wc_l("sims/distance_matrix.csv");        
    }
    else
    {
      if(gridded){
        d_file.open("../2010_bytho_data/distance_matrix_grd.csv"); //distance_matrix.csv
        n_sources = wc_l("../2010_bytho_data/distance_matrix_grd.csv");
      }
      else{
         d_file.open("../2010_bytho_data/cj_roadconn.tab"); //distance_matrix.csv
         n_sources = wc_l("../2010_bytho_data/cj_roadconn.tab");
      }
    }

    d_matrix.redim(1,n_sources,1,n_lakes);
    sources.redim(1,n_sources);

    for(int i = 1;i<=n_sources;i++)
	{
        sources(i).Dij.redim(1,n_lakes);
        sources(i).Gij.redim(1,n_lakes);
        sources(i).Xit.redim(from_year-1,to_year);
        sources(i).Xit(from_year-1)=0;
        for(int j = 1;j<=n_lakes;j++)
        {
            d_file >> sources(i).Dij(j);
        }
    }
    d_file.close();

    // Oi //
    ifstream o_file;
    if(sim)
        o_file.open("sims/Oi.csv");        
    else
    {
         if(gridded)
            o_file.open("../2010_bytho_data/Oi_grd.csv"); //gridded 
         else
            o_file.open("../2010_bytho_data/Oi.csv"); // (non-gridded)
    }
    for(int i = 1;i<=n_sources;i++)
	{
         o_file >> sources(i).Oi;
         sources(i).Oi=sources(i).Oi;   
         //cout <<   sources(i).Oi << "\n";
    }
    o_file.close();

    //init_chem_pars
    //if(!ll)
    //{
        for(int ch=1;ch<=n_chem_var+1;ch++)
            chem_pars(ch)=0; //chem_pars(ch)=runif(-0.01,0.01);
        //chem_pars(1)=-5;
        if(sim)
            chem_pars(1)=0;
    //}

   // Static pp test //
    if(FALSE)
    {
        ifstream pp_file("sims/pp.dat");
        for(int i = 1;i<=n_lakes;i++)
	    {
            for(int year=from_year;year<=to_year;year++)
                pp_file >> lakes(i).pp(year);
        }
        pp_file.close();
    }
}
void init_state()
{
    state.redim(from_year,to_year);
    for(int t=from_year;t<=to_year;t++)
    {
        state(t).inv.redim(1,n_lakes);
        state(t).u_inv.redim(1,n_lakes);
        state(t).new_inv.redim(1,n_lakes);
        state(t).calc_pp_index.redim(1,n_lakes);
    }
}

//////////////////////////////////////////////////////////////
/////////////////// FUNCTIONS ////////////////////////////////
float grav_func(int i, int k)
{
    float tmp;
    tmp = pow(lakes(k).Wj,e_par) * pow(sources(i).Dij(k),-d_par);
    return tmp;
}
void calc_traf()
{
    //#pragma omp parallel for
    for(int i=1;i<=n_sources;i++)
    {
        sources(i).Ai=0;
        for(int k=1;k<=n_lakes;k++)
        {   
            sources(i).Gij(k) = grav_func(i,k);
            sources(i).Ai += sources(i).Gij(k);
        }
        for(int j=1;j<=n_lakes;j++)
        {
            sources(i).Gij(j) = sources(i).Gij(j) / (sources(i).Ai);
            if(std::isnan(sources(i).Gij(j)))
                sources(i).Gij(j) = 0;
        }
    }
}
void calc_traf_mat()
{
   //This is the bottleneck of the whole program.
   //(n_lakes*n_lakes - n_lakes)/2 * n_sources 
   // ~ 460million calcs every call.
   
   cout << "Calculating traf mat...\n";
   //#pragma omp parallel for
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j=1;j<=i-1;j++) //start at 0 since the diagonal of traf_mat is not needed.
        {
            traf_mat(i,j)=0;            
            for(int s=1;s<=n_sources;s++)
                traf_mat(i,j)+=sources(s).Gij(i)*sources(s).Gij(j)*sources(s).Oi;

            traf_mat(j,i)=traf_mat(i,j);
        }
    }
   cout << "Finished calculating traf mat...\n";
}
void calc_traf_mat_part(_vbc_vec<int> tracked)
{
   //Calculate rows of the traf_mat only as needed.
    for(int i=1;i<=n_lakes;i++)
    {
        if(tracked(i) == 1)
        {
            for(int j=1;j<=i-1;j++)
            {
                if(tracked(j) == 1)
                {
                    if(new_tr_par || traf_mat(i,j)==0) //only calc if new pars or not already calc'd.
                    {
                        traf_mat(i,j)=0;            
                        for(int s=1;s<=n_sources;s++)
                            traf_mat(i,j)+=sources(s).Gij(i)*sources(s).Gij(j)*sources(s).Oi;

                        traf_mat(j,i)=traf_mat(i,j);
                    }
                    //cout << i << " " << j << ": " << traf_mat(j,i) << "\n";
                }
            }
        }
    }
}
void clear_traf_mat()
{
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j=1;j<=n_lakes;j++)
        {
            traf_mat(i,j) = 0;
            traf_mat(j,i) = 0;
        }
    }
}

// Uj is the maximum propagule pressure to lake j (ie when everything is invaded).
void calc_Uj()
{
    for(int j=1;j<=n_lakes;j++)
    {
        lakes(j).Uj = 0;
        for(int i=1;i<=n_sources;i++)
        {   
            lakes(j).Uj += sources(i).Gij(j) * sources(i).Oi;
        }
    }
}
// *************************************************

void calc_pp()
{
    int lake_to_index, lake_from_index;
    for(int t=from_year+1;t<=to_year;t++)
    {
        for(int j=1;j<=state(t).n_calc_pp;j++) //n_calc_pp=the number of uninvaded+newly invaded sites
        {
            lake_to_index=state(t).calc_pp_index(j);
            lakes(lake_to_index).pp(t)=lakes(lake_to_index).pp(t-1);
            for(int i=1;i<=state(t-1).n_new_inv;i++)
            {
                lake_from_index=state(t-1).new_inv(i);
                lakes(lake_to_index).pp(t) += traf_mat(lake_from_index,lake_to_index);
            }
        }
    }
}
void calc_pp_validation(_vbc_vec<int> indicies) //for calculated pp in each relevant year for the validation lakes.
{
   int lake_to_index, lake_from_index;
   for(int l=1;l<=indicies.UBound();l++)
   {
      lake_to_index = indicies(l);
      for(int t=from_year+1;t<=to_year;t++)
      {      
         lakes(lake_to_index).pp(t) = lakes(lake_to_index).pp(t-1);
         for(int i=1;i<=state(t-1).n_new_inv;i++)
         {
            lake_from_index=state(t-1).new_inv(i);
            lakes(lake_to_index).pp(t) += traf_mat(lake_from_index,lake_to_index);
         }
         //cout << lake_to_index << "\t" << lakes(lake_to_index).pp(t) << "\n";
      }
   }
}

void calc_state()
{
    for(int t=from_year;t<=to_year;t++)
    {
        state(t).n_inv=0;
        state(t).n_new_inv=0;
        state(t).n_u_inv=0;
        state(t).n_calc_pp=0;
        for(int i=1;i<=n_lakes;i++)
        {
            if(t_vec(i)<=t)
            {
                if(t_vec(i)==t)
                {
                    state(t).n_new_inv++;
                    state(t).new_inv(state(t).n_new_inv)=i;
                    state(t).n_calc_pp++;
                    state(t).calc_pp_index(state(t).n_calc_pp)=i;
                }
                state(t).n_inv++;
                state(t).inv(state(t).n_inv)=i;
            }
            else
            {
                state(t).n_u_inv++;
                state(t).u_inv(state(t).n_u_inv)=i;
                state(t).n_calc_pp++;
                state(t).calc_pp_index(state(t).n_calc_pp)=i;
            }
        }
            //cout << "year\t" << t<< "\t#uninv\t" << state(t).n_u_inv << "\n";
    }   
}
void init_t()
{
    for(int i=1;i<=n_lakes;i++)
    {
        if(lakes(i).last_abs != 0 && lakes(i).invaded == 0) //observed uninvaded
        {
            t_vec(i)=runif(lakes(i).last_abs+1,to_year+2);
        }
        else if(lakes(i).last_abs != 0 && lakes(i).invaded == 1) //observed uninvaded at last_abs and observed invaded at discovered
        {
            t_vec(i)=runif(lakes(i).last_abs+1,lakes(i).discovered+1);
        }
        else if(lakes(i).discovered==from_year)
        {
            t_vec(i)=from_year;
        }
        else if(lakes(i).invaded == 1) // never observed uninvaded and observed invaded at discovered
        {
             t_vec(i)=runif(from_year+1,lakes(i).discovered+1);
        }
        else
        {
            //t_vec(i)=runif(from_year+1,to_year+(to_year-from_year)); //allow for statespace in the future.
            //t_vec(i)=runif(from_year+1,to_year+2);
            t_vec(i)=2011;
        }
    }
}
void sample_t()
{
    int t_var=2; //+- for uniform proposal of the state variable
    for(int i=1;i<=n_lakes;i++)
    {
        if(lakes(i).last_abs != 0 && lakes(i).invaded == 0) //observed uninvaded
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1));
        }
        else if(lakes(i).last_abs != 0 && lakes(i).invaded == 1) //observed uninvaded at last_abs and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1) );
        }
        else if(lakes(i).invaded == 1) // never observed uninvaded and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(from_year,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1));
        }
        else
        {
            t_vec(i)=runif(which_max(from_year+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1) ); //allow for statespace in the future.      
        }
    }
}
void sample_t_l(int i)
{
        int t_var=4; //+- for uniform proposal of the state variable
        if(lakes(i).last_abs != 0 && lakes(i).invaded == 0) //observed uninvaded
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1));
        }
        else if(lakes(i).last_abs != 0 && lakes(i).invaded == 1) //observed uninvaded at last_abs and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1) );
        }
        else if(lakes(i).invaded == 1) // never observed uninvaded and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(from_year,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1) );
        }
        else
        {
            t_vec(i)=runif(which_max(from_year+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1) ); //allow for statespace in the future.      
        }
}
void sim_spread()
{
    _vbc_vec<float> alpha(1,n_lakes);
    int lake_index;
    // determine which lakes we need to fill in traf_mat for:
    _vbc_vec<float> unifs(1,n_lakes,from_year+1,to_year);
    _vbc_vec<int> test_lakes(1,n_lakes);
    for(int i=1;i<=n_lakes;i++)
    {
        float min_unifs = 1;
        test_lakes(i) = 0;
        alpha(i)=calc_alpha(i);
        for(int t=from_year+1;t<=to_year;t++)
        {
            unifs(i,t) = runif(0,1);
            if( unifs(i,t) < min_unifs ) min_unifs = unifs(i,t);
        }
        if( 1-exp(- pow(alpha(i) * lakes(i).Uj + gamma_par, c_par ) ) > min_unifs || lakes(i).discovered != 0 || lakes(i).last_abs != 0 )
        {
            //cout << lakes(i).last_abs << "++++\n";
            test_lakes(i) = 1;
            //cout << i << "\n";
        }
    }
    calc_traf_mat_part(test_lakes);
            //cout << "------------------------------------\n";
    update_sim_pp(from_year);
    for(int t=from_year+1;t<=to_year;t++)
    {
        state(t).n_inv=state(t-1).n_inv;
        state(t).n_u_inv=0;
        state(t).n_new_inv=0;
        for(int i=1;i<=state(t-1).n_u_inv;i++)
        {
            lake_index=state(t-1).u_inv(i);
            if(test_lakes(lake_index) == 1 || lakes(lake_index).discovered==t)
            {
                //cout << "*********" << lakes(lake_index).pp(t) << "\n";
                //alpha(lake_index)=calc_alpha(lake_index);
                if(  ( 1-exp(- pow(alpha(lake_index) *lakes(lake_index).pp(t) + gamma_par, c_par ) )  >= unifs(lake_index,t) 
                         || lakes(lake_index).discovered==t ) 
                         && lakes(lake_index).last_abs <= t ) // invade stochastically and restrict to observed pattern
                {
                    state(t).n_inv++;
                    state(t).n_new_inv++;
                    state(t).new_inv(state(t).n_new_inv)=lake_index;
                    state(t).inv(state(t).n_inv)=lake_index;
                    t_vec(lake_index)=t;
                    
                } else {
                    state(t).n_u_inv++;
                    state(t).u_inv(state(t).n_u_inv)=lake_index;
                    t_vec(lake_index)=to_year+1;
                }
            }else{
                state(t).n_u_inv++;
                state(t).u_inv(state(t).n_u_inv)=lake_index;
                t_vec(lake_index)=to_year+1;
            }
        }
        if(t<to_year)
            update_sim_pp(t);
      //cout << state(t).n_inv << " ";
    }
    cout << "\r" << to_year << "=====" << float(state(to_year).n_inv)/float(n_lakes);
}
void update_sim_pp(int t)
{
    //add pp from each newly invaded lake at time t to each uninvaded lake
    int to_lake_index;
    int from_lake_index;
    //could just do this over those lakes that we will either test or need in likelihood.
    for(int i=1;i<=state(t).n_u_inv;i++)
    {
        to_lake_index=state(t).u_inv(i);
        for(int j=1;j<=state(t).n_new_inv;j++)
        {
            from_lake_index=state(t).new_inv(j);
            lakes(to_lake_index).pp(t+1)=lakes(to_lake_index).pp(t)+traf_mat(from_lake_index,to_lake_index);
        }
    }
}
void update_pp_l_hood(int lake_index,int t) //called from l_hood to calc pp where it is missing from sim
{
    int to_lake_index=lake_index;
    int from_lake_index;
    lakes(to_lake_index).pp(t)=lakes(to_lake_index).pp(t-1);
    for(int j=1;j<=state(t-1).n_new_inv;j++)
    {
        from_lake_index=state(t-1).new_inv(j);
        if(from_lake_index != to_lake_index) //should not receive pp from itself!
            lakes(to_lake_index).pp(t)+=traf_mat(from_lake_index,to_lake_index);
    }
}


float l_hood()
{
    float lh(0); float tmp_lh;
    int lake_index;
    for(int i=1;i<=n_sampled;i++)
    {
        lake_index=sampled_index(i);
        float alpha=calc_alpha(lake_index); //insert this inside lake loop.    
        if(lakes(lake_index).last_abs != 0) // observed uninvaded
        {
            for(int t=from_year+1;t<=lakes(lake_index).last_abs;t++)
            {
                if(t_vec(lake_index)<t)
                    update_pp_l_hood(lake_index,t);
                lh+= - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); //Prob uninv to that year
            }
            if(lakes(lake_index).discovered != 0) // discovered invaded
            {
                tmp_lh=0;
                for(int t=lakes(lake_index).last_abs+1;t<=lakes(lake_index).discovered;t++)
                {
                   if(t_vec(lake_index)<t)
                       update_pp_l_hood(lake_index,t);

                    tmp_lh += - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par);
                }
                lh+= log(1-exp(tmp_lh));
            }
        }else{
            tmp_lh=0;
            for(int t=from_year+1;t<=lakes(lake_index).discovered;t++) //never observed uninvaded, but discovered at T
            {
               if(t_vec(lake_index)<t)
                  update_pp_l_hood(lake_index,t);

                tmp_lh += - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par);
            }
            lh += log(1-exp(tmp_lh));
        }
/*    if(lake_index ==134) //pp only getting calc'd in 1990? wtf?
    {
        cout << lake_index << ": " << lakes(lake_index).discovered << ":: " << lakes(lake_index).last_abs << ":: "<< lh << "\n";
        cout << lakes(lake_index).pp(lakes(lake_index).discovered) << " :: " << lakes(lake_index).pp(lakes(lake_index).last_abs-2) << "\n";
        cout << t_vec(134) << "***\n";
        for(int i=1;i<=state(2009).n_u_inv;i++)
            cout << state(2009).u_inv(i) << "---\n";
        for(int i=1;i<=state(2009).n_inv;i++)
            cout << state(2009).inv(i) << "+++\n";

    }
*/
    }

    return lh;
}

// Likelihood for the pressence only model //
float l_hood_detp()
{
    float lh(0); float tmp_lh;
    int lake_index;
    for(int i=1;i<=n_sampled;i++)
    {
        lake_index=sampled_index(i);
        if(lakes(lake_index).discovered != 0)
        {
//cout << lake_index << "\t" << lakes(lake_index).discovered << "\n";
        float alpha=calc_alpha(lake_index); //insert this inside lake loop.    
        /*if(lakes(lake_index).last_abs != 0) // observed uninvaded
        {
            for(int t=from_year+1;t<=lakes(lake_index).last_abs;t++)
            {
                if(t_vec(lake_index)<t)
                    update_pp_l_hood(lake_index,t);
                lh+= - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); //Prob uninv to that year
            }
            if(lakes(lake_index).discovered != 0) // discovered invaded
            {
                tmp_lh=0;
                for(int t=lakes(lake_index).last_abs+1;t<=lakes(lake_index).discovered;t++)
                {
                   if(t_vec(lake_index)<t)
                       update_pp_l_hood(lake_index,t);

                    tmp_lh += - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par);
                }
                lh+= log(1-exp(tmp_lh));
            }
        }else{*/
            tmp_lh=0;
            for(int t=from_year+1;t<=lakes(lake_index).discovered;t++) //never observed uninvaded, but discovered at T
            {
                // The strategy here is that we use (1-poe) for t>t_vec(lake_index), poe for t = t_vec(lake_index) 
                // and (1-pdet) from t=t_vec(lake_index) to t<lakes(lake_index).discovered, pdet for 
                // t<lakes(lake_index).discovered.
                // might want to sample pdet in log space.

               if(t < t_vec(lake_index))
               {
                    //cout << lake_index << "\n";
                    //update_pp_l_hood(lake_index,t);
                    tmp_lh += - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); //Prob uninv that year
               }
               if(t == t_vec(lake_index))
               {
                    tmp_lh += log(1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) )); //Prob inv IN that year
//cout << lakes(lake_index).pp(t) << "\t" << log(1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) )) << "\t"<< tmp_lh << "\t"<< alpha <<"\n";
                    tmp_lh += log(1 - pdet); //Prob of not detecting
               }
               if(t > t_vec(lake_index))
                    tmp_lh += log(1 - pdet); //Prob of not detecting

            }
                lh += tmp_lh; //+ log(pdet); //Prob of detecting in the year of detection.
        //}
          }
    }

    return lh;
}


float calc_alpha(int i)
{
   if(env)
   {
       float z=chem_pars(1); //intercept
       if(est_env > 0)
       {
           for(int ch=2;ch<=est_env+1;ch++)
           {
               z += chem_pars(ch)*lakes(i).chem(ch-1);
           }
       }
       float alpha = -log(1-(1/(1+exp(-z)))); //see notebook for derivation of functional form
       //float alpha = z;
       //cout << "lake "<<i << "\t" << alpha <<"\n";
       return alpha;
   }else{
      return glb_alpha;
   }
}


/// Wrapper for l_hood that sims spread n_sim times
/// and to smooth a quasi-likelihood function for fitting
/// MLE estimates.
float MLE_l_hood(_vbc_vec<float> * pars, _vbc_vec<float> * dat)
{
   _vbc_vec<float> tmplhood(1,n_sim_for_smooth);
   _vbc_vec<float> params = *pars;
   link_pars(parnames,(*pars));
   // Parameter bounds //
   if(d_par <=0 || e_par < 0 || c_par < 0 || gamma_par < 0 || glb_alpha < 0 )
      return(10000000);

   if(fixed_d == 0)
   {
      calc_traf();
      //calc_traf_mat();
      //calc_pp();
   }else{
      //*pars(1)=fixed_d;
      d_par=fixed_d;
   }
   new_tr_par=TRUE;
   clear_traf_mat();
   
   for(int i=1;i<=n_sim_for_smooth;i++)
   {
        sim_spread();
        tmplhood(i) = l_hood();
        //cout << i << "\t" << tmplhood(i) <<"\n";
        new_tr_par=FALSE;
   }

   float qll=average(tmplhood);
   for(int i=1;i<=params.UBound();i++)
      cout << setw(15) << params(i) <<"\t";

   cout << qll <<"\n";
   if(std::isnan(qll))
      return(10000000);

    // Print out distribution of alpha values at MLE
    /*
    ofstream alphas_file;
    alphas_file.open("output/alphas.tab",std::fstream::app);
    for(int i=1;i<=n_lakes;i++)
        alphas_file << calc_alpha(i) << "\t";
    alphas_file << "\n";
    alphas_file.close();
    */

   return(-qll);//-tve because simplex is a minimizer
}


/// MCMC LIB ////
void likelihood_wrapperMCMC(_vbc_vec<float> * params, float * l,int dim)
{
	*l = - MLE_l_hood(params, params);
}
void likelihood_wrapperMCMC_MD(_vbc_vec<float> * pars, float * l,int dim)
{

   float llmd;
   link_pars(parnames,(*pars));

   calc_traf();

   //calc_traf_mat();
   //calc_pp();

   //for(int i =1;i<=3;i++)
   //   sim_spread();

   //llmd=l_hood();
   _vbc_vec<float> tmplhood(1,n_sim_for_smooth);
   for(int i=1;i<=n_sim_for_smooth;i++)
   {
        sim_spread();
        if(fit_pdet)
            tmplhood(i) = l_hood_detp();
        else
            tmplhood(i) = l_hood();
        //cout << i << "\t" << tmplhood(i) <<"\n";
   }
   
   llmd = average(tmplhood);
   //if(env)
   //    llmd = llmd + dnorm(chem_pars(1),-5,5,true);
   
   int n_par=(*pars).UBound();

   for(int i=1;i<=n_par;i++)
      cout<< (*pars)(i) <<"\t";

   cout << llmd <<"\n";

   *l = llmd;
}
bool restrict_MCMC(float param,int dim)
{
   //par order: d, e, c, alpha
   // or: d, e, c, B_0:B_13
/*   if(dim==1)
   {
      if(param <=0 || param >= 4)
         return TRUE;
   }
   if(dim==2)
   {
      if(param <=0 || param >=2)
         return TRUE;
   }
*/
/*   if(dim==3)
   {
     if(param <=0)
         return TRUE;
   }
   if(dim==4)
   {
      if(param <=0)
         return TRUE;
   }
*/
   return FALSE;
}
float prior_MD(_vbc_vec<float> x, int dim)
{
   if(!env)
      return -log(glb_alpha); //prior on alpha
   else
	   return 0; //uninformative prior (log)
}
bool restrict_MCMC_MD(_vbc_vec<float> param)
{
   //par order: d, e, c, alpha
   // or: d, e, c, B_0:B_13
      if(param(1) <=0 || param(1) >= 4)
         return TRUE;

      if(param(2) <=0 || param(2) >= 2)
         return TRUE;

      if(param(3) < 0 )
         return TRUE;

      if(!env && param(4) <= 0)
         return TRUE;

   return FALSE;
}
float prior(float x, int dim)
{
	return 0; //uninformative prior (log)
}



///////////////// Posterior Spead Sim ///////////////
void sim_spread_posterior()
{
   ifstream post_file("output/thinned.mcmc");
   int index;
   for(int i=1;i<=post_length;i++)
   {
      post_file >> index;
      post_file >> d_par;
      post_file >> e_par;      
      post_file >> c_par;
      if(env)
      {
         for(int i=1;i<=n_chem_var+1;i++)
            post_file >> chem_pars(i);
      }
      calc_traf();
      calc_traf_mat();

      sim_spread();
      for(int i=1;i<=10;i++)
      {
         sim_spread();
         write_t();
      }
      cout << i << " of "<< post_length<< "\n";
   }

   post_file.close();
   t_file.close();
}


/////////////// Prediction probability vectors //////////////

_vbc_vec<float> predict_p(_vbc_vec<float> params,_vbc_vec<int> indicies,int m_pars) //pars,indicies
{
   _vbc_vec<float> pred_p(1,m_pars,1,n_val_lakes);
   ofstream val_sim_file;
   if(sim)
        val_sim_file.open("sims/gb_output/val_sim_props.tab");
   else
        val_sim_file.open("output/val_sim_props.tab");
   for(int m=1;m<=m_pars;m++)
   {
      cout << m << " of " << m_pars << "\n";
      //for(int i=1;i<=n_pars,)
      //link_pars(parnames,(*parval));
      d_par = params(m,1);
      e_par = params(m,2);
      //e_par = 1;
      c_par = params(m,3);
      //gamma_par = params(m,4);
      //glb_alpha = params(m,5);
      calc_traf();
      clear_traf_mat();
         //calc_traf_mat();
         //calc_pp();
      
      if(env)
      {
         for(int i=1;i<=n_chem_var+1;i++){
            chem_pars(i)=params(m,3+i);
            cout << chem_pars(i) << "\n";}
      }


      // -- Using the simulation method -- //
      // prob of inv is calculated as the proportion of times invaded via simulation // 
      _vbc_vec<float> prop_val_invaded(1,n_val_lakes);
      for(int i=1;i<=n_val_lakes;i++)
         prop_val_invaded(i) = 0;

      int n_sims = 1000;
      for(int s=1;s<=n_sims;s++)
      {
         
         sim_spread();
         new_tr_par = FALSE;
         for(int i=1;i<=n_val_lakes;i++)
         {
            if(t_vec(indicies(i)) != to_year+1)
               prop_val_invaded(i) += 1;
         }
      }
      cout << "\n";
      for(int i=1;i<=n_val_lakes;i++)
      {
            prop_val_invaded(i) = float(prop_val_invaded(i))/float(n_sims);
            val_sim_file << prop_val_invaded(i) << "\t";
            cout << prop_val_invaded(i) << "\t";
      }
      cout << "\n";
      val_sim_file << "\n";
      // --------------------------------------------- //


      calc_pp_validation(indicies); // fill in all pp values for validation lakes
                                    // since calc_pp() is optimized to only calc pp for uninvaded lakes.
      int cal_from;
      for(int i=1;i<=n_val_lakes;i++)
      {
         float alpha=calc_alpha(indicies(i)); //insert this inside lake loop.
         // 1 - prob of remaining uninvaded during the period since last observed absence //
         if(lakes(indicies(i)).last_abs != 0)
            cal_from = lakes(indicies(i)).last_abs;
         else
            cal_from = from_year;
   
         float log_p_uninv = 0;
         for(int t=cal_from+1;t<=2010;t++)
            log_p_uninv += -pow(alpha*lakes(indicies(i)).pp(t)+gamma_par,c_par);

          pred_p(m,i) = 1 - exp(log_p_uninv);
      }
      
   }
   val_sim_file.close();
   return(pred_p);
}





///////////////// Utilities ///////////////////

float average(_vbc_vec<float> x)
{
   int n=x.UBound();
   float tmp_sum=0;
   for(int i=1;i<=n;i++)
      tmp_sum += x(i);

   float nf=n;
   float avg= (1/nf)*tmp_sum;
   return(avg);
}
int which_min(int a,int b)
{
   if(a<=b)
       return a;
   else
       return b;
}
int which_max(int a,int b)
{
   if(a>=b)
       return a;
   else
       return b;
}
void write_t()
{
    for(int i=1;i<=n_lakes;i++)
        t_file<< t_vec(i) << "\t";

    t_file << "\n";
}
void write_par()
{
        for(int i=1;i<=n_chem_var+1;i++)
            par_file<< chem_pars(i) << "\t";

        par_file << d_par << "\t" << e_par;

    par_file << "\n";

    //for(int i=1;i<=n_lakes;i++)
    //    alpha_mcmc_file << calc_alpha(i) << "\t";

    //alpha_mcmc_file << "\n";
}
void write_traf_mat()
{
    ofstream traf;
    if(sim)
        traf.open("sims/gb_output/traf_mat.dat");    
    else
        traf.open("output/traf_mat.dat");
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j=1;j<=n_lakes;j++)
            traf << traf_mat(i,j) << ",";
        traf << "\n";
    }
    traf.close();
}
void write_inv_stat()
{
    ofstream inv_stat;
    if(sim)
       inv_stat.open("sims/gb_output/inv_stat.dat");
    else
       inv_stat.open("output/inv_stat.dat");
   for(int i=1;i<=n_lakes;i++)
   {
      inv_stat << t_vec(i) << "\n";
   }

   inv_stat.close();
}
void write_pp()
{
    ofstream pp_file;
    if(sim)
        pp_file.open("sims/gb_output/pp.dat");
    else
        pp_file.open("output/pp.dat");
    int lake_to_index;
    for(int t=from_year+1;t<=to_year;t++)
    {
        for(int j=1;j<=state(t).n_calc_pp;j++) //n_calc_pp=the number of uninvaded+newly invaded sites
        {
            lake_to_index=state(t).calc_pp_index(j);
            pp_file << t << "\t" << lake_to_index << "\t" << lakes(lake_to_index).pp(t)<<"\n";
        }
    }
   pp_file.close();
}

// Sample an integer vector with replacement. -- For generating bootstrap samples of lake indicies.
_vbc_vec<int> sample_w_replace(_vbc_vec<int> vec)
{
   _vbc_vec<int> s_vec(1,vec.UBound());
   int tmp_index;
   for(int i=1;i<=vec.UBound();i++)
   {
      tmp_index = (int) runif(1,vec.UBound()+1);
      s_vec(i) = vec(tmp_index);
   }

   return(s_vec);
}
// Sample an integer vector without replacement.
_vbc_vec<int> sample_wo_replace(_vbc_vec<int> vec,int n)
{
   int tmp_index,tmp;
   _vbc_vec<int> s_vec(1,n);
   //cout<< "In sample_wo_replace XX\n";
   for(int i=1;i<=n;i++)
   {
      tmp_index = (int) runif(1,vec.UBound() + 2 - i);
      tmp = vec(tmp_index);
      vec(tmp_index) = vec(vec.UBound() + 1 - i);
      vec(vec.UBound() + 1 - i) = tmp;
      s_vec(i)=tmp;
   }
   //cout<< "In sample_wo_replace YY\n";
   return(s_vec);
}

// Emulates wc -l (number of lines in a file)
int wc_l(string path_to_file)
{
   int number_of_lines=0;
   std::string line;
   std::ifstream myfile(path_to_file.c_str());
   while (std::getline(myfile, line))
        ++number_of_lines;
   myfile.close();
   return(number_of_lines);
}

// Parse params
// Read from a file defining the parameters to be used for fitting and simulating
void parse_par_seeds(string param_def_file,_vbc_vec<float> *parval)
{
    ifstream parseeds(param_def_file.c_str());
    int n_par = wc_l(param_def_file);
    (*parval).redim(1,n_par);
    parnames.redim(1,n_par); 
    for(int i=1;i<=n_par;i++)
    {
        parseeds >> parnames(i);
        parseeds >> (*parval)(i);
        cout <<  parnames(i) << "\t" << (*parval)(i) << "\n";        
    }
    link_pars(parnames,(*parval));
    parseeds.close();
    //Test
    //cout << glb_alpha << " :: " << d_par << " :: " << e_par << " :: " << c_par << "\n";
}


//Link par is called from parse_params
void link_pars(_vbc_vec<string> parname_str,_vbc_vec<float> parval)
{
    for(int i = 1; i<= parval.UBound(); i++)
    {
        if(parname_str(i) == "alpha")
            glb_alpha = parval(i);
        else if(parname_str(i) == "d")
            d_par = parval(i);
        else if(parname_str(i) == "e")
            e_par = parval(i);
        else if(parname_str(i) == "c")
            c_par = parval(i);
        else if(parname_str(i) == "gamma")
            gamma_par = parval(i);
        else if(parname_str(i) == "env0") //intercept
        {
            chem_pars(1) = parval(i);
            env = TRUE;
        }
        else if(parname_str(i) == "env1")
            chem_pars(2) = parval(i);
        else if(parname_str(i) == "env2")
            chem_pars(3) = parval(i);
        else if(parname_str(i) == "env3")
            chem_pars(4) = parval(i);
        else if(parname_str(i) == "env4")
            chem_pars(5) = parval(i);
        else if(parname_str(i) == "env5")
            chem_pars(6) = parval(i);
        else if(parname_str(i) == "env6")
            chem_pars(7) = parval(i);
        else if(parname_str(i) == "env7")
            chem_pars(8) = parval(i);
        else if(parname_str(i) == "env8")
            chem_pars(9) = parval(i);
        else if(parname_str(i) == "env9")
            chem_pars(10) = parval(i);
        else if(parname_str(i) == "env10")
            chem_pars(11) = parval(i);
        else if(parname_str(i) == "env11")
            chem_pars(12) = parval(i);
        else if(parname_str(i) == "env12")
            chem_pars(13) = parval(i);
        else if(parname_str(i) == "env13")
            chem_pars(14) = parval(i);
        else if(parname_str(i) == "env14")
            chem_pars(15) = parval(i);
        else if(parname_str(i) == "env15")
            chem_pars(16) = parval(i);
        else if(parname_str(i) == "env16")
            chem_pars(17) = parval(i);
        else if(parname_str(i) == "env17")
            chem_pars(18) = parval(i);
        else if(parname_str(i) == "env18")
            chem_pars(19) = parval(i);
        else if(parname_str(i) == "env19")
            chem_pars(20) = parval(i);
        else if(parname_str(i) == "env20")
            chem_pars(21) = parval(i);
        else
        {
            cout << "Malformed parameter name\n"; 
            exit(1);
        }
    }
}








////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  SIMS  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Globals for sim //
float glb_alpha_min, glb_alpha_max, c_par_min, c_par_max;
float pdet_min, pdet_max;
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


// Initiallize (seed) the largest lake //


// Simulate the invasion //


// Simulate the detection process //


// Fit the presence only model //




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
