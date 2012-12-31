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
    bool no_env=FALSE; //set to true for reproducing Gertzen.
    int n_sim_for_smooth=100; //number of sims for smoothing the likelihood surface
    bool ll=FALSE;
    int est_env=13;
    bool sim=FALSE;
    float d_par=2; //0.288344;   //to be fit
    float e_par=0.8;        //to be fit
    float c_par=1;          //to be fit
    float gamma_par=0;      //to be fit
    float glb_alpha=0.00001;          //to be fit (as a function of chem_vars)
    int n_lakes=1646;
    int n_sources=213;      //4496; with non grided Oi
    int n_val_lakes=88;    //Should be 100 -- need to review data parsing code.
    int n_chem_var=13;
    int from_year=1989;
    int to_year=2009;
    int n_sampled=0;
    _vbc_vec<float> d_matrix(1,n_sources,1,n_lakes);
    _vbc_vec<float> chem_pars(1,n_chem_var+1); //to be fit
    _vbc_vec<float> traf_mat(1,n_lakes,1,n_lakes);
    _vbc_vec<int> t_vec(1,n_lakes);
    _vbc_vec<int> sampled_index;
    _vbc_vec<int> val_lakes_index(1,n_val_lakes);

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
		float Wj,x,y,alpha;
		bool invaded;
        int discovered;
        int last_abs;
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
        _vbc_vec<float> Xit;
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
void calc_pp_switch(int,int,int);
void update_sim_pp(int);
void update_pp_l_hood(int, int);
int which_min(int,int);
int which_max(int,int);
void write_t();
void write_par();
void write_traf_mat();
void write_inv_stat();
void write_pp();

float average(_vbc_vec<float> );
float MLE_l_hood(_vbc_vec<float> *,_vbc_vec<float> *);

void likelihood_wrapperMCMC(_vbc_vec<float> *, float *,int );
void likelihood_wrapperMCMC_MD(_vbc_vec<float> *, float *,int );
float prior(float,int);
bool restrict_MCMC(float,int);

float prior_MD(_vbc_vec<float>,int);
bool restrict_MCMC_MD(_vbc_vec<float>);
_vbc_vec<float> predict_p(_vbc_vec<float>,_vbc_vec<int>,int); //pars,indicies

#endif;

////////////////////////////////////////////////////////////////////
void args(int argc,char *argv[])
{
	// decode arguments //
	int optind=1;
	while ((optind < argc) && (argv[optind][0]=='-')) 
	{
        string sw = argv[optind]; //-s flag for simulated data
      if (sw=="-s") 
		{    
			optind++;
			sim = true;
            n_sources = atoi(argv[optind]);
            optind++;
            n_lakes = atoi(argv[optind]);
	
            cout << "Using sims\n";
	    }
	   else if (sw=="-d")  //-d flag for real data
		{
			cout << "Using Data\n";
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
		   	<< "Options:\n\t-s: Simulated Data \n\t-d: real data\n" 
            << "Usage: ./gb -s <n_sources> <n_lakes>\n";
		 }
       optind++;
	} // end of arg-decoding
}

void inits()
{
   string tmp;
   ifstream init_file("inits.ini");

   init_file >> tmp;
   init_file >> run_type;
   init_file >> tmp;
   init_file >> post_length;
   init_file >> tmp;
   init_file >> no_env;
   init_file >> tmp;
   init_file >> n_sim_for_smooth;

   init_file.close();
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
    // Distance matrix //
    ifstream d_file;
    if(sim)
        d_file.open("sims/distance_matrix.csv"); //distance_matrix.csv        
    else
        d_file.open("../2010_bytho_data/distance_matrix_grd.csv"); //distance_matrix.csv
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
        o_file.open("sims/Oi.csv"); //Oi.csv (non-gridded)          
    else
        o_file.open("../2010_bytho_data/Oi_grd.csv"); //Oi.csv (non-gridded)
    for(int i = 1;i<=n_sources;i++)
	{
         o_file >> sources(i).Oi;
         sources(i).Oi=sources(i).Oi;     
    }
    o_file.close();

    // Lakes and chemistry //
    ifstream l_file;
    if(sim)
        l_file.open("sims/simmed_lakes.csv");
    else    
        l_file.open("../2010_bytho_data/lakes_processed.csv");
    /* File structure:
    "Hectares",
    "UTM_X",
    "UTM_Y",
    "invaded",
    "B_DISC_YEAR",
    "LAST_ABS",
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

    _vbc_vec<int> tmp_index_sampled(1,n_lakes);
    int n_validation = 0;
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

        // Catch for to_year < 2010.
        // modify lakes(i) to reflect missing information.
         if(lakes(i).discovered > 2009)
         {
            lakes(i).invaded = 0;
            lakes(i).discovered=0;
            n_validation++;
            val_lakes_index(n_validation)=i;
         }
         if(lakes(i).last_abs > 2009)
         {
            lakes(i).last_abs = 0;
            n_validation++;
            val_lakes_index(n_validation)=i;
         }



        if( (lakes(i).invaded == 1 || lakes(i).last_abs != 0) && lakes(i).discovered != from_year) //sampled lakes (excluding the seed(s))
        {
            n_sampled++;
            tmp_index_sampled(n_sampled)=i;
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

    l_file.close();

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
    #pragma omp parallel for
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
            sources(i).Gij(j) = sources(i).Gij(j) / sources(i).Ai;
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

    #pragma omp parallel for
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
}

// *************************************************

void calc_pp_switch(int l, int before, int after) //calculate the pp by only accounting for the change of a single lake status. STATE SPACE MODEL
{
    before=which_min(before,to_year);
    after=which_min(after,to_year);

    int lake_to_index,lake_from_index;
    if(before > after) //moved back in time (added propagules to the system)
    {
        for(int t = after; t<=before-1; t++)
        {
            for(int j=1;j<=state(t+1).n_calc_pp;j++) //n_calc_pp=the number of uninvaded+newly invaded sites
            {
                lake_to_index=state(t+1).calc_pp_index(j);
                lakes(lake_to_index).pp(t+1) += traf_mat(l,lake_to_index);
            }
        }
    }else if(after > before) //moved forward in time (subtracted propagules from the system)
    {
        for(int t = before; t<=after-1; t++)
        {
            //subtract propagules coming from l
            for(int j=1;j<=state(t+1).n_calc_pp;j++) //n_calc_pp=the number of uninvaded+newly invaded sites
            {
                lake_to_index=state(t+1).calc_pp_index(j);
                lakes(lake_to_index).pp(t+1) -= traf_mat(l,lake_to_index);
            }
            //calc pp to l for the years before+1 to after
            for(int i=1;i<=state(t).n_new_inv;i++)
            {
                lake_from_index=state(t).new_inv(i);
                lakes(l).pp(t+1) += traf_mat(lake_from_index,l);                
            }
        }
    }   //stayed the same
        //do nothing
}

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
void calc_pp_slow()
{
    for(int t=from_year; t<=to_year; t++)
    {
        for(int s=1;s<=n_sources;s++)
        {
            sources(s).Xit(t)=sources(s).Xit(t-1);
            for(int i=1;i<=state(t).n_new_inv;i++)
                sources(s).Xit(t) += sources(s).Gij(state(t).new_inv(i));
        }
    }
    for(int t=from_year+1;t<=to_year;t++)
    {
        for(int i=1;i<=state(t).n_calc_pp;i++) //n_calc_pp=the number of uninvaded+newly invaded sites
        {
            int lake_index=state(t).calc_pp_index(i);
            lakes(lake_index).pp(t)=0;
            for(int s=1;s<=n_sources;s++)
                lakes(lake_index).pp(t) =+ sources(s).Xit(t-1) * sources(s).Oi * sources(s).Gij(lake_index);
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
    float alpha;
    int lake_index;
    for(int t=from_year+1;t<=to_year;t++)
    {
        state(t).n_inv=state(t-1).n_inv;
        state(t).n_u_inv=0;
        state(t).n_new_inv=0;
        for(int i=1;i<=state(t-1).n_u_inv;i++)
        {
            lake_index=state(t-1).u_inv(i);

             alpha=calc_alpha(lake_index);

            if(  ( 1-exp(- pow(alpha *lakes(lake_index).pp(t), c_par ) )  >= runif(0,1) || lakes(lake_index).discovered==t ) 
                    && lakes(lake_index).last_abs < t ) // invade stochastically and restrict to observed pattern
            {
                state(t).n_inv++;
                state(t).n_new_inv++;
                state(t).new_inv(state(t).n_new_inv)=lake_index;
                t_vec(lake_index)=t;
                
            }else{
                state(t).n_u_inv++;
                state(t).u_inv(state(t).n_u_inv)=lake_index;
                t_vec(lake_index)=to_year+1;
            }
        }
        if(t<to_year)
            update_sim_pp(t);
    }
}
void update_sim_pp(int t)
{
    //add pp from each newly invaded lake at time t to each uninvaded lake
    int to_lake_index;
    int from_lake_index;
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
//         cout << "A\t"<<lake_index << "\t" <<  pow(lakes(lake_index).pp(t),c_par) * -alpha << "\n";
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
//         cout << "B\t"<< lake_index << "\t" <<  log(1-exp(tmp_lh)) << "\n";
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
//         cout << "C\t"<< lake_index << "\t" <<  log(1-exp(tmp_lh)) << "\n";
        }
    }

    //cout << lh << "\n";
    return lh;
}
float calc_alpha(int i)
{
   if(!no_env)
   {
       float z=chem_pars(1); //intercept
       if(est_env > 0)
       {
           for(int ch=2;ch<=est_env+1;ch++)
           {
               z += chem_pars(ch)*lakes(i).chem(ch-1);
           }
       }
       float alpha= -log(1-(1/(1+exp(-z)))); //see notebook for derivation of functional form
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
   d_par=params(1);
   e_par= params(2);   
   c_par=params(3);
   gamma_par=params(4);
   if(!no_env)
   {
      for(int i=1;i<=n_chem_var+1;i++)
         chem_pars(i)=params(4+i);
   }
   else
      glb_alpha=params(5);

   calc_traf();
   calc_traf_mat();
   calc_pp();

   sim_spread();
   sim_spread();
   sim_spread();

   for(int i=1;i<=n_sim_for_smooth;i++)
   {
        sim_spread();
        tmplhood(i) = l_hood();
        cout << i << "\t" << tmplhood(i) <<"\n";
   }

   float qll=average(tmplhood);
   
   for(int i=1;i<=params.UBound();i++)
      cout<< params(i) <<"\t";

   cout << qll <<"\n";
   if(std::isnan(qll))
      return(10000000);

   return(-qll);//-tve because simplex is a minimizer
}


/// MCMC LIB ////
void likelihood_wrapperMCMC(_vbc_vec<float> * params, float * l,int dim)
{
	*l = - MLE_l_hood(params, params);
}
void likelihood_wrapperMCMC_MD(_vbc_vec<float> * pars, float * l,int dim)
{
   _vbc_vec<float> params = *pars;
   float llmd;
   d_par=params(1);
   e_par= params(2);
   c_par=params(3);

   if(no_env)
      glb_alpha=params(4);
   else
   { 
      for(int i=1;i<=n_chem_var+1;i++)
         chem_pars(i)=params(3+i);
   }

   calc_traf();
   calc_traf_mat();
   calc_pp();

   for(int i =1;i<=3;i++)
      sim_spread();

   llmd=l_hood();

   int n_par=params.UBound();

   for(int i=1;i<=n_par;i++)
      cout<< params(i) <<"\t";

   cout << llmd <<"\n";

   *l = llmd;
}
bool restrict_MCMC(float param,int dim)
{
   //par order: d, e, c, alpha
   // or: d, e, c, B_0:B_13
   if(dim==1)
   {
      if(param <=0 || param >= 4)
         return TRUE;
   }
   if(dim==2)
   {
      if(param <=0 || param >=2)
         return TRUE;
   }
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
   if(no_env)
      return -log(x(4)); //prior on alpha
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

      if(param(3) <=0 )
         return TRUE;

      if(no_env && param(4) <= 0)
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
      if(!no_env)
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
   _vbc_vec<float> pred_p(1,m_pars,1,indicies.UBound());
   for(int m=1;m<=m_pars;m++)
   {
      d_par=params(m,1);
      e_par= params(m,2);   
      c_par=params(m,3);
      gamma_par=params(m,4);
      glb_alpha=params(m,5);

      calc_traf();
      calc_traf_mat();
      calc_pp();

      sim_spread();
      sim_spread();
      
      for(int i=1;i<=indicies.UBound();i++)
      {
         pred_p(m,i) = 1 - exp(-pow(glb_alpha*lakes(indicies(i)).pp(2010)+gamma_par,c_par));
      }
   }
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
    ofstream traf("output/traf_mat.dat");
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
   ofstream inv_stat("output/inv_stat.dat");
   for(int i=1;i<=n_lakes;i++)
   {
      inv_stat << t_vec(i) << "\n";
   }

   inv_stat.close();
}
void write_pp()
{
    ofstream pp_file("output/pp.dat");
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
