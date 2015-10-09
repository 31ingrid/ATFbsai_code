#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <atf.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("atfbsai_is2014_np.mcmc.out");;
  styr.allocate("styr");
  endyr.allocate("endyr");
  styr_fut.allocate("styr_fut");
  endyr_fut.allocate("endyr_fut");
  nsurv.allocate("nsurv");
  median_rec.allocate("median_rec");
  first_age.allocate("first_age");
  last_age.allocate("last_age");
  nages.allocate("nages");
  first_length.allocate("first_length");
  last_length.allocate("last_length");
cout<<"nages"<<nages<<std::endl;
  nsurv_aged.allocate("nsurv_aged");
  phase_F40.allocate("phase_F40");
  phase_logistic_sel.allocate(1,2*nsurv,"phase_logistic_sel");
  phase_fishery_sel.allocate(1,2,"phase_fishery_sel");
  phase_alphabeta.allocate("phase_alphabeta");
  q_Phase.allocate(1,nsurv,"q_Phase");
  phase_selcoffs.allocate("phase_selcoffs");
  nselages.allocate("nselages");
  nselages_srv.allocate(1,nsurv,"nselages_srv");
  fishsel_LB_f.allocate(1,2,"fishsel_LB_f");
  fishsel_LB_m.allocate(1,2,"fishsel_LB_m");
  fishsel_UB_f.allocate(1,2,"fishsel_UB_f");
  fishsel_UB_m.allocate(1,2,"fishsel_UB_m");
  fishsel_prior_f.allocate(1,2,"fishsel_prior_f");
  fishsel_prior_m.allocate(1,2,"fishsel_prior_m");
  nsel_params.allocate(1,nsurv,"nsel_params");
cout<<"nsel_params"<<nsel_params<<std::endl;
  sel_prior_f.allocate(1,2*nsurv,"sel_prior_f");
  sel_prior_m.allocate(1,2*nsurv,"sel_prior_m");
cout<<"sel_prior_f"<<sel_prior_f<<std::endl;
  sel_LB_f.allocate(1,2*nsurv,"sel_LB_f");
  sel_LB_m.allocate(1,2*nsurv,"sel_LB_m");
  sel_UB_f.allocate(1,2*nsurv,"sel_UB_f");
  sel_UB_m.allocate(1,2*nsurv,"sel_UB_m");
cout<<"max(nsel_params)"<<max(nsel_params)<<std::endl;     
  sel1_desc_prior_f.allocate(1,2,"sel1_desc_prior_f");
  sel1_desc_prior_m.allocate(1,2,"sel1_desc_prior_m");
  sel1_desc_LB_f.allocate(1,2,"sel1_desc_LB_f");
  sel1_desc_LB_m.allocate(1,2,"sel1_desc_LB_m");
  sel1_desc_UB_f.allocate(1,2,"sel1_desc_UB_f");
  sel1_desc_UB_m.allocate(1,2,"sel1_desc_UB_m");
  nlen.allocate("nlen");
  nobs_fish.allocate("nobs_fish");
  yrs_fish.allocate(1,nobs_fish,"yrs_fish");
cout<<"yrs_fish"<<yrs_fish<<std::endl; 
  nsamples_fish.allocate(1,2,1,nobs_fish,"nsamples_fish");
  nobs_srv.allocate(1,nsurv,"nobs_srv");
  yrs_srv.allocate(1,nsurv,1,nobs_srv,"yrs_srv");
  nobs_srv_length.allocate(1,nsurv,"nobs_srv_length");
cout<<"nobs_srv_length"<<nobs_srv_length<<std::endl; 
  yrs_srv_length.allocate(1,nsurv,1,nobs_srv_length,"yrs_srv_length");
  nsamples_srv_length_fem.allocate(1,nsurv,1,nobs_srv_length,"nsamples_srv_length_fem");
  nsamples_srv_length_mal.allocate(1,nsurv,1,nobs_srv_length,"nsamples_srv_length_mal");
  nsamples_srv_length.allocate(1,nsurv,1,2,1,nobs_srv_length,"nsamples_srv_length");
cout<<"nsamples_srv_length"<<nsamples_srv_length<<std::endl;     
  obs_p_fish.allocate(1,2,1,nobs_fish,1,nlen,"obs_p_fish");
cout<<"obs_p_fish(1)"<<obs_p_fish(1)<<std::endl;  
  obs_p_srv_length_fem.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"obs_p_srv_length_fem");
  obs_p_srv_length_mal.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"obs_p_srv_length_mal");
  catch_bio.allocate(styr,endyr,"catch_bio");
  obs_srv.allocate(1,nsurv,1,nobs_srv,"obs_srv");
  obs_srv_sd.allocate(1,nsurv,1,nobs_srv,"obs_srv_sd");
  wt.allocate(1,2,1,nages,"wt");
cout<<"wt"<<wt<<std::endl;  
  maturity.allocate(1,nages,"maturity");
cout<<"maturity"<<maturity<<std::endl;
  lenage.allocate(1,2,1,nages,1,nlen,"lenage");
cout<<"lenage"<<lenage<<std::endl; 
 nyrs_temps = 33;//nobs_srv(1); change back for BSAI assessment
cout<<"nyrs_temps"<<nyrs_temps<<std::endl; 
  bottom_temps.allocate(1,nyrs_temps,"bottom_temps");
cout<<"nyrs_temps"<<nyrs_temps<<std::endl;
  monot_sel.allocate("monot_sel");
  wt_like.allocate(1,8,"wt_like");
cout<<"wt_like"<<wt_like<<std::endl;               
  nobs_srv_age.allocate(1,nsurv_aged,"nobs_srv_age");
  yrs_srv_age.allocate(1,nsurv_aged,1,nobs_srv_age,"yrs_srv_age");
  nsamples_srv_age.allocate(1,nsurv_aged,1,2,1,nobs_srv_age,"nsamples_srv_age");
  obs_p_srv_age_fem.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"obs_p_srv_age_fem");
  obs_p_srv_age_mal.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"obs_p_srv_age_mal");
  M.allocate(1,2,"M");
cout<<"M"<<M<<std::endl;            
  offset_const.allocate("offset_const");
  q_Lower_bound.allocate(1,nsurv,"q_Lower_bound");
  q_Upper_bound.allocate(1,nsurv,"q_Upper_bound");
  q_surv_prior_mean.allocate(1,nsurv,"q_surv_prior_mean");
  nparams_srv.allocate(1,nsurv,"nparams_srv");
cout<<"q_surv_prior_mean"<<q_surv_prior_mean<<std::endl;
  assess.allocate("assess");
  mean_log_rec_prior.allocate("mean_log_rec_prior");
  log_avg_fmort_prior.allocate("log_avg_fmort_prior");
  fish_sel50_f_bound.allocate("fish_sel50_f_bound");
  fish_slope_m_bound.allocate("fish_slope_m_bound");
  srv1_sel50_m_bound.allocate("srv1_sel50_m_bound");
  obs_p_fish2.allocate(1,2,1,nobs_fish,first_length,last_length);
  obs_p_srv_age_fem2.allocate(1,nsurv_aged,1,nobs_srv_age,first_age,last_age);
  obs_p_srv_age_mal2.allocate(1,nsurv_aged,1,nobs_srv_age,first_age,last_age);
  obs_p_srv_length_fem2.allocate(1,nsurv,1,nobs_srv_length,first_length,last_length);
  obs_p_srv_length_mal2.allocate(1,nsurv,1,nobs_srv_length,first_length,last_length);
  wt2.allocate(1,2,first_age,last_age);
  maturity2.allocate(first_age,last_age);
cout<<"sel_prior_f(1,1)"<<sel_prior_f(1,1)<<std::endl; 
cout<<"assess"<<assess<<std::endl;  
  cv_srv.allocate(1,nsurv,1,nobs_srv);
   styr_rec=styr-(last_age-first_age+1)+1;
   if(nselages>nages) 
   {nselages=nages;  
   cout<<"Warning selectivity: is set to be estimated on more ages than are in the model."<<std::endl;  }
   for (i=1; i<= nsurv; i++){
   if(nselages_srv(i)>nages) nselages_srv(i)=nages;
   }
   //calculate cv for surveys
   for (int j=1;j<=nsurv;j++)
   {
     for (i=1;i<=nobs_srv(j);i++)
     { 
   cv_srv(j,i)=obs_srv_sd(j,i)/(double)obs_srv(j,i); 
     }
   } 
   //change weights to tons  
   wt=wt*.001;               
   for (int i=first_age;i<=last_age;i++)
   {
	  for(int j=1;j<=2;j++)
	  {
      wt2(j,i)=wt(j,i); 
      }  
   maturity2(i)=maturity(i);
   }     
  cout<<"wt2 "<<wt2<<"wt "<<wt<<std::endl;
  cout<<"maturity2 "<<maturity2 <<"maturity"<<maturity<<std::endl;
  obs_sexr.allocate(1,nobs_fish);
  obs_sexr_srv_2.allocate(1,nsurv,1,nobs_srv_length);
  pred_sexr.allocate(styr,endyr);
}

void model_parameters::initializationfunction(void)
{
  F40.set_initial_value(.20);
  F35.set_initial_value(.21);
  F30.set_initial_value(.23);
  mean_log_rec.set_initial_value(mean_log_rec_prior);
  log_avg_fmort.set_initial_value(log_avg_fmort_prior);
  q_surv.set_initial_value(q_surv_prior_mean);
  fmort_dev.set_initial_value(0.00001);
  fishsel_params_f.set_initial_value(fishsel_prior_f);
  fishsel_params_m.set_initial_value(fishsel_prior_m);
  srv_params_f.set_initial_value(sel_prior_f);
  srv_params_m.set_initial_value(sel_prior_m);
  srv1desc_params_f.set_initial_value(sel1_desc_prior_f);
  srv1desc_params_m.set_initial_value(sel1_desc_prior_m);
  alpha.set_initial_value(1.);
  beta.set_initial_value(0.);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  q_surv.allocate(1,nsurv,q_Lower_bound,q_Upper_bound,q_Phase,"q_surv");
  fishsel_params_f.allocate(1,2,fishsel_LB_f,fishsel_UB_f,phase_fishery_sel,"fishsel_params_f");
  fishsel_params_m.allocate(1,2,fishsel_LB_m,fishsel_UB_m,phase_fishery_sel,"fishsel_params_m");
cout<<"sel_prior_f"<<sel_LB_f<<std::endl;
cout<<"sel_LB_m"<<sel_LB_m<<std::endl;
  srv_params_f.allocate(1,2*nsurv,sel_LB_f,sel_UB_f,phase_logistic_sel,"srv_params_f");
  srv_params_m.allocate(1,2*nsurv,sel_LB_m,sel_UB_m,phase_logistic_sel,"srv_params_m");
cout<<"srv_params_f"<<srv_params_f(1)<<std::endl;    
  srv1desc_params_f.allocate(1,2,sel1_desc_LB_f,sel1_desc_UB_f,phase_fishery_sel,"srv1desc_params_f");
  srv1desc_params_m.allocate(1,2,sel1_desc_LB_m,sel1_desc_UB_m,phase_fishery_sel,"srv1desc_params_m");
  alpha.allocate(phase_alphabeta,"alpha");
  beta.allocate(phase_alphabeta,"beta");
  mean_log_rec.allocate(1,"mean_log_rec");
  rec_dev.allocate(styr_rec,endyr-1,-15,15,2,"rec_dev");
  log_avg_fmort.allocate(2,"log_avg_fmort");
  fmort_dev.allocate(styr,endyr,-3,3,1,"fmort_dev");
  log_selcoffs_fish.allocate(1,2,1,nselages,phase_selcoffs,"log_selcoffs_fish");
  sexr_param_fish.allocate(1.0,1.0,-5,"sexr_param_fish");
  F40.allocate(0.01,1.,phase_F40,"F40");
  F35.allocate(0.01,1.,phase_F40,"F35");
  F30.allocate(0.01,1.,phase_F40,"F30");
  log_sel_fish.allocate(1,2,1,nages,"log_sel_fish");
  #ifndef NO_AD_INITIALIZE
    log_sel_fish.initialize();
  #endif
  sel.allocate(1,2,1,nages,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  sel_srv.allocate(1,2,1,nsurv,1,nages,"sel_srv");
  #ifndef NO_AD_INITIALIZE
    sel_srv.initialize();
  #endif
  avgsel_fish.allocate(1,2,"avgsel_fish");
  #ifndef NO_AD_INITIALIZE
    avgsel_fish.initialize();
  #endif
  popn.allocate(1,2,styr,endyr,"popn");
  #ifndef NO_AD_INITIALIZE
    popn.initialize();
  #endif
  totn_srv.allocate(1,nsurv,1,2,styr,endyr,"totn_srv");
  #ifndef NO_AD_INITIALIZE
    totn_srv.initialize();
  #endif
  temp1.allocate(1,nages,"temp1");
  #ifndef NO_AD_INITIALIZE
    temp1.initialize();
  #endif
  temp2.allocate(1,nages,"temp2");
  #ifndef NO_AD_INITIALIZE
    temp2.initialize();
  #endif
  explbiom.allocate(styr,endyr,"explbiom");
  #ifndef NO_AD_INITIALIZE
    explbiom.initialize();
  #endif
  pred_bio.allocate(styr,endyr,"pred_bio");
  #ifndef NO_AD_INITIALIZE
    pred_bio.initialize();
  #endif
  fspbio.allocate(styr,endyr,"fspbio");
  #ifndef NO_AD_INITIALIZE
    fspbio.initialize();
  #endif
  pred_srv.allocate(1,nsurv,styr,endyr,"pred_srv");
  #ifndef NO_AD_INITIALIZE
    pred_srv.initialize();
  #endif
  pred_p_fish.allocate(1,2,styr,endyr,1,nlen,"pred_p_fish");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish.initialize();
  #endif
  pred_p_srv_age_fem.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"pred_p_srv_age_fem");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_age_fem.initialize();
  #endif
  pred_p_srv_age_mal.allocate(1,nsurv_aged,1,nobs_srv_age,1,nages,"pred_p_srv_age_mal");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_age_mal.initialize();
  #endif
  pred_p_srv_len_fem.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"pred_p_srv_len_fem");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_len_fem.initialize();
  #endif
  pred_p_srv_len_mal.allocate(1,nsurv,1,nobs_srv_length,1,nlen,"pred_p_srv_len_mal");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv_len_mal.initialize();
  #endif
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  natage.allocate(1,2,styr,endyr,1,nages,"natage");
  #ifndef NO_AD_INITIALIZE
    natage.initialize();
  #endif
  totalbiomass.allocate(styr,endyr,"totalbiomass");
  catage.allocate(1,2,styr,endyr,1,nages,"catage");
  #ifndef NO_AD_INITIALIZE
    catage.initialize();
  #endif
  Z.allocate(1,2,styr,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(1,2,styr,endyr,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(1,2,styr,endyr,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  fmort.allocate(styr,endyr,"fmort");
  #ifndef NO_AD_INITIALIZE
    fmort.initialize();
  #endif
  rbar.allocate("rbar");
  #ifndef NO_AD_INITIALIZE
  rbar.initialize();
  #endif
  surv.allocate(1,2,"surv");
  #ifndef NO_AD_INITIALIZE
    surv.initialize();
  #endif
  offset.allocate(1,6,"offset");
  #ifndef NO_AD_INITIALIZE
    offset.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  catch_like.allocate("catch_like");
  #ifndef NO_AD_INITIALIZE
  catch_like.initialize();
  #endif
  sexr_like.allocate("sexr_like");
  #ifndef NO_AD_INITIALIZE
  sexr_like.initialize();
  #endif
  age_like.allocate(1,nsurv_aged,"age_like");
  #ifndef NO_AD_INITIALIZE
    age_like.initialize();
  #endif
  length_like2.allocate(1,nsurv+1,"length_like2");
  #ifndef NO_AD_INITIALIZE
    length_like2.initialize();
  #endif
  sel_like.allocate(1,4,"sel_like");
  #ifndef NO_AD_INITIALIZE
    sel_like.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  surv_like.allocate(1,nsurv,"surv_like");
  #ifndef NO_AD_INITIALIZE
    surv_like.initialize();
  #endif
  endbiom.allocate("endbiom");
  depletion.allocate("depletion");
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  tmp.allocate("tmp");
  #ifndef NO_AD_INITIALIZE
  tmp.initialize();
  #endif
  pred_sexr.allocate(styr,endyr,"pred_sexr");
  #ifndef NO_AD_INITIALIZE
    pred_sexr.initialize();
  #endif
  sigmar.allocate("sigmar");
  #ifndef NO_AD_INITIALIZE
  sigmar.initialize();
  #endif
  ftmp.allocate("ftmp");
  #ifndef NO_AD_INITIALIZE
  ftmp.initialize();
  #endif
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  SBF40.allocate("SBF40");
  #ifndef NO_AD_INITIALIZE
  SBF40.initialize();
  #endif
  SBF35.allocate("SBF35");
  #ifndef NO_AD_INITIALIZE
  SBF35.initialize();
  #endif
  SBF30.allocate("SBF30");
  #ifndef NO_AD_INITIALIZE
  SBF30.initialize();
  #endif
  sprpen.allocate("sprpen");
  #ifndef NO_AD_INITIALIZE
  sprpen.initialize();
  #endif
  Nspr.allocate(1,4,1,nages,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  nage_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"nage_future");
  #ifndef NO_AD_INITIALIZE
    nage_future.initialize();
  #endif
  fspbiom_fut.allocate(1,4,styr_fut,endyr_fut,"fspbiom_fut");
  F_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"F_future");
  #ifndef NO_AD_INITIALIZE
    F_future.initialize();
  #endif
  Z_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"Z_future");
  #ifndef NO_AD_INITIALIZE
    Z_future.initialize();
  #endif
  S_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"S_future");
  #ifndef NO_AD_INITIALIZE
    S_future.initialize();
  #endif
  catage_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"catage_future");
  #ifndef NO_AD_INITIALIZE
    catage_future.initialize();
  #endif
  avg_rec_dev_future.allocate("avg_rec_dev_future");
  #ifndef NO_AD_INITIALIZE
  avg_rec_dev_future.initialize();
  #endif
  avg_F_future.allocate(1,4,"avg_F_future");
  #ifndef NO_AD_INITIALIZE
    avg_F_future.initialize();
  #endif
  catch_future.allocate(1,3,styr_fut,endyr_fut,"catch_future");
  future_biomass.allocate(1,4,styr_fut,endyr_fut,"future_biomass");
  explbiom_fut.allocate(styr_fut,endyr_fut,"explbiom_fut");
  #ifndef NO_AD_INITIALIZE
    explbiom_fut.initialize();
  #endif
  maxsel_fish.allocate("maxsel_fish");
  #ifndef NO_AD_INITIALIZE
  maxsel_fish.initialize();
  #endif
  maxsel_srv.allocate(1,nsurv,"maxsel_srv");
  #ifndef NO_AD_INITIALIZE
    maxsel_srv.initialize();
  #endif
  mlike.allocate("mlike");
  #ifndef NO_AD_INITIALIZE
  mlike.initialize();
  #endif
  qlike.allocate("qlike");
  #ifndef NO_AD_INITIALIZE
  qlike.initialize();
  #endif
  flike.allocate("flike");
  #ifndef NO_AD_INITIALIZE
  flike.initialize();
  #endif
  qtime.allocate(styr,endyr,"qtime");
  #ifndef NO_AD_INITIALIZE
    qtime.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  obs_mean_sexr=0.34;  //initial value for avg proportion of male population estimated from shelf surveys; calculated below
  obs_SD_sexr=0.0485;  //initial value for standard deviation of mean male population proportion: calculated below
  
  for(i=1; i<=nobs_fish;i++)
  {  
	 male_len=0; fem_len=0; 
  for(j=first_length;j<=last_length;j++)
    {
    male_len+=obs_p_fish(1,i,j); 
    fem_len+=obs_p_fish(2,i,j); 
    }  
  obs_sexr(i)=male_len/(male_len+fem_len);
  }      
 
  for(i=1;i<=nsurv;i++)
  {
	for (j=1;j<=nobs_srv_length(i);j++)
	{    
	  male_len=0; fem_len=0;
		for (k=first_length;k<=last_length;k++)
		{     
		male_len+=obs_p_srv_length_mal(i,j,k);
		fem_len+=obs_p_srv_length_fem(i,j,k);	
		} 
		obs_sexr_srv_2(i,j)=male_len/(male_len+fem_len);  
    //	obs_sexr_srv_2(i,j)=sum(obs_p_srv_length_mal(i,j)/
    //	      (sum(obs_p_srv_length_mal(i,j))+sum(obs_p_srv_length_fem(i,j))));
	}
  }
  obs_mean_sexr=mean(obs_sexr_srv_2(1)); //previously was just estimated from shelf survey data so kept that here.
  obs_SD_sexr=std_dev(obs_sexr_srv_2(1));
 //Compute offset for multinomial and length bin proportions
 // offset is a constant nplog(p) is added to the likelihood     
 // magnitude depends on nsamples(sample size) and p's_
  offset.initialize();
  
  for (i=1; i <= nobs_fish; i++)
  {
    double sumtot ; 
    fem_len=0;
    male_len=0; 
    for(j=first_length;j<=last_length;j++)
    {
    fem_len+=obs_p_fish(1,i,j);
    male_len+=obs_p_fish(2,i,j);
    }
    sumtot=fem_len+male_len;    
  //  sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i));  
     obs_p_fish(1,i) = obs_p_fish(1,i) / sumtot; 
    obs_p_fish(2,i) = obs_p_fish(2,i) / sumtot; 
    for(j=first_length;j<=last_length;j++)
    {           
	  for (k=1;k<=2;k++){  
      obs_p_fish2(k,i,j)=obs_p_fish(k,i,j);
  //    obs_p_fish2(k,i,j-first_length+1)=value(obs_p_fish(k,i,j)); 
   // 			cout<<"obs_p_fish2(1,1) "<<obs_p_fish2(1,1)<<std::endl; 
   // 	cout<<"j "<<j<<std::endl; 
   // 			cout<<"j-first_length_+1 "<<j-first_length+1<<std::endl;
      }
    }
    for(k=1; k<=2;k++) 
    {
    //  offset(1) -= nsamples_fish(k,i)*obs_p_fish(k,i) * log(obs_p_fish(k,i)+.0001); //this multiplies elements together then sums them up.
      offset(1) -= nsamples_fish(k,i)*obs_p_fish2(k,i) * log(obs_p_fish2(k,i)+.0001);
    } 
  }  
  //this loops over all surveys and makes sure all proportions sum to 1.
  for(i=1;i<=nsurv;i++){
	for(j=1;j<=nobs_srv_length(i);j++){    
		double sumtot;   
		
		fem_len=0;
	    male_len=0; 
	    for(k=first_length;k<=last_length;k++)
	    {
	    fem_len+=obs_p_srv_length_fem(i,j,k);
	    male_len+=obs_p_srv_length_mal(i,j,k);
	    }
	    sumtot=fem_len+male_len;
		sumtot=sum(obs_p_srv_length_fem(i,j)+obs_p_srv_length_mal(i,j));
    obs_p_srv_length_mal(i,j)=obs_p_srv_length_mal(i,j)/sumtot;  //changing these to proportions rather than numbers
    obs_p_srv_length_fem(i,j)=obs_p_srv_length_fem(i,j)/sumtot; 
       for (k=first_length;k<=last_length;k++)
       {
       obs_p_srv_length_fem2(i,j,k)= obs_p_srv_length_fem(i,j,k);
       obs_p_srv_length_mal2(i,j,k)= obs_p_srv_length_mal(i,j,k);
       }
  //      offset(i+1)-= nsamples_srv_length_fem(i,j)*obs_p_srv_length_fem(i,j)*log(obs_p_srv_length_fem(i,j)+.0001)
  //                 +nsamples_srv_length_mal(i,j)*obs_p_srv_length_mal(i,j)*log(obs_p_srv_length_mal(i,j)+.0001); 
	    offset(i+1)-= nsamples_srv_length_fem(i,j)*obs_p_srv_length_fem2(i,j)*log(obs_p_srv_length_fem2(i,j)+.0001)
		           +nsamples_srv_length_mal(i,j)*obs_p_srv_length_mal2(i,j)*log(obs_p_srv_length_mal2(i,j)+.0001); 
	}
  } 
 
  //survey age offsets
  //survey age offsets
  for (i=1;i<=nsurv_aged;i++)
  {
    for (j=1;j<=nobs_srv_age(i);j++)
    {
	double sumtot; 
    fem_len=0;
    male_len=0; 
    for(k=first_age;k<=last_age;k++)
    {
    fem_len+=obs_p_srv_age_fem(i,j,k);
    male_len+=obs_p_srv_age_mal(i,j,k);
    }
    sumtot=fem_len+male_len;
   //cout<<"sumtot1 "<<sumtot<<std::endl;
  //  sumtot=sum(obs_p_srv_age_fem(i,j)+obs_p_srv_age_mal(i,j));   
   //    cout<<"sumtot2 "<<sumtot<<std::endl;  
    obs_p_srv_age_fem(i,j)=obs_p_srv_age_fem(i,j)/sumtot;
    obs_p_srv_age_mal(i,j)=obs_p_srv_age_mal(i,j)/sumtot; 
    for (k=first_age;k<=last_age;k++)
    {
       obs_p_srv_age_fem2(i,j,k)= obs_p_srv_age_fem(i,j,k);
       obs_p_srv_age_mal2(i,j,k)= obs_p_srv_age_mal(i,j,k);
    }   
  //  offset(i+nsurv+1)-=nsamples_srv_age(i,1,j)*obs_p_srv_age_fem(i,j)*log(obs_p_srv_age_fem(i,j)+.0001)+
  //               nsamples_srv_age(i,2,j)*obs_p_srv_age_mal(i,j)*log(obs_p_srv_age_mal(i,j)+.0001);
 //   cout<<"offset1 "<<offset<<std::endl;
    offset(i+nsurv+1)-=nsamples_srv_age(i,1,j)*obs_p_srv_age_fem2(i,j)*log(obs_p_srv_age_fem2(i,j)+.0001)+
    			             nsamples_srv_age(i,2,j)*obs_p_srv_age_mal2(i,j)*log(obs_p_srv_age_mal2(i,j)+.0001);
    }  
  }
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  ofstream& evalout= *pad_evalout;
   get_selectivity();     
   get_mortality();  
    surv(1)=mfexp(-1.0* M(1));
    surv(2)=mfexp(-1.0* M(2));   
   get_numbers_at_age();
   get_catch_at_age();
  if (active(F40))
    compute_spr_rates();
  if (last_phase())
  {
    Future_projections();      
  }
  if (sd_phase() || mceval_phase()) 
   { 
	Do_depend();
    if (mceval_phase())
    {
      evalout << obj_fun << " " ;
      for (i=styr;i<=endyr;i++)
      evalout<<  fspbio(i) << " ";
      for (i=styr;i<=endyr;i++)    
      evalout<< natage(1,i)*wt(1) + natage(2,i)*wt(2) <<" ";
      for (i=styr;i<=endyr;i++)
      evalout << 2*natage(1,i,1) <<" ";
      evalout <<  endl;
    }
  }
    evaluate_the_objective_function();  
}

void model_parameters::get_selectivity(void)
{
  ofstream& evalout= *pad_evalout;
  if(active(log_selcoffs_fish))// 
 {
    for(k=1;k<=2;k++)
    {
     for (j=1;j<=nselages;j++)
      {
        log_sel_fish(k,j)=log_selcoffs_fish(k,j);
      }
    }
   //sets selectivity of ages older than nselages to selectivity at nselages  
   for(k=1;k<=2;k++)
   {
     for (j=nselages+1;j<=nages;j++)
     {
       log_sel_fish(k,j)=log_sel_fish(k,j-1);
     }
   }
 for(k=1;k<=2;k++)
   {
     avgsel_fish(k)=log(mean(mfexp(log_selcoffs_fish(k))));
   }
 //vector=vector-scalar same as  vector-=scalar  
 //scaling selectivities by subracting the mean so exp(mean(s))=1.   
 //selectivities can be greater than 1 but mean is 1. 
 
  for(k=1;k<=2;k++)
    {
      log_sel_fish(k)-=log(mean(mfexp(log_sel_fish(k))));
      sel(k)=mfexp(log_sel_fish(k));
      if(k==2)   //males
       {
         sel(k)=sel(k)*sexr_param_fish; //fixed at 1 in GOA model not BSAI model
       }    
    }     
 }//  end if(active(log_selcoffs_fish))
  else
    {
     //logistic selectivity curve
          for (j=1;j<=nages;j++)
          { 
            if(j<=nselages)
             {
               sel(1,j)=1./(1.+mfexp(-1.*fishsel_params_f(1)*(double(j)-fishsel_params_f(2))));
               sel(2,j)=1./(1.+mfexp(-1.*fishsel_params_m(1)*(double(j)-fishsel_params_m(2))));
             }
            else
            {
             sel(1,j)=sel(1,j-1);
             sel(2,j)=sel(2,j-1);
            }    				
          } 
     }
  for(i=1;i<=nsurv;i++)
  {
     if(nsel_params(i)==4)
     {
	 sel_srv(1,i) = get_sel(srv_params_f(i),srv_params_f(i+1),srv1desc_params_f(1),srv1desc_params_f(2));
	 sel_srv(2,i) = get_sel(srv_params_m(i),srv_params_m(i+1),srv1desc_params_m(1),srv1desc_params_m(2));         
     }
       else
       {
       sel_srv(1,i) = get_sel(srv_params_f((2*i)-1),srv_params_f(2*i));
       sel_srv(2,i) = get_sel(srv_params_m((2*i)-1),srv_params_m(2*i)); 
       }        
  }
  
}

dvar_vector model_parameters::get_sel(const dvariable& slp, const dvariable& a50)
{
  ofstream& evalout= *pad_evalout;
   {
	dvar_vector sel_tmp(1,nages);
    for (j=1;j<=nages;j++)  //this is selectivity for the surveys
    {  
    sel_tmp(j)=1./(1.+mfexp(-1.*slp*(double(j)-a50))); 
    }          
    return(sel_tmp);
   }          
         
}

dvar_vector model_parameters::get_sel(const dvariable& slp, const dvariable& a50, const dvariable& dslp, const dvariable& d50)
{
  ofstream& evalout= *pad_evalout;
   {
	dvar_vector sel_tmp(1,nages);
   for (j=1;j<=nages;j++)  //this is selectivity for the surveys         
   {
	  sel_tmp(j) = 1./(1.+mfexp(-1.*slp*(double(j)-a50)));           
      sel_tmp(j) *= 1./(1.+mfexp(dslp*(double(j)-d50)));
   }
 	return(sel_tmp);
  }          
}

void model_parameters::get_mortality(void)
{
  ofstream& evalout= *pad_evalout;
  maxsel_fish=max(sel(1));     //1 is females
  if(maxsel_fish<max(sel(2)))  //if highest female selectivity is > male selectivity, make maxsel_fish=male high selectivity
      maxsel_fish=max(sel(2));
  fmort = mfexp(log_avg_fmort+fmort_dev); 
  for(k=1;k<=2;k++)
  {
    for (i=styr;i<=endyr;i++)
    {
      F(k,i)=(sel(k))*fmort(i);
      Z(k,i)=F(k,i) + M(k); 
    }
  } 
  S = mfexp(-1.*Z); 
}

void model_parameters::get_numbers_at_age(void)
{
  ofstream& evalout= *pad_evalout;
  maxsel_fish=max(sel(1));   
  if(maxsel_fish<max(sel(2)))//if females greater than males, then set the max to the females.
    maxsel_fish=max(sel(2)); //set max to whichever sex is larger
  for(i=1;i<=nsurv;i++)
  {
   maxsel_srv(i)=max(sel_srv(1,i));
   if(maxsel_srv(i)<max(sel_srv(2,i)))
   maxsel_srv(i)=max(sel_srv(2,i));
  }     
  int itmp;
 //calc initial population  
  for (j=1;j<nages;j++)
    {
      itmp=styr+1-j;
      natage(1,styr,j)=mfexp(mean_log_rec-(M(1)*double(j-1))+rec_dev(itmp));
      natage(2,styr,j)=mfexp(mean_log_rec-(M(2)*double(j-1))+rec_dev(itmp));
    }
    itmp=styr+1-nages;
  //last age    
    natage(1,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(1)*(nages-1)))/(1.- surv(1));
    natage(2,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(2)*(nages-1)))/(1.- surv(2));
 // Now do for next several years----------------------------------
  for (i=styr+1;i<=endyr;i++)
  {
    //for age 1 recruits in the last year use value read in from data file
    if(i<=(endyr-1))
    {
      natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
      natage(2,i,1)=natage(1,i,1);
    }
    else
    {
      natage(1,i,1)=median_rec;
      natage(2,i,1)=natage(1,i,1);
    }
  }
 //numbers at age 
  for(k=1;k<=2;k++)
  {
    for (i=styr;i< endyr;i++)
    {
      for(j=1;j<nages;j++)
      {
        natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
      }
      natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
      popn(k,i)= natage(k,i)*sel(k);
    }
    popn(k,endyr)=natage(k,endyr)*sel(k);
  }
  for (i=styr;i<=endyr;i++)
  {
      pred_sexr(i)=sum(natage(2,i))/(sum((natage(1,i)+natage(2,i))));  //calculation of prop. of males in pred. population 
  }
  //predicted survey values
  fspbio.initialize(); 
  qtime=q_surv(1); 
  
  for (j=1;j<=nsurv;j++)
 {
    for (i=styr;i<=endyr;i++)
   {
  fspbio(i) = natage(1,i)*elem_prod(wt(1),maturity);
  explbiom(i)=0.;
  pred_bio(i)=0.; 
  pred_srv(j,i)=0.;
  //catchability calculation for survey years
  if (assess==1 && (i>=1982) && (i-1981 <= nobs_srv(1)))      //JNI catchability calculation for survey years    
   {   
   qtime(i)=q_surv(1)*mfexp(-alpha+beta*bottom_temps(i-1981));
   }
  for(k=1;k<=2;k++)
    {
    if (j==1 && assess==1)
      {             
    pred_srv(j,i) += qtime(i)*(natage(k,i)*elem_prod(sel_srv(k,j)/maxsel_srv(j),wt(k)));maxsel_srv(j);   //shelf survey, dividing by the maxsel constrains female selectivity to be 1.0
      } 
    else 
      {
    pred_srv(j,i) += q_surv(j)*(natage(k,i)*elem_prod(sel_srv(k,j),wt(k)));///maxsel_srv(j);         //slope survey JNI  do not need to divide by maxsel_srv if it is logistic but does not hurt
      } 
     
       //Aleutian Islands survey JNI
    //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
    // are set different in the tpl file the program will take to value from the bin file and use that 
    explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
    pred_bio(i)+=natage(k,i)*wt(k);
      }
   }  
 }     
  //Fit survey length compositions
  for (i=1;i<=nsurv;i++)
  {
	for (j=1;j<=nobs_srv_length(i);j++)
	{
			ii=yrs_srv_length(i,j);
			pred_p_srv_len_fem(i,j)=q_surv(i)*elem_prod(sel_srv(1,i),natage(1,ii))*lenage(1);
			pred_p_srv_len_mal(i,j)=q_surv(i)*elem_prod(sel_srv(2,i),natage(2,ii))*lenage(2);
			dvariable sum_tot=sum(pred_p_srv_len_fem(i,j)+pred_p_srv_len_mal(i,j));
			pred_p_srv_len_fem(i,j)/=sum_tot;
			pred_p_srv_len_mal(i,j)/=sum_tot;
	 }
   }
     
  //Fit survey age composition
   for (i=1;i<=nsurv_aged;i++)
  {
	for (j=1;j<=nobs_srv_age(i);j++)
	{
			ii=yrs_srv_age(i,j);
			pred_p_srv_age_fem(i,j)=q_surv(i)*elem_prod(sel_srv(1,i),natage(1,ii));
			pred_p_srv_age_mal(i,j)=q_surv(i)*elem_prod(sel_srv(2,i),natage(2,ii));
			dvariable sum_tot=sum(pred_p_srv_age_fem(i,j)+pred_p_srv_age_mal(i,j));
			pred_p_srv_age_fem(i,j)/=sum_tot;
			pred_p_srv_age_mal(i,j)/=sum_tot;
	 }
   }
  depletion=pred_bio(endyr)/pred_bio(styr);
  endbiom=pred_bio(endyr);
 
}

void model_parameters::get_catch_at_age(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr; i<=endyr; i++)
  {
    pred_catch(i)=0.;
    for(k=1;k<=2;k++)
    {      
      //--Baranov's equation here-----------------------------------
      for (j = 1 ; j<= nages; j++)
      {
        catage(k,i,j) = natage(k,i,j)*F(k,i,j)*(1.-S(k,i,j))/Z(k,i,j);
        pred_catch(i) += catage(k,i,j)*wt(k,j);
      }
      pred_p_fish(k,i)=elem_prod(sel(k),natage(k,i))*lenage(k)/(popn(1,i)+popn(2,i));
    }
  }   
    
}

void model_parameters::Future_projections(void)
{
  ofstream& evalout= *pad_evalout;
  for(k=1;k<=2;k++)
  {
    nage_future(k,styr_fut)(2,nages)=++elem_prod(natage(k,endyr)(1,nages-1),S(k,endyr)(1,nages-1));
    nage_future(k,styr_fut,nages)+=natage(k,endyr,nages)*S(k,endyr,nages);
   }
    future_biomass.initialize();
    catch_future.initialize();
    for (int l=1;l<=4;l++)
    {
      switch (l)
      {
        case 1:
          ftmp=F40;
          break;
        case 2:
          ftmp=F35;
          break;
        case 3:
          ftmp=F30;
          break;
        case 4:
          ftmp.initialize();
          break;
      }
     for(k=1;k<=2;k++)
     {
      for (i=endyr+1;i<=endyr_fut;i++)
      {
        for (j=1;j<=nages;j++)
        {
          F_future(k,i,j) = (sel(k,j)/maxsel_fish)*ftmp;
          Z_future(k,i,j) = F_future(k,i,j)+M(k);
          S_future(k,i,j) = exp(-1.*Z_future(k,i,j));
        }
      }
      for (i=styr_fut;i<endyr_fut;i++)
      {
        nage_future(k,i,1)  = median_rec;
       // Now graduate for the next year....
        nage_future(k,i+1)(2,nages) = ++elem_prod(nage_future(k,i)(1,nages-1),S_future(k,i)(1,nages-1));
        nage_future(k,i+1,nages)   += nage_future(k,i,nages)*S_future(k,i,nages);
      }
      nage_future(k,endyr_fut,1)  = median_rec;
      // Now get catch at future ages
      for (i=styr_fut; i<=endyr_fut; i++)
      {
        for (j = 1 ; j<= nages; j++)
        {
          catage_future(k,i,j) = nage_future(k,i,j) * F_future(k,i,j) * ( 1.- S_future(k,i,j) ) / Z_future(k,i,j);
         if(k==1)
          {
          fspbiom_fut(l,i) += nage_future(1,i,j)*wt(1,j)*maturity(j);
          }
        }
        if (l!=4) catch_future(l,i)   += catage_future(k,i)*wt(k);
        future_biomass(l,i) += nage_future(k,i)*wt(k);
 
      }   //end loop over future years
     }   //end loop over sex
     fspbiom_fut(l)=0.;
     for(i=styr_fut;i<=endyr_fut;i++)
       fspbiom_fut(l,i) = elem_prod(nage_future(1,i),wt(1)) * maturity;
    }   //End of loop over F's
}

void model_parameters::compute_spr_rates(void)
{
  ofstream& evalout= *pad_evalout;
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0.initialize();
  SBF40.initialize();
  SBF35.initialize();
  SBF30.initialize();
  // Initialize the recruit (1) for each F  (F40 etc)
  for (i=1;i<=3;i++)
    Nspr(i,1)=1.;
  for (j=2;j<nages;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*exp(-1.*M(1));
    Nspr(2,j)=Nspr(2,j-1)*exp(-1.*(M(1)+F40*sel(1,j-1)/maxsel_fish));
    Nspr(3,j)=Nspr(3,j-1)*exp(-1.*(M(1)+F35*sel(1,j-1)/maxsel_fish));
    Nspr(4,j)=Nspr(4,j-1)*exp(-1.*(M(1)+F30*sel(1,j-1)/maxsel_fish));
  }
 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*exp(-1.*M(1))/(1.-exp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*exp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*exp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
  Nspr(4,nages)=Nspr(4,nages-1)*exp(-1.* (M(1)+F30*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F30*sel(1,nages)/maxsel_fish)));
  for (j=1;j<=nages;j++)
  {
   // Kill them off till april (0.25) atf spawn in winter so put in 0.0
   //         Number   ProportMat  Wt    Amount die off prior to spawning (within that year)
    SB0    += Nspr(1,j)*maturity(j)*wt(1,j)*exp(-0.0*M(1));
    SBF40  += Nspr(2,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F40*sel(1,j)/maxsel_fish));
    SBF35  += Nspr(3,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F35*sel(1,j)/maxsel_fish));
    SBF30  += Nspr(4,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F30*sel(1,j)/maxsel_fish));
  }
  sprpen    = 200.*square((SBF40/SB0)-0.4);
  sprpen   += 200.*square((SBF35/SB0)-0.35);
  sprpen   += 200.*square((SBF30/SB0)-0.30);
}

void model_parameters::Do_depend(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr;  i<=endyr;  i++) 
  totalbiomass(i)=natage(1,i)*wt(1) + natage(2,i)*wt(2);
  obj_fun += 1.*sexr_like;             // male proportion prior, emphasis factor = 1
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& evalout= *pad_evalout;
  length_like2.initialize();
  age_like.initialize();
  fpen.initialize();
  rec_like.initialize();
  surv_like.initialize();
  catch_like.initialize();
  sexr_like.initialize();
  obj_fun.initialize();
  if (active(rec_dev))
  {
  length_like2.initialize();
  int ii;
    //recruitment likelihood - norm2 is sum of square values   
    rec_like = norm2(rec_dev);
    for(k=1;k<=2;k++)
    {
      for (i=1; i <= nobs_fish; i++)
      {
        ii=yrs_fish(i);
        //fishery length likelihood fitting
          length_like2(1) -= nsamples_fish(k,i)*(1e-5+obs_p_fish(k,i))*log(pred_p_fish(k,ii)+1e-5);
      }
    }
    //add the offset to the likelihood   
    length_like2(1)-=offset(1);
  //survey length composition fitting
   for (i=1;i<=nsurv;i++)
   {
     for (j=1;j<=nobs_srv_length(i);j++) 
     {    
	   length_like2(i+1)-=((nsamples_srv_length_fem(i,j)*(1e-3+obs_p_srv_length_fem(i,j))*log(pred_p_srv_len_fem(i,j)+1e-3))
	                      +(nsamples_srv_length_mal(i,j)*(1e-3+obs_p_srv_length_mal(i,j))*log(pred_p_srv_len_mal(i,j)+1e-3)));
	 } 
	  length_like2(i+1)-=offset(i+1); 
    } 
 
  for (i=1;i<=nsurv_aged;i++)
  {
	for (j=1;j<=nobs_srv_age(i);j++)
	{
   age_like(i)-=nsamples_srv_age(i,1,j)*(1e-3+obs_p_srv_age_fem(i,j))*log(pred_p_srv_age_fem(i,j)+1e-3)+
                  nsamples_srv_age(i,2,j)*(1e-3+obs_p_srv_age_mal(i,j))*log(pred_p_srv_age_mal(i,j)+1e-3); 	
	}	
	age_like(i)-=offset(i+nsurv+1);
  }
  //end of if(active (rec_dev))
  }
  // Fit to indices (lognormal)  
  //weight each years estimate by 1/(2*variance) - use cv as an approx to s.d. of log(biomass) 
  for (i=1;i<=nsurv;i++)
  {   
  surv_like(i) = norm2(elem_div(log(obs_srv(i))-log(pred_srv(i)(yrs_srv(i))),sqrt(2)*cv_srv(i)));
  }  
   
   catch_like=norm2(log(catch_bio+.000001)-log(pred_catch+.000001));
   // sex ratio likelihood
   sexr_like=0.5*norm2((obs_mean_sexr-pred_sexr)/obs_SD_sexr); 
 //selectivity likelihood is penalty on how smooth selectivities are   
 //here are taking the sum of squares of the second differences  
  if(active(log_selcoffs_fish))
  {  
    sel_like(1)=wt_like(1)*norm2(first_difference(first_difference(log_sel_fish(1)))); //fishery females
    sel_like(3)=wt_like(3)*norm2(first_difference(first_difference(log_sel_fish(2)))); //fishery males 
   for (j=1;j<nages;j++)
   {
    if(monot_sel==1)
    { 
        if (log_sel_fish(1,j)>log_sel_fish(1,j+1))
        sel_like(1)+=wt_like(5)*square(log_sel_fish(1,j)-log_sel_fish(1,j+1));   //monotonicity constraint fishery females
        if (log_sel_fish(2,j)>log_sel_fish(2,j+1))
        sel_like(3)+=wt_like(6)*square(log_sel_fish(2,j)-log_sel_fish(2,j+1));  //monotonicity constraing fishery males
    }
   } 
    obj_fun+=1.*sum(sel_like);    
    obj_fun+=1.*square(avgsel_fish(1));
    obj_fun+=1.*square(avgsel_fish(2));  
  } //end if active(log_selcoffs_fish)
  // Phases less than 3, penalize High F's
    if (current_phase()<2)
    {
       //F's are low for arrowtooth changed the value to compare from .2 to .001
       //don't know if makes any difference since the penalty is reduced at the end
       fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-.01);
    }
    else
    {
      fpen=.001*norm2(mfexp(fmort_dev+log_avg_fmort)-.01);
    }
    if (active(fmort_dev))
    {
      fpen+=.01*norm2(fmort_dev);
    }
 
  obj_fun += rec_like;
  obj_fun += 1.*sum(length_like2);  
  obj_fun+= 1.*sum(age_like); 
  obj_fun += 1.*sum(surv_like);     
  obj_fun += 300*catch_like;        // large emphasis to fit observed catch 
  obj_fun += fpen;   
  obj_fun += sprpen;     
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Estimated numbers of female fish at age " << endl;
    for (i=styr; i<=endyr;i++)
      report << i <<natage(1,i)<< endl;
   
  report << "Estimated numbers of male fish at age " << endl;
    for (i=styr; i<=endyr;i++)
      report << i  <<natage(2,i)<<endl;
  report << "Estimated catch numbers for females at age " << endl;
    for (i=styr; i<=endyr;i++)
      report << i << catage(1,i)<< endl;
   
  report << "Estimated catch numbers for males at age  " << endl;
    for (i=styr; i<=endyr;i++)
      report << i <<catage(2,i)<<endl;
  report << "Estimated female F mortality at age " << endl;
    for (i=styr; i<=endyr;i++)
      report <<  i <<F(1,i)<< endl;
   
  report << "Estimated male F mortality at age " << endl;
    for (i=styr; i<=endyr;i++)
      report <<  i <<F(2,i)<<endl;
  report << "Estimated fishery selectivity for females at age " << endl;
    for (j=1; j<=nages;j++)
      report << j<<" " <<sel(1,j)/maxsel_fish<< endl;
    
  report << "Estimated fishery selectivity for males at age " << endl;
    for (j=1; j<=nages;j++)
      report <<  j<<" "  <<sel(2,j)/maxsel_fish<<endl;
  report << "Estimated survey selectivity for females at age " << endl;
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
  for (j=1; j<=nages;j++)
      report << j<<" " <<sel_srv(1,k,j)/maxsel_srv(1)<< endl;
  }  
 
  report << "Estimated survey selectivity for males at age " << endl;    
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
    for (j=1; j<=nages;j++)
      report <<  j<<" "  <<sel_srv(2,k,j)/maxsel_srv(1)<<endl;
   }
  report << endl << "Survey biomass (Year, Obs_biomass, Pred_biomass) " << endl; 
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv(k);i++)
      report << yrs_srv(k,i) << ","<< obs_srv(k,i) << "," << pred_srv(k,yrs_srv(1,i)) << endl; 
  }
 
  report <<" Observed female survey length composition " << endl; 
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_length(k); i++)
      report << yrs_srv_length(k,i) << obs_p_srv_length_fem(k,i) << endl; 
   }
   
  report <<" Predicted female survey length composition " << endl; 
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_length(k); i++)
      report << yrs_srv_length(k,i) << pred_p_srv_len_fem(k,i) << endl; 
  }
  report <<" Observed male survey length composition " << endl; 
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_length(k); i++)
      report << yrs_srv_length(k,i) << obs_p_srv_length_mal(k,i) << endl;  
  }
  report <<" Predicted male survey length composition " << endl; 
  for (k=1;k<=nsurv;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_length(k); i++)
      report << yrs_srv_length(k,i)  << pred_p_srv_len_mal(k,i) << endl;  
  }
  report <<" Observed female survey age composition " << endl; 
  for (k=1;k<=nsurv_aged;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_age(k); i++)
      report << yrs_srv_age(k,i) << obs_p_srv_age_fem(k,i) << endl;
   }
  report <<" Predicted female survey age composition " << endl;  
  for (k=1;k<=nsurv_aged;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_age(k); i++)
      report << yrs_srv_age(k,i)  << pred_p_srv_age_fem(k,i) << endl; 
  }     
  report <<" Observed male survey age composition " << endl; 
  for (k=1;k<=nsurv_aged;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_age(k); i++)
      report << yrs_srv_age(k,i) << obs_p_srv_age_mal(k,i) << endl;
   }
  report <<" Predicted male survey age composition " << endl;  
  for (k=1;k<=nsurv_aged;k++)
  {  
  report<<"survey"<<k<<endl;
    for (i=1; i<=nobs_srv_age(k); i++)
      report << yrs_srv_age(k,i)  << pred_p_srv_age_mal(k,i) << endl; 
  }
  report <<" Observed Catch by year " << endl;
    for (i=styr;  i<=endyr; i++)
       report << i<<" " << catch_bio(i) << endl;
  report  << "Predicted Catch by year "<< endl;
    for (i=styr;  i<=endyr;  i++)
        report  <<  i <<" "<< pred_catch(i) << endl; 
  report  << "Female spawning biomass , Total biomass" <<  endl;
      for (i=styr;  i<=endyr;  i++)
         report  << i<<" " << fspbio(i) <<" "<<natage(1,i)*wt(1) + natage(2,i)*wt(2)<< endl;
 
  report << "$estimated.annual.fishing.mortality" << endl;
	  report << mfexp(log_avg_fmort+fmort_dev) << endl;
	
  report << endl<<endl;
  report << "F40" << endl;
  report << F40 << endl;
  report << "F35=  " << endl; 
  report << F35 << endl;
  report << "F30=  " << endl;
  report << F30 << endl;
  report << "spawning biomass per recruit at F40 harvest rate " << endl;
  report << SBF40<< endl;
  report << "spawning biomass per recruit at F35 harvest rate " << endl;
  report << SBF35 << endl;
  report << "spawning biomass per recruit at F30 harvest rate " << endl;
  report << SBF30 << endl;
  report << "spawning biomass per recruit with no fishing " << endl;
  report  << SB0 << endl;
  report << "Likelihood components" << endl;
  report << "Survey likelihood component " << endl;  
  for (k=1;k<=nsurv;k++)
  {
	report<<"survey"<<k<<endl;
  report << surv_like(k) << endl; 
  }
  report << "fishery length composition likelihood " << endl;
  report << length_like2(1) << endl;
  report << "recruitment likelihood component est.  " << endl;
  report << rec_like << endl;
  report << "catch likelihood component est.  " << endl;
  report << catch_like << endl;
  report << "sex ratio likelihood component " << endl;
  report << sexr_like << endl;
  report << "shelf survey age composition  " << endl; 
  report << age_like<< endl;    
  for (k=1;k<=nsurv_aged;k++)
  {
	report<<"survey"<<k<<endl;
  report << age_like<< endl; 
  }
  report << "Projected biomass" << endl;
  report << future_biomass << endl;
  report <<"projected future female spawning biomass " << endl;
  report <<fspbiom_fut << endl;
  report << "Future yield " << endl;
  report << catch_future << endl;  
  report << "survey q =" << endl; 
  for(k=1;k<=nsurv;k++)
  {
	 report<<"survey"<<k<<endl;
  report << q_surv(k) << endl;
  }
  report << " female natural mortality for this run" << endl;
  report << M(1) << endl;
  report << " male natural mortality for this run" << endl;
  report << M(2) << endl;
  report << "predicted male proportion in population" << endl;
   for (i=styr;i<=endyr;i++)
      report << i << " " << pred_sexr(i) << endl;
  report << "mean observed prop. male in shelf surveys = "<< endl;
  report << obs_mean_sexr << endl;
  report << "stdev of mean observed prop. male in shelf surveys = " << endl;
  report << obs_SD_sexr << endl;
  report << "standard error of biomass in surveys = " << endl; 
  for(k=1;k<=nsurv;k++)
  {
	 report<<"survey"<<k<<endl;
  report << obs_srv_sd(k) << endl;     
  }
  report << " recruits" << endl;
    for (i=styr;  i<=endyr;  i++)
         report  << i <<" " << 2*natage(1,i,1) <<" "<< endl;   
		
		report<<"mean_log_rec"<<mean_log_rec<<std::endl;
	
  report << " Go drink coffee " << endl; 
  
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{4000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3 1e-4 1e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(300);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
