//                      atfbsai2.tpl
//           Uses the shelf survey, slope survey and the Aleutian Islands (fits biomass and length comps)
//           modified version of the GOA model with domed shaped selectivity for the shelf survey males,
//           this run tries to estimate the temperature
//           effect on shelf survey catchability  
//			1 is female, 2 is male in model (not in survey data)   
			//this is how to change the name of the data file
// when you run it add 			./atfbsai_is2014_6 -ind atfbsai_is2014_4.dat
DATA_SECTION
!!CLASS ofstream evalout("atfbsai_is2014_np.mcmc.out");
  init_int styr         //(1) start year of model
  init_int endyr        //(2) end year
  init_int styr_fut     //(3) start year of projections (endyr+1) 
  init_int endyr_fut    //(4) end year of projections
  init_int phase_F40      //(5) phase F40 is estimated
  init_number median_rec  //(6) median recruit value to use for the last 3 years
  init_int nages          //(7) # of ages in the model      
 //selectivity is set to the selectivity at nselages-1 after age nselages 
  init_int nselages       //(8) fishery (for asymptotic selectivity) set to 19
  init_int nsurv     //number of surveys
  init_vector nselages_srv(1,nsurv) //vector or values to calculate selectivity for each survey
  init_int nselages_srv1  //(9) slope survey (for asymptotic selectivity) 
  init_int nselages_srv2  //(10) shelf survey (for asymptotic selectivity)
  init_int nselages_srv3  //(11) Aleutian Islands survey (for asymptotic selectivity)
  init_int monot_sel     //(51) selectivity smoothing function for fishery  
  init_int phase_logistic_sel //(12)
  init_int phase_selcoffs      //(52) generally set to phase 4 phase for smooth selectivity curve fishery       

 //sample size for length comps for weighting likelihoods  
  init_vector wt_like(1,8)    //(53) 
  init_int nlen             //(13) # of length bins
  init_int nobs_fish          //(14) # of years of fishery data
  init_ivector yrs_fish(1,nobs_fish)   //(15) years with fishery data
  init_matrix nsamples_fish(1,2,1,nobs_fish)  //(16) sample size (weights) for each sex and yr of fishery data
  init_int nobs_srv1          //(17) # of years of shelf survey biomass data
  init_int nobs_srv2          //(18) # of years of slope survey data
  init_int nobs_srv3          //(19) # of years of Aleutian Islands data
  init_ivector yrs_srv1(1,nobs_srv1)   //(20) years with shelf survey biomass data
  init_ivector yrs_srv2(1,nobs_srv2)   //(21) years with slope survey biomass data
  init_ivector yrs_srv3(1,nobs_srv3)   //(22) years with Aleutian Islands survey bioamass data
  init_int nobs_srv1_length          //(23) # yrs with shelf survey length data
  init_int nobs_srv2_length          //(24) # yrs with slope survey length data
  init_int nobs_srv3_length          //(25) # yrs with Aleutian Islands survey length data
  init_ivector yrs_srv1_length(1,nobs_srv1_length)    //(26) yrs with shelf survey length data
  init_ivector yrs_srv2_length(1,nobs_srv2_length)    //(27) yrs with slope survey length data
  init_ivector yrs_srv3_length(1,nobs_srv3_length)    //(28) yrs with Aleutian Islands survey length data
  init_int nobs_srv1_age                   //(47) # of years with shelf survey ages
  init_int nobs_srv3_age        //(54) # of years with ai survey ages
  init_ivector yrs_srv1_age(1,nobs_srv1_age)  //(48) years of shelf survey with ages
  init_ivector yrs_srv3_age(1,nobs_srv3_age)  //(55) years of ai survey with ages
  init_matrix nsamples_srv1_age(1,2,1,nobs_srv1_age)   //(49) sample size of ages read in each year, by sex
  init_matrix nsamples_srv3_age(1,2,1,nobs_srv3_age)   //(56) sample size of ages read in each year, by sex
  init_3darray obs_p_srv1_age(1,2,1,nobs_srv1_age,1,nages)  //(50) shelf survey age comps by sex and year
  init_3darray obs_p_srv3_age(1,2,1,nobs_srv3_age,1,nages)  //(57) AI survey age comps by sex and year
  init_matrix nsamples_srv1_length(1,2,1,nobs_srv1_length)  // (29) sample size for each length comp by sex and year from shelf survey
  init_matrix nsamples_srv2_length(1,2,1,nobs_srv2_length)  // (30) sample size for each length comp by sex and year from slope survey
  init_matrix nsamples_srv3_length(1,2,1,nobs_srv3_length)  //(31) sample size for each length comp by sex and year from Aleutian I survey
  init_3darray obs_p_srv1_length(1,2,1,nobs_srv1_length,1,nlen) //(32) shelf survey length comps by bin, sex and yr 
  init_3darray obs_p_srv2_length(1,2,1,nobs_srv2_length,1,nlen) //(33) slope survey length comps by bin, sex and yr
  init_3darray obs_p_srv3_length(1,2,1,nobs_srv3_length,1,nlen) //(34) Aleutian Islands survey length comps by bin, sex and yr
  init_3darray obs_p_fish(1,2,1,nobs_fish,1,nlen)  //(35) fishery length comps
  init_vector catch_bio(styr,endyr)    //(36) catch by year
  init_vector obs_srv1(1,nobs_srv1)    //(37) shelf survey biomass by year
  init_vector obs_srv1_sd(1,nobs_srv1) //(38) shelf survey SE by year
  init_vector obs_srv2(1,nobs_srv2)    //(39) slope survey biomass by year
  init_vector obs_srv2_sd(1,nobs_srv2) //(40) slope survey SE by year
  init_vector obs_srv3(1,nobs_srv3)    //(41) Aleutian Islands survey biomass by year 
  init_vector obs_srv3_sd(1,nobs_srv3) //(42) Aleutian Islands survey SE by year 
  init_matrix wt(1,2,1,nages)          //(43) weight-at-age by sex
  init_vector maturity(1,nages)        //(44) female prop. mature-at-age
  init_3darray lenage(1,2,1,nages,1,nlen)  //(45) length-age transition matrix
  init_vector bottom_temps(1,nobs_srv1)    //(46) shelf survey bottom temperatures
  int styr_rec; 
  vector cv_srv1(1,nobs_srv1);      //shelf survey CV
  vector cv_srv2(1,nobs_srv2);      //slope survey CV
  vector cv_srv3(1,nobs_srv3);      //Aleutian Islands survey CV 
 //year
  int i
//age
  int j
//sex
  int k
//
  int ii
  int m

 LOCAL_CALCS
   styr_rec=styr-nages+1;
   if(nselages>nages) nselages=nages;
   if(nselages_srv1>nages) nselages_srv1=nages;  
  if(nselages_srv2>nages) nselages_srv2=nages; 
  for (i=0;i<3;i++)
  {if (nselages_srv(i)>=nages)
  nselages_srv(i)=nages;
  }
   //calculate cv for surveys
    cv_srv1=elem_div(obs_srv1_sd,obs_srv1);   //shelf survey CV
    cv_srv2=elem_div(obs_srv2_sd,obs_srv2);   //slope survey CV
    cv_srv3=elem_div(obs_srv3_sd,obs_srv3);   //Aleutian Island survey CV
   //change weights to tons
   wt=wt*.001;
 END_CALCS

  vector obs_sexr(1,nobs_fish)  // prop. females in fishery length data
  vector obs_sexr_srv1_2(1,nobs_srv1_length) // prop. males in shelf survey length data
  vector obs_sexr_srv2_2(1,nobs_srv2_length) // prop. males in slope survey length data
  vector obs_sexr_srv3_2(1,nobs_srv3_length) // prop. males in Aleutian Islands length data
  number obs_mean_sexr    //average proportion of males in shelf survey population estimates
  number obs_SD_sexr      //standard deviation from male prop. in shelf survey population estimates
  vector pred_sexr(styr,endyr)   //proportion of males in num at age matrix to be calculated
  4darray test(1,3,1,2,1,4,1,5)
  !! cout <<"test"<<endl; 
  !! cout <<test<<endl; 
  matrix testmat(4,5);    
  !! cout <<testmat<<endl;

INITIALIZATION_SECTION
  //can have different mortality for males and females
  F40 .20
  F35 .21
  F30 .23
  mean_log_rec 10.
  log_avg_fmort -5.
  q1 .75   //proportion of the estimated BSAI survey biomass that is on the shelf
  q2 .10   //proportion of the estimated BSAI survey biomass that is on the slope
  q3 .14   //proportion of the estimated BSAI survey biomass that is in the Aleutian Islands
  fmort_dev 0.00001
  fish_slope_f .4
  fish_sel50_f  5.
  fish_slope_m  .1
  fish_sel50_m  8
  srv1_slope_f1  .8
  srv1_slope_f2  .8
  srv1_slope_m1  .4
  srv1_sel50_f1  4.
  srv1_sel50_f2  4.
  srv1_slope_m2 .4
  srv1_sel50_m1  8.
  srv1_sel50_m2  8.
  srv2_slope_f  .4
  srv2_sel50_f  8.
  srv2_slope_m  .4
  srv2_sel50_m  8.
  srv3_slope_f  .4
  srv3_sel50_f  8.
  srv3_slope_m  .4
  srv3_sel50_m  8.
  alpha 1.
  beta 0.

PARAMETER_SECTION
 //parameters to be estimated are all ones that begin with init_ and have a positive
 //phase, negative phase means are fixed.
 //phase of 8 is greater than last phase so does q1 in last phase  
  // init_bounded_number q1(.5,2,8)
 //fix q1 to be 1 otherwise it went to lower bound of .5
  init_bounded_number q1(0.5,2.0,-4)
  init_bounded_number q2(0.05,1.5,-4)
  init_bounded_number q3(0.05,1.5,-4)
  init_number alpha(4)       // used to estimate temperature effect on shelf survey catchability
  init_number beta(4)  // used to estimate temperature effect on shelf survey catchability
 //phase of -1 means M is fixed   
 // init_bounded_vector M(1,2,.02,.8,-3)
  init_number mean_log_rec(1)
  // init_bounded_dev_vector rec_dev(styr_rec,endyr-2,-15,15,2)
  init_bounded_dev_vector rec_dev(styr_rec,endyr-1,-15,15,2) //JNI
//  init_vector rec_dev_future(styr_fut,endyr_fut,phase_F40);
  init_number log_avg_fmort(2)
  !! cout <<"log_avg_fmort"<<endl;
  !! cout <<log_avg_fmort<<endl;
  init_bounded_dev_vector fmort_dev(styr,endyr,-3,3,1)
 
//  Selectivity parameters from the GOA version of the model

  init_matrix log_selcoffs_fish(1,2,1,nselages,phase_selcoffs)
  init_bounded_number fish_slope_f(.1,5.,phase_logistic_sel)
  init_bounded_number fish_sel50_f(1.,15.,phase_logistic_sel)
  init_bounded_number fish_slope_m(.05,.8,phase_logistic_sel)
  init_bounded_number fish_sel50_m(1.,25.,phase_logistic_sel)
  init_bounded_number srv1_slope_f1(.1,5.,phase_logistic_sel)
  init_bounded_number srv1_sel50_f1(1.,10.,phase_logistic_sel)
  init_bounded_number srv1_slope_f2(.1,5.,phase_logistic_sel)
  init_bounded_number srv1_sel50_f2(1.,10.,phase_logistic_sel)

  init_bounded_number srv1_slope_m1(.01,.5,phase_logistic_sel)
  init_bounded_number srv1_sel50_m1(1.,12.,phase_logistic_sel)
  init_bounded_number srv1_slope_m2(.01,.5,phase_logistic_sel)
  init_bounded_number srv1_sel50_m2(1.,12.,phase_logistic_sel)

  init_bounded_number srv2_slope_f(.1,5.,phase_logistic_sel)
  init_bounded_number srv2_sel50_f(1.,10.,phase_logistic_sel)
  init_bounded_number srv2_slope_m(.01,.5,phase_logistic_sel)
  init_bounded_number srv2_sel50_m(1.,12.,phase_logistic_sel)

  init_bounded_number srv3_slope_f(.1,5.,phase_logistic_sel)
  init_bounded_number srv3_sel50_f(1.,10.,phase_logistic_sel)
  init_bounded_number srv3_slope_m(.01,.5,phase_logistic_sel)
  init_bounded_number srv3_sel50_m(1.,12.,phase_logistic_sel)
  init_bounded_number sexr_param_fish(1.0,1.0,-5)  //this was hitting bound of 1.0 so fixed it - should free up to check

// Parameters for computing SPR rates 
  init_bounded_number F40(0.01,1.,phase_F40)
  init_bounded_number F35(0.01,1.,phase_F40)
  init_bounded_number F30(0.01,1.,phase_F40)

  matrix log_sel_fish(1,2,1,nages)
  matrix sel(1,2,1,nages) //fishery selectivity
  matrix sel_srv1(1,2,1,nages)
  matrix sel_srv2(1,2,1,nages)
  matrix sel_srv3(1,2,1,nages)
  vector avgsel_fish(1,2)
  matrix popn(1,2,styr,endyr)
  matrix totn_srv1(1,2,styr,endyr)
  matrix totn_srv2(1,2,styr,endyr)
  matrix totn_srv3(1,2,styr,endyr)
  vector M(1,2)
  vector temp1(1,nages)
  vector temp2(1,nages)
  vector explbiom(styr,endyr)
  vector pred_bio(styr,endyr)
  vector fspbio(styr,endyr)
  vector pred_srv1(styr,endyr)
  vector pred_srv2(styr,endyr)
  vector pred_srv3(styr,endyr)
  3darray pred_p_fish(1,2,styr,endyr,1,nlen)
  3darray pred_p_srv1_age(1,2,1,nobs_srv1_age,1,nages)
  3darray pred_p_srv3_age(1,2,1,nobs_srv3_age,1,nages) 
  3darray pred_p_srv1_len(1,2,1,nobs_srv1_length,1,nlen)
  3darray pred_p_srv2_len(1,2,1,nobs_srv2_length,1,nlen)
  3darray pred_p_srv3_len(1,2,1,nobs_srv3_length,1,nlen)
  vector pred_catch(styr,endyr)
  3darray natage(1,2,styr,endyr,1,nages) 
  sdreport_vector totalbiomass(styr,endyr)
  3darray catage(1,2,styr,endyr,1,nages)
  3darray Z(1,2,styr,endyr,1,nages)
  3darray F(1,2,styr,endyr,1,nages)
  3darray S(1,2,styr,endyr,1,nages)
  vector fmort(styr,endyr)
  number rbar
  vector surv(1,2)  //survival for each sex
  vector offset(1,7)
  number rec_like
  number catch_like
  number sexr_like
  vector age_like(1,4) //really only need shelf and AI but may need other elements later
  vector length_like(1,4) 
  vector sel_like(1,4)
  number fpen    
  number surv1_like
  number surv2_like
  number surv3_like
  sdreport_number endbiom
  sdreport_number depletion
  objective_function_value obj_fun
  number tmp
  vector pred_sexr(styr,endyr)
 // Stuff for SPR and yield projections
  number sigmar
  number ftmp
  number SB0
  number SBF40
  number SBF35
  number SBF30
  number sprpen
  matrix Nspr(1,4,1,nages)
  3darray nage_future(1,2,styr_fut,endyr_fut,1,nages)
  sdreport_matrix fspbiom_fut(1,4,styr_fut,endyr_fut)
  3darray F_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray Z_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray S_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray catage_future(1,2,styr_fut,endyr_fut,1,nages)
  number avg_rec_dev_future
  vector avg_F_future(1,4)
  sdreport_matrix catch_future(1,3,styr_fut,endyr_fut) // Note, don't project for F=0 (it 
  sdreport_matrix future_biomass(1,4,styr_fut,endyr_fut)
  vector explbiom_fut(styr_fut,endyr_fut)
  number maxsel_fish
  number maxsel_srv1
  number maxsel_srv2
  number maxsel_srv3
  number mlike
  number qlike
  number flike
  vector qtime(styr,endyr)

PRELIMINARY_CALCS_SECTION  
  test(1,1,2)=6;
  test(1,2)=5;
  test(2,1)=3;

  cout <<"test"<<endl; 
  cout <<test<<endl;
  obs_mean_sexr=0.34;  //initial value for avg proportion of male population estimated from shelf surveys; calculated below
  obs_SD_sexr=0.0485;  //initial value for standard deviation of mean male population proportion: calculated below
//compute sex ratio in catch: proportion of females?
  for(i=1; i<=nobs_fish;i++)
  {
    obs_sexr(i) = sum(obs_p_fish(1,i))/sum(obs_p_fish(1,i) + obs_p_fish(2,i)); 
  }

//length obs sex ratio in surveys: proportion of males
  for(i=1; i<=nobs_srv1_length;i++)
    obs_sexr_srv1_2(i) = (sum(obs_p_srv1_length(2,i)))/
                         (sum(obs_p_srv1_length(1,i)) + sum(obs_p_srv1_length(2,i)));
    obs_mean_sexr=mean(obs_sexr_srv1_2);  
  cout <<"obs_sexr_srv1_2"<<endl; 
  cout <<obs_sexr_srv1_2<<endl;
    obs_SD_sexr=std_dev(obs_sexr_srv1_2);

  for(i=1; i<=nobs_srv2_length;i++)
    obs_sexr_srv2_2(i) = (sum(obs_p_srv2_length(2,i)))/
                         (sum(obs_p_srv2_length(1,i)) + sum(obs_p_srv2_length(2,i))); 

  for(i=1; i<=nobs_srv3_length;i++)
    obs_sexr_srv3_2(i) = (sum(obs_p_srv3_length(2,i)))/
                         (sum(obs_p_srv3_length(1,i)) + sum(obs_p_srv3_length(2,i))); 

 //Compute offset for multinomial and length bin proportions
 // offset is a constant nplog(p) is added to the likelihood     
 // magnitude depends on nsamples(sample size) and p's_
  //k is sex loop
  offset.initialize(); 
  for (i=1; i <= nobs_fish; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i));
    obs_p_fish(1,i) = obs_p_fish(1,i) / sumtot; 
    obs_p_fish(2,i) = obs_p_fish(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(1) -= nsamples_fish(k,i)*obs_p_fish(k,i) * log(obs_p_fish(k,i)+.0001);
  }
 //shelf survey length offset and bin proportions
  for (i=1; i <= nobs_srv1_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv1_length(1,i)+obs_p_srv1_length(2,i));
    obs_p_srv1_length(1,i) = obs_p_srv1_length(1,i) / sumtot; 
    obs_p_srv1_length(2,i) = obs_p_srv1_length(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(2) -= nsamples_srv1_length(k,i)*obs_p_srv1_length(k,i) * log(obs_p_srv1_length(k,i)+.0001);
  }

//slope survey length offset and bin proportions
  for (i=1; i <= nobs_srv2_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv2_length(1,i)+obs_p_srv2_length(2,i));
    obs_p_srv2_length(1,i) = obs_p_srv2_length(1,i) / sumtot; 
    obs_p_srv2_length(2,i) = obs_p_srv2_length(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(3) -= nsamples_srv2_length(k,i)*obs_p_srv2_length(k,i) * log(obs_p_srv2_length(k,i)+.0001);
  }

//Aleutian Islands survey length offset and bin proportions
  for (i=1; i <= nobs_srv3_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv3_length(1,i)+obs_p_srv3_length(2,i));
    obs_p_srv3_length(1,i) = obs_p_srv3_length(1,i) / sumtot; 
    obs_p_srv3_length(2,i) = obs_p_srv3_length(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(4) -= nsamples_srv3_length(k,i)*obs_p_srv3_length(k,i) * log(obs_p_srv3_length(k,i)+.0001);
  }

//shelf survey age offset 
  for (i=1; i <= nobs_srv1_age; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv1_age(1,i)+obs_p_srv1_age(2,i));
    obs_p_srv1_age(1,i) = obs_p_srv1_age(1,i) / sumtot; 
    obs_p_srv1_age(2,i) = obs_p_srv1_age(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(5) -= nsamples_srv1_age(k,i)*obs_p_srv1_age(k,i) * log(obs_p_srv1_age(k,i)+.0001);
  }   
  cout<<"obs_p_srv3_age"<<endl;
  cout<<obs_p_srv3_age<<endl;
//AI survey age offset 
  for (i=1; i <= nobs_srv3_age; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv3_age(1,i)+obs_p_srv3_age(2,i));
    obs_p_srv3_age(1,i) = obs_p_srv3_age(1,i) / sumtot; 
    obs_p_srv3_age(2,i) = obs_p_srv3_age(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(6) -= nsamples_srv3_age(k,i)*obs_p_srv3_age(k,i) * log(obs_p_srv3_age(k,i)+.0001);
  }
  cout<<"obs_p_srv3_age"<<endl;
  cout<<obs_p_srv3_age<<endl;  
  cout<<"offset"<<endl;
  cout<<offset<<endl;  
  cout<<"offset(6)"<<endl;
  cout<<offset(6)<<endl;
  M(1)=0.20;
  M(2)=0.35;

PROCEDURE_SECTION
//this is for bootstraping where qrun is a vector of q's from bootstrap irun is the 
//run number.  sets the q (=q1) for the run.
//   q1=qrun(irun);
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
      // loop over years and print in one long row.
      for (i=styr;i<=endyr;i++)
      evalout<<  fspbio(i) << " ";
      for (i=styr;i<=endyr;i++)    
      evalout<< natage(1,i)*wt(1) + natage(2,i)*wt(2) <<" ";
      for (i=styr;i<=endyr;i++)
      evalout << 2*natage(1,i,1) <<" ";
      // hit carriage return on file
      evalout <<  endl;
    }
  }

    evaluate_the_objective_function();

   
FUNCTION get_selectivity
//fishery selectivities
  if(active(log_selcoffs_fish))// init_matrix log_selcoffs_fish(1,2,1,nselages,phase_selcoffs) set to phase 4  
 {
//turn off logistic curve
//   phase_logistic_sel=-2;  //its already set to -4 in beginning
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
 //cout<<"calc log_sel_fish"<<endl;
  for(k=1;k<=2;k++)
    {
      log_sel_fish(k)-=log(mean(mfexp(log_sel_fish(k))));
      sel(k)=mfexp(log_sel_fish(k));
      if(k==2)
       {
         sel(k)=sel(k)*sexr_param_fish; //fixed at 1 in GOA model not BSAI model
       }
      //cout<<"sel survey"<<sel_srv1<<endl;
    } 
    for (j=1;j<=nages;j++)  //this is selectivity for the surveys
    { 
      
      //ascending limb of curve for shelf survey
		sel_srv1(1,j)=1./(1.+mfexp(-1.*srv1_slope_f1*(double(j)-srv1_sel50_f1)));
		sel_srv1(2,j)=1./(1.+mfexp(-1.*srv1_slope_m1*(double(j)-srv1_sel50_m1)));
		//decending limb of curve for shelf survey 
		temp1=1./(1.+mfexp(srv1_slope_f2*(double(j)-srv1_sel50_f2)));
		temp2=1./(1.+mfexp(srv1_slope_m2*(double(j)-srv1_sel50_m2)));
		sel_srv1(1,j)=sel_srv1(1,j)*temp1(j);
		sel_srv1(2,j)=sel_srv1(2,j)*temp2(j);
		//slope surveys
		sel_srv2(1,j)=1./(1.+mfexp(-1.*srv2_slope_f*(double(j)-srv2_sel50_f)));
		sel_srv2(2,j)=1./(1.+mfexp(-1.*srv2_slope_m*(double(j)-srv2_sel50_m)));
		//Aleutian Islands surveys
		sel_srv3(1,j) = 1./(1.+mfexp(-1.*srv3_slope_f*(double(j)-srv3_sel50_f)));
		sel_srv3(2,j) = 1./(1.+mfexp(-1.*srv3_slope_m*(double(j)-srv3_sel50_m)));							
    }
 }//  end if(active(log_selcoffs_fish))
  else
    {
     //logistic selectivity curve
          for (j=1;j<=nages;j++)
          { 
            if(j<=nselages)
             {
               sel(1,j)=1./(1.+mfexp(-1.*fish_slope_f*(double(j)-fish_sel50_f)));
               sel(2,j)=1./(1.+mfexp(-1.*fish_slope_m*(double(j)-fish_sel50_m)));
             }
            else
            {
             sel(1,j)=sel(1,j-1);
             sel(2,j)=sel(2,j-1);
            }
            //ascending limb of curve for shelf survey
			sel_srv1(1,j)=1./(1.+mfexp(-1.*srv1_slope_f1*(double(j)-srv1_sel50_f1)));
			sel_srv1(2,j)=1./(1.+mfexp(-1.*srv1_slope_m1*(double(j)-srv1_sel50_m1)));
			//decending limb of curve for shelf survey 
			temp1=1./(1.+mfexp(srv1_slope_f2*(double(j)-srv1_sel50_f2)));
			temp2=1./(1.+mfexp(srv1_slope_m2*(double(j)-srv1_sel50_m2)));
			sel_srv1(1,j)=sel_srv1(1,j)*temp1(j);
			sel_srv1(2,j)=sel_srv1(2,j)*temp2(j);
			//slope surveys
			sel_srv2(1,j)=1./(1.+mfexp(-1.*srv2_slope_f*(double(j)-srv2_sel50_f)));
			sel_srv2(2,j)=1./(1.+mfexp(-1.*srv2_slope_m*(double(j)-srv2_sel50_m)));
			//Aleutian Islands surveys
			sel_srv3(1,j) = 1./(1.+mfexp(-1.*srv3_slope_f*(double(j)-srv3_sel50_f)));
			sel_srv3(2,j) = 1./(1.+mfexp(-1.*srv3_slope_m*(double(j)-srv3_sel50_m)));							
          } 
     }
//  cout<<"end new sel here"<<endl; 

//     logistic selectivity curves, asymptotic for fishery, slope survey and the Aleutian Islands but domed shape for shelf survey  
 
FUNCTION get_mortality 
  maxsel_fish=max(sel(1));     //1 is females
  if(maxsel_fish<max(sel(2)))  //if highest female selectivity is > male selectivity, make maxsel_fish=male high selectivity
      maxsel_fish=max(sel(2));

  fmort = mfexp( log_avg_fmort+fmort_dev); 
  for(k=1;k<=2;k++)
  {
    for (i=styr;i<=endyr;i++)
    {
      F(k,i)=(sel(k))*fmort(i);
      Z(k,i)=F(k,i) + M(k); 
    }
  } 
  S = mfexp(-1.*Z); 

FUNCTION get_numbers_at_age
  maxsel_fish=max(sel(1));   
  if(maxsel_fish<max(sel(2)))//if females greater than males, then set the max to the females.
    maxsel_fish=max(sel(2)); //set max to whichever sex is larger

  maxsel_srv1=max(sel_srv1(1));
  if(maxsel_srv1<max(sel_srv1(2)))
    maxsel_srv1=max(sel_srv1(2)); 

  maxsel_srv2=max(sel_srv2(1));
  if(maxsel_srv2<max(sel_srv2(2)))
    maxsel_srv2=max(sel_srv2(2)); 
  
  maxsel_srv3=max(sel_srv3(1));
  if(maxsel_srv3<max(sel_srv3(2)))
    maxsel_srv3=max(sel_srv3(2)); 

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
      //subvector - avoids writing a j loop  =++ increments the right side 
      //(1,nages-1) to 1+1 to nages-1+1 then does the assignment x(i)(1,n) 
      //takes the ith row of x the columns 1 to n
      //      natage(k,i+1)(2,nages)=++elem_prod(natage(k,i)(1,nages-1),S(k,i)(1,nages-1));
      for(j=1;j<nages;j++)
      {
        natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
      }
      //accumulates oldest ages
      // cout<<"done with j loop"<<endl;
      natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
      // cout<<"done with natage nages"<<endl;
      //popn is exploitable numbers
      popn(k,i)= natage(k,i)*sel(k);
      // cout<<"popn "<<endl;
      // cout<<popn(k,i)<<endl;
    }
    // cout<<"to popn"<<endl; 
    popn(k,endyr)=natage(k,endyr)*sel(k);
  }
  for (i=styr;i<=endyr;i++)
  {
      pred_sexr(i)=sum(natage(2,i))/(sum((natage(1,i)+natage(2,i))));  //calculation of prop. of males in pred. population 
    for(k=1;k<=2;k++)
    {
      totn_srv1(k,i)=q1*(natage(k,i)*sel_srv1(k)); // not used in further calculations
      totn_srv2(k,i)=q2*(natage(k,i)*sel_srv2(k)); // not used in further calculations        
    }
  }   

  //predicted survey values
  fspbio.initialize(); 
  qtime=q1;
  for (i=styr;i<=endyr;i++)
  {
    fspbio(i) = natage(1,i)*elem_prod(wt(1),maturity);
    explbiom(i)=0.;
    pred_bio(i)=0.;
    pred_srv1(i)=0.;
    pred_srv2(i)=0.;
    pred_srv3(i)=0.; //JNI

    // if (i>=1982 )      //catchability calculation for survey years
    if (i>=1982 && i-1981 <= nobs_srv1)      //JNI catchability calculation for survey years
    
    qtime(i)=q1*mfexp(-alpha+beta*bottom_temps(i-1981));
    

    for(k=1;k<=2;k++)
    {
      
      pred_srv1(i) += qtime(i)*(natage(k,i)*elem_prod(sel_srv1(k),wt(k)))/maxsel_srv1;   //shelf survey, dividing by the maxsel constrains female selectivity to be 1.0
      //pred_srv1(i) += q1*(natage(k,i)*elem_prod(sel_srv1(k),wt(k)))/maxsel_srv1;   //shelf survey, without temperature q modeling
      pred_srv2(i) += q2*(natage(k,i)*elem_prod(sel_srv2(k),wt(k))); // /maxsel_srv2;         //slope survey JNI
      pred_srv3(i) += q3*(natage(k,i)*elem_prod(sel_srv3(k),wt(k))); // /maxsel_srv3;         //Aleutian Islands survey JNI

      //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
      // are set different in the tpl file the program will take to value from the bin file and use that 
      //   pred_srv1(i)=1.0*(natage(i)*elem_prod(sel_srv1,wt));
      explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
      pred_bio(i)+=natage(k,i)*wt(k);
    }
  }
      // cout <<q3<<endl<<pred_srv3<<endl;exit(1);
    //don't need to divide by max_sel because totn_srv1 is calculated using selectivities and the
    //max_sel would cancel out.

    // Fitting the survey length compositions
    for(i=1; i<=nobs_srv1_length;i++)
    {
      ii = yrs_srv1_length(i);
      pred_p_srv1_len(1,i) = q1 * elem_prod(sel_srv1(1),natage(1,ii)) * lenage(1);
      pred_p_srv1_len(2,i) = q1 * elem_prod(sel_srv1(2),natage(2,ii)) * lenage(2);
      dvariable sum_tot = sum(pred_p_srv1_len(1,i)+pred_p_srv1_len(2,i));
      pred_p_srv1_len(1,i) /= sum_tot;
      pred_p_srv1_len(2,i) /= sum_tot;
    }
   
    for(i=1; i<=nobs_srv2_length;i++)
    {
      ii = yrs_srv2_length(i);
      pred_p_srv2_len(1,i)=q2*elem_prod(sel_srv2(1),natage(1,ii))*lenage(1);
      pred_p_srv2_len(2,i)=q2*elem_prod(sel_srv2(2),natage(2,ii))*lenage(2);
      dvariable sum_tot = sum(pred_p_srv2_len(1,i)+pred_p_srv2_len(2,i));
      pred_p_srv2_len(1,i) /= sum_tot;
      pred_p_srv2_len(2,i) /= sum_tot;
    }
    for(i=1; i<=nobs_srv3_length;i++)
    {
      ii = yrs_srv3_length(i);
      pred_p_srv3_len(1,i)=q3*elem_prod(sel_srv3(1),natage(1,ii))*lenage(1);
      pred_p_srv3_len(2,i)=q3*elem_prod(sel_srv3(2),natage(2,ii))*lenage(2);
      dvariable sum_tot = sum(pred_p_srv3_len(1,i)+pred_p_srv3_len(2,i));
      pred_p_srv3_len(1,i) /= sum_tot;
      pred_p_srv3_len(2,i) /= sum_tot;
    }

    //Calculation of survey age composition

    for(i=1; i<=nobs_srv1_age;i++)
    {
      ii = yrs_srv1_age(i);
      pred_p_srv1_age(1,i) = q1 * elem_prod(sel_srv1(1),natage(1,ii));
      pred_p_srv1_age(2,i) = q1 * elem_prod(sel_srv1(2),natage(2,ii));
      dvariable sum_tot = sum(pred_p_srv1_age(1,i)+pred_p_srv1_age(2,i));
      pred_p_srv1_age(1,i) /= sum_tot;
      pred_p_srv1_age(2,i) /= sum_tot;
    } 

    for(i=1; i<=nobs_srv3_age;i++)  //LOOK BACK EHRE
    {
      ii = yrs_srv3_age(i);
      pred_p_srv3_age(1,i) = q3 * elem_prod(sel_srv3(1),natage(1,ii));
      pred_p_srv3_age(2,i) = q3 * elem_prod(sel_srv3(2),natage(2,ii));
      dvariable sum_tot = sum(pred_p_srv3_age(1,i)+pred_p_srv3_age(2,i));
      pred_p_srv3_age(1,i) /= sum_tot;
      pred_p_srv3_age(2,i) /= sum_tot;
    }

  depletion=pred_bio(endyr)/pred_bio(styr);
  endbiom=pred_bio(endyr);

FUNCTION get_catch_at_age
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

FUNCTION Future_projections
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

      // Get future F's
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
    // Future Recruitment (and spawners)
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
  
FUNCTION compute_spr_rates
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
  //cout<<F40<<" "<<F30<<" "<<Nspr<<endl; 
 // cout<<"spr calc"<<endl;
 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*exp(-1.*M(1))/(1.-exp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*exp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*exp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
  Nspr(4,nages)=Nspr(4,nages-1)*exp(-1.* (M(1)+F30*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F30*sel(1,nages)/maxsel_fish)));
 //cout<<"plus group"<<endl;
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

FUNCTION Do_depend
  for (i=styr;  i<=endyr;  i++) 
  totalbiomass(i)=natage(1,i)*wt(1) + natage(2,i)*wt(2);
  obj_fun += 1.*sexr_like;             // male proportion prior, emphasis factor = 1

FUNCTION evaluate_the_objective_function
  length_like.initialize();
  age_like.initialize();
  fpen.initialize();
  rec_like.initialize();
  surv1_like.initialize();
  surv2_like.initialize();
  surv3_like.initialize();
  catch_like.initialize();
  sexr_like.initialize();
  obj_fun.initialize();

  if (active(rec_dev))
  {
    length_like.initialize();   //length-like vector has the likelihoods for the 4 components: 1) fishery length 2) shelf survey lengths 3) slope survey lengths 4) Aleutians
    int ii;

    //recruitment likelihood - norm2 is sum of square values   
    rec_like = norm2(rec_dev);

    for(k=1;k<=2;k++)
    {
      for (i=1; i <= nobs_fish; i++)
      {
        ii=yrs_fish(i);
        //fishery length likelihood fitting
          length_like(1) -= nsamples_fish(k,i)*(1e-5+obs_p_fish(k,i))*log(pred_p_fish(k,ii)+1e-5);
      }
    }
    //add the offset to the likelihood   
    length_like(1)-=offset(1);

    //shelf survey length composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv1_length; i++)
        length_like(2)-=nsamples_srv1_length(k,i)*(1e-3+obs_p_srv1_length(k,i))*log(pred_p_srv1_len(k,i)+1e-3);
    length_like(2)-=offset(2);

    //slope survey length composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv2_length; i++)
        length_like(3)-=nsamples_srv2_length(k,i)*(1e-3+obs_p_srv2_length(k,i))*log(pred_p_srv2_len(k,i)+1e-3);
    length_like(3)-=offset(3);

    //Aleutian Island survey length composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv3_length; i++)
        length_like(4)-=nsamples_srv3_length(k,i)*(1e-3+obs_p_srv3_length(k,i))*log(pred_p_srv3_len(k,i)+1e-3);
    length_like(4)-=offset(4);

   //shelf survey age composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv1_age; i++)
        age_like(2)-=nsamples_srv1_age(k,i)*(1e-3+obs_p_srv1_age(k,i))*log(pred_p_srv1_age(k,i)+1e-3);
    age_like(2)-=offset(5);
  
   //AI survey age composition fitting   LOOK BACK HERE
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv3_age; i++)
        age_like(4)-=nsamples_srv3_age(k,i)*(1e-3+obs_p_srv3_age(k,i))*log(pred_p_srv3_age(k,i)+1e-3);
    age_like(4)-=offset(6);

  //end of if(active (rec_dev))
  }
  // Fit to indices (lognormal) 
  //weight each years estimate by 1/(2*variance) - use cv as an approx to s.d. of log(biomass) 

   surv1_like = norm2(elem_div(log(obs_srv1+.01)-log(pred_srv1(yrs_srv1)+.01),sqrt(2)*cv_srv1));
   surv2_like = norm2(elem_div(log(obs_srv2+.01)-log(pred_srv2(yrs_srv2)+.01),sqrt(2)*cv_srv2));
   surv3_like = norm2(elem_div(log(obs_srv3+.01)-log(pred_srv3(yrs_srv3)+.01),sqrt(2)*cv_srv3));
   double var_tmp; for (i=1;i<=nobs_srv3;i++) { var_tmp = 2.*square(log(obs_srv3(i))*cv_srv3(i)); surv3_like += square(log(obs_srv3(i)+.01)-log(pred_srv3(yrs_srv3(i))+.01))/var_tmp; }
   // surv_like = norm2(log(obs_srv1+.01)-log(pred_srv1(yrs_srv1)+.01));

    catch_like=norm2(log(catch_bio+.000001)-log(pred_catch+.000001));

   // sex ratio likelihood   note only calculated based on survey 1
     sexr_like=0.5*norm2((obs_mean_sexr-pred_sexr)/obs_SD_sexr); 
 //selectivity likelihood is penalty on how smooth selectivities are   
 //here are taking the sum of squares of the second differences 
  if(active(log_selcoffs_fish))
  {  
    sel_like(1)=wt_like(1)*norm2(first_difference(first_difference(log_sel_fish(1))));
//    sel_like(2)=wt_like(2)*norm2(first_difference(first_difference(log_sel_srv1(1))));
    sel_like(3)=wt_like(3)*norm2(first_difference(first_difference(log_sel_fish(2))));
//    sel_like(4)=wt_like(4)*norm2(first_difference(first_difference(log_sel_srv1(2)))); 
   for (j=1;j<nages;j++)
   {
    if(monot_sel==1)
    { 
        if (log_sel_fish(1,j)>log_sel_fish(1,j+1))
        sel_like(1)+=wt_like(5)*square(log_sel_fish(1,j)-log_sel_fish(1,j+1));
        if (log_sel_fish(2,j)>log_sel_fish(2,j+1))
        sel_like(3)+=wt_like(6)*square(log_sel_fish(2,j)-log_sel_fish(2,j+1));
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
  obj_fun += 1.*length_like(1);     //emphasis factor = 1 for fishery lengths   
  obj_fun += 1.*length_like(2);     //emphasis factor = 1
  obj_fun += 1.*length_like(3);     //emphasis factor = 1
  obj_fun += 1.*length_like(4);     //emphasis factor = 1
  obj_fun += 1.*age_like(2);           //emphasis factor = 1 
  obj_fun += 1.*age_like(4);           //emphasis factor = 1 
  obj_fun += 1.*surv1_like;         //emphasis factor = 1
  obj_fun += 1.*surv2_like;         //emphasis factor = 1
  obj_fun += 1.*surv3_like;         //emphasis factor = 1
  obj_fun += 300*catch_like;        // large emphasis to fit observed catch
  obj_fun += fpen;
  obj_fun += sprpen;

REPORT_SECTION
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

  report << "Estimated shelf survey selectivity for females at age " << endl;
    for (j=1; j<=nages;j++)
      report << j<<" " <<sel_srv1(1,j)/maxsel_srv1<< endl;
   
  report << "Estimated shelf survey selectivity for males at age " << endl;
    for (j=1; j<=nages;j++)
      report <<  j<<" "  <<sel_srv1(2,j)/maxsel_srv1<<endl;

  report << "Estimated slope survey selectivity for females at age " << endl;
    for (j=1; j<=nages;j++)
      report << j <<" " <<sel_srv2(1,j)<< endl;
   
  report << "Estimated slope survey selectivity for males at age " << endl;
    for (j=1; j<=nages;j++)
      report << j<<" "  <<sel_srv2(2,j)<<endl;


  report << "Estimated Aleutian Islands survey selectivity for females at age " << endl;
    for (j=1; j<=nages;j++)
      report <<  j<<" "  <<sel_srv3(1,j)<< endl;
   
  report << "Estimated Aleutian Islands survey selectivity for males at age " << endl;
    for (j=1; j<=nages;j++)
      report << j<<" "  <<sel_srv3(2,j)<<endl;


  report << endl << "Bering Sea shelf survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nobs_srv1;i++)
      report << yrs_srv1(i) << ","<< obs_srv1(i) << "," << pred_srv1(yrs_srv1(i)) << endl;

  report << "Bering Sea slope survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nobs_srv2;i++)
      report << yrs_srv2(i) << ","<< obs_srv2(i) << "," << pred_srv2(yrs_srv2(i)) << endl;

  report << "Aleutian Islands survey biomass (Year, Obs_biomass, Pred_biomass) "  << endl;
    for (i=1; i<=nobs_srv3;i++)
      report << yrs_srv3(i) << ","<< obs_srv3(i) << "," << pred_srv3(yrs_srv3(i)) << endl;

  report <<" Observed female shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << obs_p_srv1_length(1,i) << endl;

  report <<" Predicted female shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << pred_p_srv1_len(1,i) << endl;

  report <<" Observed male shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << obs_p_srv1_length(2,i) << endl;

  report <<" Predicted male shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i)  << pred_p_srv1_len(2,i) << endl;
  
  report <<" Observed female slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i) << obs_p_srv2_length(1,i) << endl;

  report <<" Predicted female slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i)  << pred_p_srv2_len(1,i) << endl;

  report <<" Observed male slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i) << obs_p_srv2_length(2,i) << endl;
  
  report <<" Predicted male slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i)  << pred_p_srv2_len(2,i) << endl;
  
  report <<" Observed female Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i) << obs_p_srv3_length(1,i) << endl;

  report <<" Predicted female Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i)  << pred_p_srv3_len(1,i) << endl;

  report <<" Observed male Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i) << obs_p_srv3_length(2,i) << endl;
  
  report <<" Predicted male Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i)  << pred_p_srv3_len(2,i) << endl;

  report <<" Observed female shelf survey age composition " << endl;
    for (i=1; i<=nobs_srv1_age; i++)
      report << yrs_srv1_age(i) << obs_p_srv1_age(1,i) << endl;

  report <<" Predicted female shelf survey age composition " << endl;
    for (i=1; i<=nobs_srv1_age; i++)
      report << yrs_srv1_age(i)  << pred_p_srv1_age(1,i) << endl;

  report <<" Observed male shelf survey age composition " << endl;
    for (i=1; i<=nobs_srv1_age; i++)
      report << yrs_srv1_age(i) << obs_p_srv1_age(2,i) << endl;

  report <<" Predicted male shelf survey age composition " << endl;
    for (i=1; i<=nobs_srv1_age; i++)
      report << yrs_srv1_age(i)  << pred_p_srv1_age(2,i) << endl;

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
  report << "shelf survey like component " << endl;
  report << surv1_like << endl;
  report << "slope survey like component " << endl;
  report << surv2_like << endl;
  report << "Aleutian Islands survey lilke component "<< endl;
  report <<surv3_like << endl;
  report << "shelf survey length composition " << endl;
  report << length_like(2) << endl;
  report << "slope survey length composition " << endl;
  report << length_like(3) << endl;
  report << "Aleutian Islands survey length composition " << endl;
  report << length_like(4) << endl;
  report << "fishery length composition likelihood " << endl;
  report << length_like(1) << endl;
  report << "recruitment likelihood component est.  " << endl;
  report << rec_like << endl;
  report << "catch likelihood component est.  " << endl;
  report << catch_like << endl;
  report << "sex ratio likelihood component " << endl;
  report << sexr_like << endl;
  report << "shelf survey age composition  " << endl;
  report << age_like << endl;
  
  report << "Projected biomass" << endl;
  report << future_biomass << endl;
  report <<"projected future female spawning biomass " << endl;
  report <<fspbiom_fut << endl;
  report << "Future yield " << endl;
  report << catch_future << endl;
  report << "shelf survey q =" << endl;
  report << q1 << endl;
  report << "slope survey q = " << endl;
  report << q2 << endl;
  report << "Aleutian Islands survey q = " << endl;
  report << q3 << endl;
  report << " female natural mortality for this run" << endl;
  report << M(1) << endl;
  report << " male natural mortality for this run" << endl;
  report << M(2) << endl;

  report <<endl << "temperature effect (q) for the shelf survey "<< endl;
   for (i=1;i<=nobs_srv1;i++)
     report <<yrs_srv1(i)<<","<<bottom_temps(i)<<","<<qtime(yrs_srv1(i))<<endl;

  report << "predicted male proportion in population" << endl;
   for (i=styr;i<=endyr;i++)
      report << i << " " << pred_sexr(i) << endl;

  report << "mean observed prop. male in shelf surveys = "<< endl;
  report << obs_mean_sexr << endl;
  report << "stdev of mean observed prop. male in shelf surveys = " << endl;
  report << obs_SD_sexr << endl;

  report <<"alpha = "<< endl;
  report <<alpha<<endl;
  report << "beta= "<< endl;
  report  << beta << endl;

  report << "standard error of biomass in shelf surveys = " << endl;
  report << obs_srv1_sd << endl;
  report << "standard error of biomass in slope surveys = " << endl;
  report << obs_srv2_sd << endl;
  report << "standard error of biomass in AI surveys = " << endl;
  report << obs_srv3_sd << endl;
  report << " recruits" << endl;
    for (i=styr;  i<=endyr;  i++)
         report  << i <<" " << 2*natage(1,i,1) <<" "<< endl;   

	  report <<" Observed male Aleutians survey age composition " << endl;
	    for (i=1; i<=nobs_srv3_age; i++)
	      report << yrs_srv3_age(i) << obs_p_srv3_age(2,i) << endl;

	  report <<" Predicted male Aleutians survey age composition " << endl;
	    for (i=1; i<=nobs_srv3_age; i++)
	      report << yrs_srv3_age(i)  << pred_p_srv3_age(2,i) << endl;
	
		  report <<" Observed female Aleutians survey age composition " << endl;
		    for (i=1; i<=nobs_srv3_age; i++)
		      report << yrs_srv3_age(i) << obs_p_srv3_age(1,i) << endl;

		  report <<" Predicted female Aleutians survey age composition " << endl;
		    for (i=1; i<=nobs_srv3_age; i++)
		      report << yrs_srv3_age(i)  << pred_p_srv3_age(1,i) << endl;
	
  report << " Go drink coffee " << endl; 

  report << "SARA form for Angie Grieg" << endl;

  report << "ATF        # stock  " << endl;
  report << "BSAI       # region     (AI AK BOG BSAI EBS GOA SEO WCWYK)" << endl;
  report << "2013       # ASSESS_YEAR - year assessment is presented to the SSC" << endl;
  report << "3a         # TIER  (1a 1b 2a 2b 3a 3b 4 5 6) " << endl;
  report << "none       # TIER2  if mixed (none 1a 1b 2a 2b 3a 3b 4 5 6)" << endl;
  report << "partial    # UPDATE (new benchmark full partial)" << endl;
  report << "2          # LIFE_HIST - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "2          # ASSES_FREQ - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "5          # ASSES_LEV - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "5          # CATCH_DAT - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "3          # ABUND_DAT - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "567000     # Minimum B  Lower 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "665000     # Maximum B  Upper 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "202138     # BMSY  is equilibrium spawning biomass at MSY (Tiers 1-2) or 7/8 x B40% (Tier 3)" << endl;
  report << "ADMB       # MODEL - Required only if NMFS toolbox software used; optional otherwise " << endl;
  report << "           # VERSION - Required only if NMFS toolbox software used; optional otherwise" << endl;
  report << "2          # number of sexes  if 1 sex=ALL elseif 2 sex=(FEMALE, MALE) " << endl;
  report << "1          # number of fisheries" << endl;
  report << "1          # multiplier for recruitment, N at age, and survey number (1,1000,1000000)" << endl;
  report << "1          # recruitment age used by model or size" << endl;
  report << "1          # age+ or mmCW+ used for biomass estimate" << endl;
  report << "Single age        # Fishing mortality type such as \"Single age\" or \"exploitation rate\"" << endl;
  report << "Age model         # Fishing mortality source such as \"Model\" or \"(total catch (t))/(survey biomass (t))\"" << endl;
  report << "Age of maximum F  # Fishing mortality range such as \"Age of maximum F\"" << endl; 
  report << "#FISHERYDESC -list of fisheries (ALL TWL LGL POT FIX FOR DOM TWLJAN LGLMAY POTAUG ...)" << endl; 
  report << "ALL" << endl; 

  report <<"#FISHERYYEAR - list years used in the model " << endl;
    for (i=styr;  i<=endyr; i++)
       report << i << "	";
       report<<endl;  

  report<<"#AGE - list of ages used in the model"<<endl;
    for (i=3; i<=15;i++)
       report << i << "	";
       report<<endl;    

  report <<"#RECRUITMENT - Number of recruits by year " << endl;
    for (i=styr;  i<=endyr;  i++)
	   report  << 2*natage(1,i,1) << "	";
	   report<<endl;     
	
  report <<"#SPAWNBIOMASS - Spawning biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << fspbio(i) << "	";
       report<<endl;  

  report <<"#TOTALBIOMASS - Total biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << natage(1,i)*wt(1) + natage(2,i)*wt(2) << "	";
       report<<endl;
 	
  report <<"#TOTFSHRYMORT - Fishing mortality rate by year " << endl;
	for (i=styr;  i<=endyr;  i++)
	   report  << (F(1,i,21)+ F(2,i,21))/2<< "	";
	   report<<endl;
	  
  report <<"#TOTALCATCH - Total catch by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << catch_bio(i) << "	";
       report<<endl;
   
 report <<"#MATURITY - Maturity ratio by age (females only)" << endl;  
       for (i=1;  i<=21;  i++) 
       report  << maturity(i) <<"	";
       report<< endl; 
 
 report <<"#SPAWNWT - Average spawning weight (in kg) by age"<< endl; 
       report <<"0.019	0.041	0.113	0.224	0.376	0.566	0.784	1.028	1.292	1.569	1.855	2.142	2.417	2.667	2.881	3.057	3.198	3.308	3.393"<<endl;                              
       report<<endl;
                     
 report <<"#NATMORT - Natural mortality rate for females then males"<< endl; 
 for (i=1;  i<=21;  i++) 
 report  << 0.2 <<"	";
 report<< endl;   
 for (i=1;  i<=21;  i++) 
 report  << 0.35 <<"	";
 report<< endl;

 report << "#N_AT_AGE - Estimated numbers of female (first) then male (second) fish at age " << endl;
   for (i=styr; i<=endyr;i++)
     report <<natage(1,i)<< "	";
     report<<endl;

   for (i=styr; i<=endyr;i++)
     report <<natage(2,i)<< "	";
     report<<endl;
 
 report <<"#FSHRY_WT_KG - Fishery weight at age (in kg) females (first) males (second), only one fishery"<< endl;   
    report <<wt(1)*1000  << "	";
    report<<endl; //1 is females        
  
    report <<wt(2)*1000  << "	";
    report<<endl; //2 is males


  report << "#SELECTIVITY - Estimated fishery selectivity for females (first) males (second) at age " << endl;
    for (j=1; j<=nages;j++)
      report <<" " <<sel(1,j)<< "	";
      report<<endl;

    for (j=1; j<=nages;j++)
      report <<" "  <<sel(2,j)<< "	";
      report<<endl;
 report << "#SURVEYDESC"<<endl;
 report<<"EBS_trawl_survey BS_slope_trawl_survey AI_trawl_survey"<<endl;
 
 report<<"SURVEYMULT"<<endl;
 report<<"1 1 1"<<endl;
 report << "#EBS_trawl_survey - Bering Sea shelf survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nobs_srv1;i++)
      report << yrs_srv1(i) << "	";
      report<<endl;
    for (i=1; i<=nobs_srv1;i++) 
      report<< obs_srv1(i)<< "	";
      report<< endl;

 report << "#BS_slope_trawl_survey - Bering Sea slope survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nobs_srv2;i++)
      report << yrs_srv2(i) << "	";
      report<<endl;
    for (i=1; i<=nobs_srv2;i++)
      report << obs_srv2(i) << "	";
      report<<endl;

 report << "#AI_trawl_survey - Aleutian Islands survey biomass (Year, Obs_biomass, Pred_biomass) "  << endl;
    for (i=1; i<=nobs_srv3;i++)
      report << yrs_srv3(i) << "	"; 
      report<<endl;
    for (i=1; i<=nobs_srv3;i++)
      report << obs_srv3(i) << "	";
      report<<endl;	

 report<<"#STOCKNOTES"<<endl;
 report<<"SAFE report indicates that this stock was not subjected to overfishing in 2012 and is neither overfished nor approaching a condition of being overfished in 2013."<<endl;
RUNTIME_SECTION
  maximum_function_evaluations 4000
  convergence_criteria 1e-3 1e-4 1e-7

TOP_OF_MAIN_SECTION
  arrmblsize = 20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(300);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);

