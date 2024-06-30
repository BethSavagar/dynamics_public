dynmod_wrapper <- function(
    imm_decay_corrected, # immunity decay rate
    demographics_df_full, # state vars
    # f_list, # initial state of female population
    # m_list, # initial state of male population
    TimeStop_dynamics, # 1 year, weekly timestep for demographic component
    # TimeStop_transmission, # 1 day, hourly timestep for transission component
    output, # model output: tracker or summary stats
    summary_df, # empty output dataframe
    Vstart, # time of vaccination
    Vprog,
    birthpeak_w,
    seasonal
){
  
  #########################
  ## 1. PARAMETERISATION ##
  #########################
  
  ###
  # -------------------
  ## State Variables ##
  # -------------------
  
  ## Flock Size ##
  N_tot <- demographics_df %>% filter(parameter=="N_tot") %>% pull(default) # number of animals in population (flock or village)?
  
  ## Flock Structure ##
  # set proportion of male/female: kid (<6m), sub (<12m), adu (12m+)
  
  kid_f_prop <- demographics_df %>% filter(parameter=="kid_f_prop") %>% pull(default) 
  sub_f_prop <- demographics_df %>% filter(parameter=="sub_f_prop") %>% pull(default) 
  adu_f_prop <- demographics_df %>% filter(parameter=="adu_f_prop") %>% pull(default) 
  
  kid_m_prop <- demographics_df %>% filter(parameter=="kid_m_prop") %>% pull(default) 
  sub_m_prop <- demographics_df %>% filter(parameter=="sub_m_prop") %>% pull(default) 
  adu_m_prop <- demographics_df %>% filter(parameter=="adu_m_prop") %>% pull(default) 
  
  age_p_f <- c("kid_f_p"=kid_f_prop,"sub_f_p"=sub_f_prop,"adu_f_p"=adu_f_prop) # sex-age group proportions (INITIAL)
  age_p_m <- c("kid_m_p"=kid_m_prop,"sub_m_p"=sub_m_prop,"adu_m_p"=adu_m_prop)
  
  ## Births in peak month ##
  p_births <- demographics_df %>% filter(parameter=="pBirths") %>% pull(default)  # proportion initially in R state
  
  ## Age Group Limits ##
  
  wk2mnth <- 4.345 # weeks per month
  
  # kids
  kid_max <- demographics_df %>% filter(parameter=="kid_max") %>% pull(default) 
  kid_max_wks <- round(kid_max*wk2mnth)
  # sub-adults
  sub_max <- demographics_df %>% filter(parameter=="sub_max") %>% pull(default) 
  sub_max_wks <- round(sub_max*wk2mnth)
  # adults
  max_age_F <- demographics_df %>% filter(parameter == "adu_f_max_yrs") %>% pull(default) *52
  max_age_M <- demographics_df %>% filter(parameter == "adu_m_max_yrs") %>% pull(default) *52
  # weeks in age class
  Kid <- 1:kid_max_wks # Kid: 1-6m (1-26w)
  Sub <- (kid_max_wks+1):sub_max_wks # Sub: 6-12m (26-52w)
  Adu_F <- (sub_max_wks+1):max_age_F # Adult F: 12m+
  Adu_M <- (sub_max_wks+1):max_age_M # Adult M: 12m+
  
  ## Maternal Immunity ## (Bodjo et al) 
  Imm_b <- imm_decay_corrected %>% filter(wk ==0) %>% pull(imm_corrected) # Imm_b = # proportion of young born to immune mothers that gain maternal antibodies
  ## Immune Fraction ##
  pR <- demographics_df %>% filter(parameter=="pR") %>% pull(default)  # proportion initially in R state
  # -------------------
  ## Demographics ##
  # -------------------
  
  ## ppr_mortality not incorporated
  ppr_mort_1 <- 0
  ppr_mort_2 <- 0
  
  ## Pull Demographics ##
  off_1 <- demographics_df %>% filter(parameter == "NET_offtake_y") %>% pull(default)
  off_F <- demographics_df %>% filter(parameter == "NET_offtake_f") %>% pull(default)
  off_M <- demographics_df  %>% filter(parameter == "NET_offtake_m") %>% pull(default)
  off_M2 <- demographics_df  %>% filter(parameter == "NET_offtake_m2") %>% pull(default)
  mort_1 <- demographics_df  %>% filter(parameter == "mortality_y") %>% pull(default)
  mort_2 <- demographics_df  %>% filter(parameter == "mortality_a") %>% pull(default)
  # birth_r <- demographics_df  %>% filter(parameter == "birth_rate") %>% pull(default)
  nBirths <- demographics_df  %>% filter(parameter == "nBirths_yr") %>% pull(default)
  npeak_wks <- demographics_df %>% filter(parameter == "peak_wks") %>% pull(default)
  min_age_offtake <- demographics_df  %>% filter(parameter == "min_age_offtake") %>% pull(default)
  min_age_repro <- demographics_df  %>% filter(parameter == "min_age_repro") %>% pull(default)
  
  
  ## Seasonal Births ##
  nbirthH <- (nBirths*p_births)/npeak_wks # high birth rate for 1 month (4wks)
  nbirthL <- (nBirths*(1-p_births))/(52-npeak_wks) # low birth rate for 11months
  nbirthW <-  nBirths/52  # weekly births
  
  ## Convert to weekly rates ##
  # (if yearly):

  off_1 <- 1-((1-off_1)^(1/52))
  off_F <- 1-((1-off_F)^(1/52))
  off_M <- 1-((1-off_M)^(1/52))
  off_M2 <- 1-((1-off_M2)^(1/52))
  mort_1 <- 1-((1-mort_1)^(1/52))
  mort_2 <- 1-((1-mort_2)^(1/52))
  
  # Convert ages to wks
  min_offtake_wks <- round(min_age_offtake*wk2mnth)
  max_M2off_wks <- round(24*wk2mnth) # most male offtake before 12-18months 
  min_repro_wks <- round(min_age_repro*wk2mnth)
  
  
  ## Construct Demographics DF ##
  demographic_pars <- data.frame(
    age_weeks = 1:max_age_F) %>% # nrow = age in weeks
    mutate(age_cat = ifelse(age_weeks <= kid_max_wks, "Kid",
                            ifelse(age_weeks <= sub_max_wks, "Sub", "Adu"))) %>%
    # fill in maternal immunity
    left_join(imm_decay_corrected %>% select(wk, "imm" = imm_corrected), c("age_weeks" = "wk")) %>%
    
    mutate(imm = ifelse(is.na(imm),0,imm)) %>%
    ## Join Demographics Data >>>
    mutate(net_off_F = ifelse(age_weeks<min_offtake_wks, off_1, off_F),
           net_off_M = ifelse(age_weeks<min_offtake_wks, off_1,
                              ifelse(age_weeks<max_M2off_wks, off_M2, off_M)), # add additional offtake rate for males <2y 
           # mort_F = ifelse(age_weeks<=kid_max, mort_1, mort_2),
           # mort_F = ifelse(age_weeks == max_age_F, 1, mort_F),
           mort_F = ifelse(age_weeks<=kid_max_wks, mort_1, 
                           ifelse(age_weeks == max_age_F, 1, mort_2)),
           mort_M = ifelse(age_weeks<=kid_max_wks, mort_1, 
                           ifelse(age_weeks == max_age_M, 1, mort_2)),
           ppr_mort = ifelse(age_weeks<=kid_max_wks, ppr_mort_1,ppr_mort_2), # ppr_mortality
           # birth = ifelse(age_weeks<min_repro_wks,0,birth_r), # birth rate (only Adu_F)
           
           age_p_F = ifelse(age_cat=="Kid", kid_f_prop,
                            ifelse(age_cat=="Sub", sub_f_prop, adu_f_prop)),
           age_p_M = ifelse(age_cat=="Kid", kid_m_prop,
                            ifelse(age_cat=="Sub", sub_m_prop, 
                                   ifelse(age_cat =="Adu" & age_weeks<=max_age_M, adu_m_prop,0))),
           n_weeks_F = ifelse(age_cat=="Kid", kid_max_wks,
                              ifelse(age_cat=="Sub", sub_max_wks - kid_max_wks, 
                                     ifelse(age_cat =="Adu" & age_weeks<=max_age_F, max_age_F - sub_max_wks,0))),
           n_weeks_M = ifelse(age_cat=="Kid", kid_max_wks,
                              ifelse(age_cat=="Sub", sub_max_wks - kid_max_wks, 
                                     ifelse(age_cat =="Adu" & age_weeks<=max_age_M, max_age_M - sub_max_wks,0)))
    ) %>%
    
    # divide initial population between age groups according to proportions in age_p_F
    mutate(pop_init_F = (age_p_F*N_tot)/n_weeks_F,
           pop_init_M = (age_p_M*N_tot)/n_weeks_M, 
           pop_init_M = ifelse(is.na(pop_init_M), 0, pop_init_M)) %>%
    # select relevant variables   
    select(age_weeks,
           age_cat,
           imm,
           # birth,
           mort_F,
           mort_M,
           ppr_mort,
           net_off_F,
           net_off_M,
           pop_init_F,
           pop_init_M
    )
  
  nBirths <- data.frame("nbirthW" = nbirthW, "nbirthH" = nbirthH, "nbirthL" = nbirthL)
  
  
  #########################
  ## 2. SIMULATION SETUP ##
  #########################
  
  #--------------------------
  ## Initialise Population ##
  #--------------------------
  
  fIm_init <- rep(0,max_age_F); mIm_init <- rep(0,max_age_F)
  # fIm_init <- demographic_pars %>% pull(pop_init_F)*pIm
  # mIm_init <- demographic_pars %>% pull(pop_init_M)*pIm
  fS_init <- demographic_pars %>% pull(pop_init_F) *(1-pR) # all S except already recovered
  mS_init <- demographic_pars %>% pull(pop_init_M) *(1-pR) # all S except already recovered
  fE_init <- rep(0,max_age_F); mE_init <- rep(0,max_age_F)
  # # fE_init <- demographic_pars %>% pull(pop_init_F)*pE
  # # mE_init <- demographic_pars %>% pull(pop_init_M)*pE
  fI_init <- rep(0,max_age_F); mI_init <- rep(0,max_age_F)
  # # fI_init <- demographic_pars %>% pull(pop_init_F)*pI 
  # # mI_init <- demographic_pars %>% pull(pop_init_M)*pI
  fR_init <- demographic_pars %>% pull(pop_init_F)*pR
  mR_init <- demographic_pars %>% pull(pop_init_M)*pR
  
  f_list <- list("fIm"=fIm_init,
                 "fS"=fS_init,
                 "fE"=fE_init,
                 "fI"=fI_init,
                 "fR"=fR_init)
  m_list <- list("mIm"=mIm_init,
                 "mS"=mS_init,
                 "mE"=mE_init,
                 "mI"=mI_init,
                 "mR"=mR_init)
  
  ###################
  ## 3. RUN DYNMOD ##
  ###################
  
  ## Update output "summary_df" ##
  # f_list/m_list is initial state for Females and Males, output defines what stats are provided in output, sumamry_df is a df to keep output in.
  summary_df <- summarydf(w = 1, f_list, m_list, output, summary_df, Kid, Sub, Adu_F)
  
  output_df <- dynmod_sim(
    f_list, # initial state of female population
    m_list, # initial state of male population
    TimeStop_dynamics, # 1 year, weekly timestep for demographic component
    # TimeStop_transmission, # 1 day, hourly timestep for transission component
    output, # model output: tracker or summary stats
    demographic_pars, # df containing demographic rates (birth rate, mortality etc) for every week-long age group
    summary_df, #
    Imm_b,
    Kid,
    Sub,
    Adu_F,
    Vstart,
    Vprog,
    birthpeak_w,
    seasonal,
    nBirths
  )
  
  # only works for summary all needs updating for summary (only)
  if(output %in% c("summary", "summary_all")){
    # res <- output_df %>%
    #   filter(w==TimeStop_dynamics) %>%
    #   mutate(pop_growth = output_df[nrow(output_df),"sum_pop"] / output_df[1,"sum_pop"],
    #          tenyr_growth = (output_df[t1, "sum_pop"]) / output_df[t2, "sum_pop"]
    #          )
    res <- GSA_output(output_df,Vstart)
  }
  
  if(output %in% c("dynamics", "pop_tracker")){
    res <- output_df
  }
  
  if(output %in% c("mort_only")){
    
    ######
    ## F Mortality
    ######
    
    # Female 0-12m Mortality
    outMort <- output_df %>% mutate(age = rep(1:(nrow(output_df)/TimeStop_dynamics), TimeStop_dynamics))
    F0init <- outMort %>% filter(w == 1, age == 1) %>% pull("fpop")
    F12end <- outMort %>% filter(w == 52, age == 52) %>% pull("fpop")
    annMortF012 <- 1 - F12end / F0init # annual mortality
    mortF012 <- annMortF012
    
    # Female 12m+ Mortality
    F12init <- outMort %>% filter(w == 53, age == 53) %>% pull("fpop")
    Fend <- outMort %>% filter(w == 104, age == 104) %>% pull("fpop")
    annMortF12end <- 1 - Fend / F12init # annual mortality
    mortF12end <- annMortF12end
    
    # Female 0-6m
    F6end <- outMort %>% filter(w == 26, age == 26) %>% pull("fpop")
    monMortF06 <- 1 - (F6end/F0init)
    annMortF06 <- 1- (1-monMortF06)^2
    mortF06 <- annMortF06
    
    # Female 6-12m
    F6init <- outMort %>% filter(w == 26, age == 26) %>% pull("fpop")
    monMortF612 <- 1 - (F12end/F6init)
    annMortF612 <- 1- (1-monMortF612)^2
    mortF612 <- annMortF612
    
    ######
    ## M Mortality
    ######
    
    # Male 0-12m Mortality
    M0init <- outMort %>% filter(w == 1, age == 1) %>% pull("mpop")
    M12end <- outMort %>% filter(w == 52, age == 52) %>% pull("mpop")
    annMortM012 <- 1 - M12end / M0init # annual mortality
    mortM012 <- annMortM012
    
    # Male 12m+ Mortality
    M12init <- outMort %>% filter(w == 53, age == 53) %>% pull("mpop")
    Mend <- outMort %>% filter(w == 104, age == 104) %>% pull("mpop")
    annMortM12end <- 1 - Mend / M12init # annual mortality
    mortM12end <- annMortM12end
    
    
    # Male 0-6m
    M6end <- outMort %>% filter(w == 26, age == 26) %>% pull("mpop")
    monMortM06 <- 1 - (M6end/M0init)
    annMortM06 <- 1- (1-monMortM06)^2
    mortM06 <- annMortM06
    
    # Male 6-12m
    M6init <- outMort %>% filter(w == 26, age == 26) %>% pull("mpop")
    monMortM612 <- 1 - (M12end/M6init)
    annMortM612 <- 1- (1-monMortM612)^2
    mortM612 <- annMortM612
    
    # Male 6m+ 
    M18end <- outMort %>% filter(w == 78, age == 78) %>% pull("mpop")
    annMortM6end <- 1 - M18end / M6init # annual 
    mortM6end <- annMortM6end
    
    # All Mortality
    mortinit <- outMort %>% filter(w == 1, age == 1) %>% pull("sum_pop")
    mortend <- outMort %>% filter(w == 52, age == 52) %>% pull("sum_pop")
    mortAll <- 1-(mortend / mortinit)
    
    # Female Offtake Annual
    Finit <- outMort %>% filter(w==1) %>% summarise(sum(fpop)) %>% pull()
    Fend <- outMort %>% filter(w==52) %>% summarise(sum(fpop)) %>% pull()
    mortF <- 1- (Fend / Finit)
    
    # Male Offtake Annual
    Minit <- outMort %>% filter(w==1) %>% summarise(sum(mpop)) %>% pull()
    Mend <- outMort %>% filter(w==52) %>% summarise(sum(mpop)) %>% pull()
    mortM <- 1- (Mend / Minit)
    
    
    res <- c("mortF012"=mortF012, "mortF12end"=mortF12end, "mortM012"=mortM012, "mortM12end"=mortM12end,
             "mortF06"=mortF06, "mortF612"=mortF612, "mortM06"=mortM06, "mortM612"=mortM612, "mortM6end" = mortM6end,
             "mortF" = mortF, "mortM" = mortM, "mortAll" = mortAll)
    
  }
  
  if(output %in% c("off_only")){
    
    ######
    ## F Offtake
    ######
    
    # Female 0-12m Offtake
    outOff <- output_df %>% mutate(age = rep(1:(nrow(output_df)/TimeStop_dynamics), TimeStop_dynamics))
    F0init <- outOff %>% filter(w == 1, age == 1) %>% pull("fpop")
    F12end <- outOff %>% filter(w == 52, age == 52) %>% pull("fpop")
    annOffF012 <- 1 - F12end / F0init # annual offtake
    offF012 <- annOffF012
    
    # Female 12m+ Offtake
    F12init <- outOff %>% filter(w == 53, age == 53) %>% pull("fpop")
    Fend <- outOff %>% filter(w == 104, age == 104) %>% pull("fpop")
    annOffF12end <- 1 - Fend / F12init # annual offtake
    offF12end <- annOffF12end
    
    # Female 0-6m
    F6end <- outOff %>% filter(w == 26, age == 26) %>% pull("fpop")
    monOffF06 <- 1 - (F6end/F0init)
    annOffF06 <- 1- (1-monOffF06)^2
    offF06 <- annOffF06
    
    # Female 6-12m
    F6init <- outOff %>% filter(w == 26, age == 26) %>% pull("fpop")
    monOffF612 <- 1 - (F12end/F6init)
    annOffF612 <- 1- (1-monOffF612)^2
    offF612 <- annOffF612
    
    ######
    ## M offtake
    ######
    
    # Male 0-12m offtake
    M0init <- outOff %>% filter(w == 1, age == 1) %>% pull("mpop")
    M12end <- outOff %>% filter(w == 52, age == 52) %>% pull("mpop")
    annOffM012 <- 1 - M12end / M0init # annual offtake
    offM012 <- annOffM012
    
    # Male 12m+ offtake
    M12init <- outOff %>% filter(w == 53, age == 53) %>% pull("mpop")
    Mend <- outOff %>% filter(w == 104, age == 104) %>% pull("mpop")
    annOffM12end <- 1 - Mend / M12init # annual offtake
    offM12end <- annOffM12end
    
    # Male 0-6m
    M6end <- outOff %>% filter(w == 26, age == 26) %>% pull("mpop")
    monOffM06 <- 1 - (M6end/M0init)
    annOffM06 <- 1- (1-monOffM06)^2
    offM06 <- annOffM06
    
    # Male 6-12m
    M6init <- outOff %>% filter(w == 26, age == 26) %>% pull("mpop")
    monOffM612 <- 1 - (M12end/M6init)
    annOffM612 <- 1- (1-monOffM612)^2
    offM612 <- annOffM612
    
    # Male 6m+ offtake
    M18end <- outOff %>% filter(w == 78, age == 78) %>% pull("mpop")
    annOffM6end <- 1 - M18end / M6init # annual offtake
    offM6end <- annOffM6end
    
    # All Offtake
    offinit <- outOff %>% filter(w == 1, age == 1) %>% pull("sum_pop")
    offend <- outOff %>% filter(w == 52, age == 52) %>% pull("sum_pop")
    offAll <- 1-(offend / offinit)
    
    # Female Offtake Annual
    Finit <- outOff %>% filter(w==1) %>% summarise(sum(fpop)) %>% pull()
    Fend <- outOff %>% filter(w==52) %>% summarise(sum(fpop)) %>% pull()
    offF <- 1- (Fend / Finit)
    
    # Male Offtake Annual
    Minit <- outOff %>% filter(w==1) %>% summarise(sum(mpop)) %>% pull()
    Mend <- outOff %>% filter(w==52) %>% summarise(sum(mpop)) %>% pull()
    offM <- 1- (Mend / Minit)
    
    res <- c("offF012"=offF012, "offF12end"=offF12end, "offM012"=offM012, "offM12end"=offM12end,
             "offF06"=offF06, "offF612"=offF612, "offM06"=offM06, "offM612"=offM612,"offM6end"=offM6end,
             "offF" = offF, "offM" = offM, "offAll"=offAll)
    
  }
  
  return(res)
  
}

