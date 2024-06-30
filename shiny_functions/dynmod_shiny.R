# 29 Aug - include vacciantion fro GSA
# 21 Nov - Produce mortality and offtake in each age group for applied example validation.

dynmod_sim <- function(
    f_list, # initial state of female population
    m_list, # initial state of male population
    TimeStop_dynamics, # 1 year, weekly timestep for demographic component
    # TimeStop_transmission, # 1 day, hourly timestep for transission component
    output, # model output: tracker or summary stats
    demographic_pars, # demographic transitions by age-sex group
    summary_df, # for storing flock dynamics
    Imm_b, # maternal immunity decay rate
    Kid, # kid age groups (<6m)
    Sub, # sub-adult age groups (6-12m)
    Adu_F, # adult age groups
    Vstart, # start date for vaccination
    Vprog, # full vaccination programme
    birthpeak_w, # weeks of birth peaks (entire programme)
    seasonal, # birth pattern: seasonal/non-seasonal
    nBirths # mean number of births per year
){
  
  fIm <- f_list[["fIm"]]
  fS <- f_list[["fS"]]
  fE <- f_list[["fE"]]
  fI <- f_list[["fI"]]
  fR <- f_list[["fR"]]
  
  mIm <- m_list[["mIm"]]
  mS <- m_list[["mS"]]
  mE <- m_list[["mE"]]
  mI <- m_list[["mI"]]
  mR <- m_list[["mR"]]
  
  ###############################
  ## SET DEMOGRAPHIC PARS ##
  ###############################
  
  # Pull demographic rates as vectors from data-frames (for vectorised equations)
  # NB: create if statements for if there's a difference between males and females?
  immunity <- demographic_pars %>% pull(imm) # no difference m/f
  # Imm_b immune birth rate is already defined (proportion of offspring born with maternal immunity)
  
  net_off_F <- demographic_pars %>% pull(net_off_F)
  net_off_M <- demographic_pars %>% pull(net_off_M)
  mort_F <- demographic_pars %>% pull(mort_F) # no m/f difference
  mort_M <- demographic_pars %>% pull(mort_M)
  # mort_M <- age_pars_M %>% pull(mort)
  ppr_mort <- demographic_pars %>% pull(ppr_mort)
  # ppr_mort_M <- age_pars_M %>% pull(ppr_mort)
  nbirthW <- nBirths %>% pull(nbirthW)
  nbirthL <- nBirths %>% pull(nbirthL)
  nbirthH <- nBirths %>% pull(nbirthH)
  
  ###############################
  ## DEMOGRAPHICS LOOP ##
  ###############################
  
  for(w in 2:TimeStop_dynamics){
    
    # w == 1 is already completed (initial conditions) 
    # week 2 is first week in simulation.
    if(w == 2){
      fIm_prev <- fIm ; mIm_prev <- mIm
      fS_prev <- fS ; mS_prev <- mS
      fE_prev <- fE ; mE_prev <- mE
      fI_prev <- fI ; mI_prev <- mI
      fR_prev <- fR ; mR_prev <- mR
    }
    
    # ---------------------
    ## BIRTHS - fixed ## 

    if(seasonal & w>= TimeStop_dynamics/4){
      if(w %in% birthpeak_w){ # birthpeak_w is a vector which contains the weeks of birth peaks (over the full simulation period i.e. 1, 53...)
        births <- nbirthH # birthsHw: no. births in peak week
      }else{
        births <- nbirthL # birthsHw: no. births in low week.
      }
    }else{
      births <- nbirthW  # birthsW: no. births per week with no seasonality:
    }

    pfR = (sum(fR_prev)/sum(fR_prev+fS_prev))
    pfS = (sum(fS_prev)/sum(fR_prev+fS_prev))

    Im_births <- births* Imm_b* pfR # immune births, pfR is proportion of females in Immune (R) state
    S_births <- births-Im_births # susceptible births = total births per week - immune births per week
    E_births <- 0; I_births <- 0; R_births <- 0 # set EIR births to 0
    
    # ---------------------
    ## SEIR DEMOGRAPHICS ##
    
    # Immune Offspring Demographics (f/m)
    fIm_new <- fIm_prev*(immunity)*(1-net_off_F)*(1-mort_F) # decline in mat immunity, offtake and mortality
    mIm_new <- mIm_prev*(immunity)*(1-net_off_M)*(1-mort_M)
    
    fIm_cur <- c(0.5*Im_births, fIm_new[-length(fIm_new)]) # 0.5*births for F/M, # subset demographics remove last age-group
    mIm_cur <- c(0.5*Im_births, mIm_new[-length(mIm_new)]) # 0.5*births for F/M, # subset demographics remove last age-group
    
    ##
    # Susceptible Demographics
    fS_new <- fS_prev*(1-net_off_F)*(1-mort_F) + # survival (those not offtake, and not died)
      fIm_prev*(1-immunity) # add offspring which have lost immunity
    
    mS_new <- mS_prev*(1-net_off_M)*(1-mort_M) +
      mIm_prev*(1-immunity) # add offspring which have lost immunity
    
    # update fS and mS matrices
    fS_cur <- c(0.5*S_births, fS_new[-length(fS_new)]) # 0.5* births are F, lose last age group
    mS_cur <- c(0.5*S_births, mS_new[-length(mS_new)])# 0.5* births are M, lose last age group
    
    ##
    # Exposed Demographics:
    fE_new <- fE_prev*(1-net_off_F)*(1-mort_F) # F offtake and mortality
    mE_new <- mE_prev*(1-net_off_M)*(1-mort_M) # M offtake and mortality
    
    fE_cur <- c(0.5*E_births,fE_new[-length(fE_new)]) # F add births (E_births = 0)  and remove last age group
    mE_cur <- c(0.5*E_births,mE_new[-length(mE_new)]) # M add births (E_births = 0) and remove last age group
    
    ##
    # Infectious Demographics:
    fI_new <- fI_prev*(1-net_off_F)*(1-mort_F) # F offtake and mortality (**ppr_mortality in disease loop?)
    mI_new <- mI_prev*(1-net_off_M)*(1-mort_M) # M offtake and mortality
    
    fI_cur <- c(0.5*I_births,fI_new[-length(fI_new)]) # F add births (I_births = 0) and remove last age group
    mI_cur <- c(0.5*I_births,mI_new[-length(mI_new)]) # M add births (I_births = 0) and remove last age group
    
    # Recovered Demographics:
    fR_new <- fR_prev*(1-net_off_F)*(1-mort_F) # offtake and mortality (**ppr_mortality in disease loop?)
    mR_new <- mR_prev*(1-net_off_M)*(1-mort_M)
    
    fR_cur <- c(0.5*R_births,fR_new[-length(fR_new)]) # F add births (R_births = 0) and remove last age group
    mR_cur <- c(0.5*R_births,mR_new[-length(mR_new)]) # M add births (R_births = 0) and remove last age group
    
    # ---------------
    ## VACCINATION ##
    # only run if week is in vaccination implementation week in Vprog:
    if(w %in% Vprog$Vweek){
      
      # which round of vaccination (determines animals selected for vac)
      round_i <- Vprog %>% filter(Vweek == w) %>% pull(Vround)
      
      # pull vaccination parameters from Vprog dataframe for current implementation week
      pV <- Vprog %>% filter(Vweek == w) %>% pull(pV)
      Vtype <- Vprog %>% filter(Vweek == w) %>% pull(Vtype)
      Vmin <- Vprog %>% filter(Vweek == w) %>% pull(Vmin)
      Vmax <- Vprog %>% filter(Vweek == w) %>% pull(Vmax)
      
      # Vcov vector contains the pV for each age group
      # Update Vcov to reflect age-groups vaccination in round_i (i.e. Partial or Full)
      
      if(Vtype=="full"){
        # Full: All animals >4m are vaccinated
        Vcov <- rep(pV, length(fS_cur))
        Vcov[0:Vmin] <- 0
      }else if(Vtype=="partial"){
        # Partial: Animals 4-12m are vacciated
        Vcov <- rep(pV, length(fS_cur))
        Vcov[0:Vmin] <- 0
        Vcov[Vmax:length(Vcov)] <- 0
      }
      
      # Update Female population, post-vaccination
      fR_cur <- fR_cur + Vcov*(fS_cur+fE_cur+fI_cur) # new recovered units are all those already immune plus newly vaccinated units from non-immune (SEI) compartments
      fS_cur <- fS_cur*(1-Vcov) 
      fE_cur <- fE_cur*(1-Vcov) 
      fI_cur <- fI_cur*(1-Vcov)
      
      # Update Male population, post-vaccination
      mR_cur <- mR_cur+Vcov*(mS_cur+mE_cur+mI_cur) # new recovered units are all those already immune plus newly vaccinated units from non-immune (SEI) compartments
      mS_cur <- mS_cur*(1-Vcov) 
      mE_cur <- mE_cur*(1-Vcov) 
      mI_cur <- mI_cur*(1-Vcov)
    }
    # ---------------
    
    # -----------------------
    ## TRANSMISSION LOOP ##
    # 
    #
    #
    # ----------------------
    
    # ---------------------
    ## UPDATE POPULATION ##
    
    f_list <- list("fIm"=fIm_cur,
                   "fS"=fS_cur,
                   "fE"=fE_cur,
                   "fI"=fI_cur,
                   "fR"=fR_cur)
    m_list <- list("mIm"=mIm_cur,
                   "mS"=mS_cur,
                   "mE"=mE_cur,
                   "mI"=mI_cur,
                   "mR"=mR_cur)
    
    # -----------------
    ## UPDATE OUTPUT ##
    
    summary_df <- summarydf( w,
                                 f_list,
                                 m_list,
                                 output,
                                 summary_df,
                                 Kid,
                                 Sub,
                                 Adu_F)
    
    # ----------- 
    ## FADEOUT ##
    
    # If there is less than 1 animal in the population at week "w" then set all groups to 0:
    min_pop <- 1
    
    if(summary_df[w,"sum_pop"]<min_pop){
      f_list <- list("fIm"=0,
                     "fS"=0,
                     "fE"=0,
                     "fI"=0,
                     "fR"=0)
      m_list <- list("mIm"=0,
                     "mS"=0,
                     "mE"=0,
                     "mI"=0,
                     "mR"=0)
      
      ## Update Output ##
      summary_df <- summary_demos( w,
                                   f_list,
                                   m_list,
                                   output,
                                   summary_df,
                                   Kid,
                                   Sub,
                                   Adu_F)
    
    } 
    #  
    # # mortality output 
    # if(output == mortSummary){
    #   fIm_mort <- fIm_prev*mort_F
    #   mIm_mort <- mIm_prev*mort_M
    #   fS_mort <- fS_prev*mort_F
    #   mS_mort <- mS_prev*mort_M
    #   fE_mort <- fE_prev*mort_F
    #   mE_mort <- mE_prev*mort_M
    #   fI_mort <- fI_prev*mort_F
    #   mI_mort <- mI_prev*mort_M
    #   fR_mort <- fR_prev*mort_F
    #   mR_mort <- mR_prev*mort_M
    # }
    
    
    ###################
    ## UPDATE STATES ##
    ###################
   
    fIm_prev <- fIm_cur ; mIm_prev <- mIm_cur
    fS_prev <- fS_cur ; mS_prev <- mS_cur
    fE_prev <- fE_cur ; mE_prev <- mE_cur
    fI_prev <- fI_cur ; mI_prev <- mI_cur
    fR_prev <- fR_cur ; mR_prev <- mR_cur
    
    # --------------------------
    ## END DEMOGRAPHICS LOOP ##
  }
  
  ###################
  ## FINAL OUTPUT ##
  ###################
  
  summary_df <- summary_df %>% replace(is.na(.), 0)
  return(summary_df)
}

  
  



