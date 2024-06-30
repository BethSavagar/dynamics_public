# summary stats function for demographics model

summarydf <- function(
    w,
    f_list,
    m_list,
    output,
    summary_df,
    Kid,
    Sub,
    Adu_F
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
  
  if(output %in% c("summary", "summary_all", "dynamics", "pop_tracker", "mort_only", "off_only")){
    # total population size
    fpop <- fIm+fS+fE+fI+fR
    mpop <- mIm+mS+mE+mI+mR
    pop_tot <- fpop+mpop
    sum_pop <- sum(fpop)+sum(mpop)
    
    # proportion immune
    fimmune <- fIm+fR
    mimmune <- mIm+mR
    sum_immune <- sum(fimmune)+sum(mimmune)
    prop_immune <- sum_immune / sum_pop
    
    # proportion immune excluding offspring
    sumR <- sum(fR)+sum(mR)
    pR_noIm <- sumR/sum_pop
    
    # proportion infectious
    sum_inf <- sum(fI)+sum(mI)
    prop_inf <- sum_inf/sum_pop
    
    # proportion in each age-group
    Kid_tot <- pop_tot[Kid] ; sum_Kid <- sum(Kid_tot) ; pKid <- sum_Kid/sum_pop
    Sub_tot <- pop_tot[Sub] ; sum_Sub <- sum(Sub_tot) ; pSub <- sum_Sub/sum_pop
    Adu_tot <- pop_tot[Adu_F] ; sum_Adu <- sum(Adu_tot) ; pAdu <- sum_Adu/sum_pop
    
    # sex-age group
    fKid_tot <- fpop[Kid] ; sum_fKid <- sum(fKid_tot) ; pfKid <- sum_fKid/sum_pop
    mKid_tot <- mpop[Kid] ; sum_mKid <- sum(mKid_tot) ; pmKid <- sum_mKid/sum_pop
    fSub_tot <- fpop[Sub] ; sum_fSub <- sum(fSub_tot) ; pfSub <- sum_fSub/sum_pop
    mSub_tot <- mpop[Sub] ; sum_mSub <- sum(mSub_tot) ; pmSub <- sum_mSub/sum_pop
    fAdu_tot <- fpop[Adu_F] ; sum_fAdu <- sum(fAdu_tot) ; pfAdu <- sum_fAdu/sum_pop
    mAdu_tot <- mpop[Adu_F] ; sum_mAdu <- sum(mAdu_tot) ; pmAdu <- sum_mAdu/sum_pop
    
    # proportion female
    pF <- sum(fpop)/sum(pop_tot)
    
    # number of births
    nBirths <- pop_tot[1]
  }
  
  # add row to summary stats dataframe
  if (output == "summary") {
    summary_df[w, ] <-c(w,sum_pop,prop_immune,pR_noIm,prop_inf,pF,pKid,pSub,pAdu)
  } else if (output %in% c("summary_all")) {
    summary_df[w, ] <- c(w,sum_pop,prop_immune,pR_noIm,prop_inf,pF,pfKid,pfSub,pfAdu,pmKid,pmSub,pmAdu)
  } else if (output %in% c("dynamics")) {
    # nBirths
    summary_df[w, ] <- c(w,sum_pop,prop_immune,pR_noIm,prop_inf,nBirths,pF,pfKid,pfSub,pfAdu,pmKid,pmSub,pmAdu)
  } else if (output %in% c("pop_tracker", "mort_only", "off_only")) {
    pop_w <- cbind(w,fpop,mpop,pop_tot,sum_pop)
    summary_df <- rbind(summary_df,pop_w)
    count_df <- c()
  }
  # Output of Model: 
  if(output == "summary"){
    return(summary_df)
  }else if(output == "summary_all"){
    return(summary_df)
  }else if(output == "dynamics"){
    return(summary_df) 
  }else if(output %in% c("pop_tracker", "mort_only", "off_only")){
    return(summary_df)
  }
}
  
  


