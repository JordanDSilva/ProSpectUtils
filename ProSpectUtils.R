#ProSpect utils

## Dust mass stuff from D'Silva+26
.New_DustMass_vDTH = function(out, qPAH_VSG = 0.035){
  
  new_dust = out$SEDout$dustlum["birth"]/Dale_M2L_variableDTH_func(out$parm["alpha_SF_birth"], qPAH_VSG = qPAH_VSG) + 
    out$SEDout$dustlum["screen"]/Dale_M2L_variableDTH_func(out$parm["alpha_SF_screen"], qPAH_VSG = qPAH_VSG)
  
  return(as.numeric(new_dust))
}
.RR14_BPL = function(Z, doDTG = FALSE){
  
  ## Remy Ruyer+14 using metallicity dependent XCO
  ## x = 12 + log(O/H)
  ## xSol = 12 + log(O/H)Sol
  
  ## Z/Zsol = (O/H)/(O/HSol)
  ## log(Z) - log(Zsol) = log(O/H) - log(O/HSol)
  ## log(Z)+12 - log(Zsol)-12 = log(O/H)+12 - log(O/Hsol)-12
  ## log(Z) - log(Zsol) = x - xSol
  ## log(Z/0.014) + xSol = log(O/H) + 12
  
  a = 2.21
  alphaH = 1.00
  b = 0.96
  alphaL = 3.10
  xt = 8.10

  xSol = 8.69
  ZOH = log10(Z / 0.014) + xSol
  
  GTD = ifelse(
    ZOH > xt,
    a + alphaH*(xSol - ZOH),
    b + alphaL*(xSol - ZOH)
  )
  DTG = (10^GTD)^-1
  
  ##Mgas = muGal * Mhydrogen
  ## DTG = Mdust / Mgas
  ## DTH = Mdust / Mhydrogen = Mdust / (Mgas / muGal) = DTG * muGal = DTG * (1 / (1 - Ysol - Zgal))
  muGal = 1 / (1 - 0.270 - Z)
  DTH = DTG*muGal
  if(doDTG){
    return(DTG)
  }else{
    return(DTH)
  }
}

# shivaei24 = fread("~/Documents/DustMassDensity/data/shivaeqPAHZ.csv")
shivaei24 = data.frame(
  "OH" = c(7.8741935483871, 7.99354838709677, 8.09516129032258, 8.24838709677419, 8.39354838709678, 8.57903225806452, 8.83064516129032),
  "q" = c(1.04568527918782, 1.01522842639594, 1.01522842639594, 2.20304568527919, 3.39086294416244, 3.39086294416244, 3.3756345177665) 
)
shivaei24_fit = approxfun(
  shivaei24$OH,
  shivaei24$q,
  rule = 2
)
.shivaei24_qPAHZ = function(Z){
  xSol = 8.69
  ZOH = log10(Z / 0.014) + xSol
  
  ff = c(shivaei24_fit(ZOH)) * 0.01
  ff[ZOH >= 8.4] =  0.035
  ff[ZOH < 8.1] =  0.01
  return(ff)
}

## Now get some statistics and astrophysics 
.UV_tophat = function(wave){
  throughput = rep(0L, length(wave))
  throughput[wave >= (1500 - 50) & wave <= (1500 + 50)] = 1L
  return(throughput)
}
.UVLum1500 = function(out){
  out$Data$fit = "check"
  
  UV1500 = photom_lum(
    wave = out$SEDout$FinalLum$wave,
    lum = out$SEDout$FinalLum$lum,
    z = 0.0, 
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref,
    filters = list(.UV_tophat)
  )
  
  return(UV1500)
}
.UVLumAGN1500 = function(out){
  out$Data$fit = "check"
  
  UV1500 = photom_lum(
    wave = out$SEDout$AGN$wave,
    lum = out$SEDout$AGN$lum,
    z = 0.0, 
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref,
    filters = list(.UV_tophat)
  )
  
  return(UV1500)
}
.stellar_mass = function(out){
  
  out$Data$fit = "check"
  
  parm = out$parm
  parm[out$Data$logged] = 10^parm[out$Data$logged]
  
  SMstar = ParmOff(
    .func = SMstarfunc, #the function we want to run
    .args = c(parm, out$Data$arglist), #the superset of potential matching parameters
    massfunc = out$Data$arglist$massfunc,
    speclib = BC03hr,
    Z = Zfunc_massmap_lin
  )
  return(as.numeric(SMstar['TotSMstar']))
}
.SFR10 = function(out){
  
  out$Data$fit = "check"
  
  parm = out$parm
  parm[out$Data$logged] = 10^parm[out$Data$logged]
  
  mass_func_args_idx = out$Data$parm.names %in% names(formals(out$Data$arglist$massfunc))
  massfunc_args = c(parm[mass_func_args_idx], "magemax" = out$Data$arglist$magemax)
  
  sfr10 = do.call('integrate', c(list(
    f = out$Data$arglist$massfunc, lower = 0, upper = 1e7
  ), massfunc_args))$value / 1e7
  
  return(sfr10)
}
.fitSEDObs = function(out, doAGN = TRUE, doDust = TRUE){
  
  waveout = out$Data$waveout
  
  fflux_StarsUnAtten = Lum2Flux(
    out$SEDout$StarsUnAtten$wave,
    out$SEDout$StarsUnAtten$lum,
    z = out$SEDout$z,
    H0 = out$SEDout$cosmo$H0,
    OmegaM = out$SEDout$cosmo$OmegaM,
    OmegaL = out$SEDout$cosmo$OmegaL,
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref
  )
  fflux_StarsUnAttenJy = CGS2Jansky(convert_wave2freq(fflux_StarsUnAtten$flux, fflux_StarsUnAtten$wave))
  
  fflux_StarsAtten = Lum2Flux(
    out$SEDout$StarsAtten$wave,
    out$SEDout$StarsAtten$lum,
    z = out$SEDout$z,
    H0 = out$SEDout$cosmo$H0,
    OmegaM = out$SEDout$cosmo$OmegaM,
    OmegaL = out$SEDout$cosmo$OmegaL,
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref
  )
  fflux_StarsAttenJy = CGS2Jansky(convert_wave2freq(fflux_StarnAtten$flux, fflux_StarnAtten$wave))
  
  if(doAGN){
    fflux_AGN = Lum2Flux(
      out$SEDout$AGN$wave,
      out$SEDout$AGN$lum,
      z = out$SEDout$z,
      H0 = out$SEDout$cosmo$H0,
      OmegaM = out$SEDout$cosmo$OmegaM,
      OmegaL = out$SEDout$cosmo$OmegaL,
      LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
      ref = out$Data$arglist$ref
    )
    fflux_AGNJy = CGS2Jansky(convert_wave2freq(fflux_AGN$flux, fflux_AGN$wave))
  }else{
    fflux_AGN = rep(0, length(waveout))
  }
  if(doDust){
    fflux_Dust = Lum2Flux(
      out$SEDout$DustEmit$wave,
      out$SEDout$DustEmit$lum,
      z = out$SEDout$z,
      H0 = out$SEDout$cosmo$H0,
      OmegaM = out$SEDout$cosmo$OmegaM,
      OmegaL = out$SEDout$cosmo$OmegaL,
      LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
      ref = out$Data$arglist$ref
    )
    fflux_DustJy = CGS2Jansky(convert_wave2freq(fflux_Dust$flux, fflux_Dust$wave))
  }else{
    fflux_Dust = rep(0, length(waveout))
  }
  
  ffluxStarsUnAttenJy_rebin = specReBin(
    wave = fflux_StarsUnAtten$wave,
    flux = fflux_StarsUnAttenJy,
    wavegrid = out$SEDout$FinalFlux$wave
  )
  ffluxStarsAttenJy_rebin = specReBin(
    wave = fflux_StarsAtten$wave,
    flux = fflux_StarsAttenJy,
    wavegrid = out$SEDout$FinalFlux$wave
  )
  ffluxAGNJy_rebin = specReBin(
    wave = fflux_AGN$wave,
    flux = fflux_AGNJy,
    wavegrid = out$SEDout$FinalFlux$wave
  )
  ffluxDustJy_rebin = specReBin(
    wave = fflux_Dust$wave,
    flux = fflux_DustJy,
    wavegrid = out$SEDout$FinalFlux$wave
  )
  return(
    list(
      "wave" = out$SEDout$FinalFlux$wave,
      "flux" = out$SEDout$FinalFlux$flux,
      "starsUnAtten" = ffluxStarsUnAttenJy_rebin$flux,
      "starsAtten" = ffluxStarsAttenJy_rebin$flux,
      "AGN" = ffluxAGNJy_rebin$flux,
      "Dust" = ffluxDustJy_rebin$flux
    )
  )
}
.fitSEDRest = function(out, doAGN = TRUE, doDust = TRUE){
  
  waveout = out$Data$waveout
  
  fflux = Lum2Flux(
    out$SEDout$FinalLum$wave,
    out$SEDout$FinalLum$lum,
    z = out$SEDout$z,
    H0 = out$SEDout$cosmo$H0,
    OmegaM = out$SEDout$cosmo$OmegaM,
    OmegaL = out$SEDout$cosmo$OmegaL,
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref
  )
  ffluxJy = CGS2Jansky(convert_wave2freq(fflux$flux, fflux$wave))
  
  fflux_StarsUnAtten = Lum2Flux(
    out$SEDout$StarsUnAtten$wave,
    out$SEDout$StarsUnAtten$lum,
    z = out$SEDout$z,
    H0 = out$SEDout$cosmo$H0,
    OmegaM = out$SEDout$cosmo$OmegaM,
    OmegaL = out$SEDout$cosmo$OmegaL,
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref
  )
  fflux_StarsUnAttenJy = CGS2Jansky(convert_wave2freq(fflux_StarsUnAtten$flux, fflux_StarsUnAtten$wave))
  
  fflux_StarnAtten = Lum2Flux(
    out$SEDout$StarsAtten$wave,
    out$SEDout$StarsAtten$lum,
    z = out$SEDout$z,
    H0 = out$SEDout$cosmo$H0,
    OmegaM = out$SEDout$cosmo$OmegaM,
    OmegaL = out$SEDout$cosmo$OmegaL,
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref
  )
  fflux_StarsAttenJy = CGS2Jansky(convert_wave2freq(fflux_StarnAtten$flux, fflux_StarnAtten$wave))
  
  if(doAGN){
    ffluxAGN_rebin = specReBin(
      wave = out$SEDout$AGN$wave,
      flux = out$SEDout$AGN$lum,
      wavegrid = out$SEDout$FinalFlux$wave
    )
  }else{
    ffluxAGN_rebin = cbind("wave" = waveout, "flux" = rep(0, length(waveout)))
  }
  if(doDust){
    ffluxDust_rebin = specReBin(
      wave = out$SEDout$DustEmit$wave,
      flux = out$SEDout$DustEmit$lum,
      wavegrid = out$SEDout$FinalFlux$wave
    )
  }else{
    ffluxDust_rebin = cbind("wave" = waveout, "flux" = rep(0, length(waveout)))
  }
  
  ffluxStarsUnAtten_rebin = specReBin(
    wave = out$SEDout$StarsUnAtten$wave,
    flux = out$SEDout$StarsUnAtten$lum,
    wavegrid = out$SEDout$FinalFlux$wave
  )
  ffluxStarsAtten_rebin = specReBin(
    wave = out$SEDout$StarsAtten$wave,
    flux = out$SEDout$StarsAtten$lum,
    wavegrid = out$SEDout$FinalFlux$wave
  )
  
  return(
    list(
      "wave" = out$SEDout$FinalLum$wave,
      "flux" = out$SEDout$FinalLum$lum,
      "starsUnAtten" = ffluxStarsUnAtten_rebin$flux,
      "starsAtten" = ffluxStarsAtten_rebin$flux,
      "AGN" = ffluxAGN_rebin$flux,
      "Dust" = ffluxDust_rebin$flux
    )
  )
}

.calc_astro = function(bestfit, highout, newDust = TRUE, specz = FALSE, cores = 1){
  ## provide a prospect and highlander fit RDS
  ## must have all of the speclib, AGN and Dale templates included in Data 
  
  Data = bestfit$Data
  
  nparm = length(bestfit$parm)
  post = cbind(highout$LD_last$Posterior1, highout$LD_last$Monitor)
  post = data.table(post)
  post = post[order(-LP),]
  
  ## calculating a chisq cut based on the number of parameters to isolate a 1sigma range (from J. Thorne Thesis Codes on GitHub)
  chisq_cut = max(post$LP) - qchisq(0.68,nparm)/2
  toppost = post[post$LP > chisq_cut,]
  toppost[,LP:=NULL]
  
  registerDoParallel(cores = cores)
  astro = foreach(jj = 1:dim(toppost)[1], .combine = rbind) %dopar% {
    if(jj %% 100 == 0){
      message("Sample: ", jj)
    }
    if(newDust){
      ## Use the most up to date dust stuff from D'Silva+26 otherwise use defaults used in the fit and stored in Data
      Data$Dale_M2l_func = function(alpha_SF){
        Dale_M2L_variableDTH_func(alpha_SF = alpha_SF, qPAH_VSG = .shivaei24_qPAHZ(10^toppost[jj, 'Zfinal'])) * (0.0073/.RR14_BPL(10^toppost[jj, 'Zfinal']))
      }
    }
    
    samp = ProSpectSEDlike(parm = unlist(toppost[jj, 1:nparm]), Data = Data)
    
    df = c(
      "StellarMass" = .stellar_mass(out = samp),
      "SFR10" = .SFR10(out = samp),
      "SFR100" = samp$SEDout$Stars$SFRburst,
      "UV1500" = .UVLum1500(out = samp),
      "UV1500AGN" = .UVLumAGN1500(out = samp),
      "DustLumBirth" = as.numeric(samp$SEDout$dustlum['birth']),
      "DustLumScreen" = as.numeric(samp$SEDout$dustlum['screen']),
      "DustLumTotal" = as.numeric(samp$SEDout$dustlum['total']),
      "DustMassBirth" = as.numeric(samp$SEDout$dustmass['birth']),
      "DustMassScreen" = as.numeric(samp$SEDout$dustmass['screen']),
      "DustMassTotal" = as.numeric(samp$SEDout$dustmass['total']),
      "NgasMassBirth" = as.numeric(samp$SEDout$dustmass['birth']) / .RR14_BPL(as.numeric(10^toppost[jj, 'Zfinal']), doDTG = TRUE),
      "NgasMassScreen" = as.numeric(samp$SEDout$dustmass['screen']) / .RR14_BPL(as.numeric(10^toppost[jj, 'Zfinal']), doDTG = TRUE),
      "NgasMassTotal" = as.numeric(samp$SEDout$dustmass['total']) / .RR14_BPL(as.numeric(10^toppost[jj, 'Zfinal']), doDTG = TRUE),
      "LP" = samp$LP
    )
    extra_parm = samp$parm 
    extra_parm[samp$Data$logged] = 10^extra_parm[samp$Data$logged]
    c(extra_parm, df)
  }
  stopImplicitCluster()
  
  astro_quantiles = colQuantiles(
    as.matrix(astro), probs = c(0.5, 0.16, 0.84)
  )
  params = rownames(astro_quantiles)
  
  summary_dt = list(
    astro_quantiles[, "50%"],
    astro_quantiles[, "16%"],
    astro_quantiles[, "84%"],
    astro_quantiles[, "50%"] - astro_quantiles[, "16%"],  # lower uncertainty
    astro_quantiles[, "84%"] - astro_quantiles[, "50%"]   # upper uncertainty
  )
  
  Qnames = c("Q50", "Q16", "Q84", "_lo", "_hi")
  for (x in seq_along(summary_dt)) {
    names(summary_dt[[x]]) = paste0(names(summary_dt[[x]]), Qnames[x])
  }
  
  if(specz){
    summary_dt_nice = as.data.frame(t(c("z" = Data$arglist$z, unlist(summary_dt))))
  }else{
    summary_dt_nice = as.data.frame(t(c(unlist(summary_dt))))
  }
  
  returnDF = list(
    "summary" = summary_dt_nice,
    "Posterior" = astro
  )
  
  return(returnDF)
}
