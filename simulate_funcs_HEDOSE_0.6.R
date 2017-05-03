library(assertthat)

simulate_oneSNP <- function(Xm, delmaf) { ## generate paternal SNP by permuting the maternal SNP
  N <- length(Xm)
  Xf <- Xm[sample(N)] ### But, check expected correlations in this case
  ### DELTA-MAF section: how many paternal transmission to change to make different allele freqs from M and from F?
  ntoswitch <- abs(N * delmaf)
  
  if (delmaf  > 0) {
    switchrows <- sample(which(Xf == 0), ntoswitch)
  } else {
    switchrows <- sample(which(Xf == 1), ntoswitch)
  }
  Xf[switchrows] <- 1 - Xf[switchrows]
  return(Xf)
}

simulate_corHap <- function(Xj, Xk, lcor){
  
  N <- length(Xj)
  assert_that(length(Xj) == length(Xk))
  
  Xjstar <- Xj
  Xkstar <- Xk
  localcor <- lcor
  localdeltacor <- localcor - cor(Xj, Xk)
  nchng <- floor((N - 1) * abs(localdeltacor) * sd(Xj) * sd(Xk) )
  
  # Sanity check 1: Correlation can't be adjusted to be above 1 or below -1
  assert_that( abs(cor(Xj, Xk) + localdeltacor) < 1)
  
  if (nchng > 0) {
    # Which are 1s in both, which are 0s in J, 1s in K, which are 1s in J, 0s in K
    # Xj == Xk
    AkAj <- which(Xj * Xk > 0) ## agreement as 1, can be changed to decrease correlation
    Akaj <- which(Xj < Xk) ## 1s in k, 0s in j, can be changed to increase correlation,
    akAj <- which(Xj > Xk) ## 0s in k, 1s in j, change j-1s to maintain constant pi if correlation increased
    akaj <- which((Xk == 0) & (Xj == 0)) ## 0s in both, change j-0s to maintain constant pi if correlation decreased
    
    if (localdeltacor > 0) {
      #     Increase correlation by changing some of the changeUp values to 1, and maintain pi by changing changeDowns  to 0
      min_len <- min(length(Akaj), length(akAj))
      to_chng <- ifelse(nchng > min_len, min_len, nchng) ## can't have better-than-perfect association at this allele frequency
      
      Xjstar[sample(Akaj, to_chng)] <- 1
      Xjstar[sample(akAj, to_chng)] <- 0
      
    } else {
      # Decrease correlation by changing some agg1 values to 0, maintain pi by changing agg0s to 1s
      min_len <- min(length(akaj), length(AkAj))
      to_chng <- ifelse(nchng > min_len, min_len, nchng) ## can't have better-than-perfect association at this allele frequency
      Xjstar[sample(akaj, to_chng)] <- 1
      Xjstar[sample(AkAj, to_chng)] <- 0
    }
    
    
  }
  
  
  second_deltacor <- localcor - cor(Xk, Xjstar) ## Is it necessary to change both SNPs to achieve desired LD?
  second_nchng <-  floor((N - 1) * abs(second_deltacor) * sd(Xjstar) * sd(Xkstar) )
  
  if (second_nchng > 0) {
    
    AkAjs <- which(Xk * Xjstar > 0) ## individuals sharing the reference allele at Xk and the modified Xj; change to decrease correlation
    akAjs <- which(Xk < Xjstar) ## individuals with nonref at K and ref at modified Xj; convert a_k to A_k to increase correlation
    Akajs <- which(Xk > Xjstar) ## individuals with ref at K and nonref at modified Xj; convert when increasing correlation to maintain equal MAF
    akajs <- which(Xk == 0) && which(Xjstar == 0) ## change to 1s to maintain MAF if correlation must decrease
    
    if (second_deltacor > 1/N) {
      min_len <- min(length(akAjs), length(Akajs)) 
      to_chng <- ifelse(second_nchng > min_len, min_len, nchng) 
      
      Xkstar[sample(akAjs, to_chng)] <- 1
      Xkstar[sample(Akajs, to_chng)] <- 0
      
    } else if (second_deltacor < 1/N) { 
      min_len <- min(length(akajs), length(AkAjs))
      to_chng <- ifelse(second_nchng > min_len, min_len, nchng) 
      
      Xkstar[sample(akajs, to_chng)] <- 1
      Xkstar[sample(AkAjs, to_chng)] <- 0
      
      
    }
  }
  return(list(Xjstar = Xjstar, Xkstar = Xkstar))
  
}

simulate_twoSNP <- function(Xmj, Xmk, hcor, deltamaf = 0){
  Xfj <- simulate_oneSNP(Xm = Xmj, delmaf = deltamaf)
  Xfk <- simulate_oneSNP(Xm = Xmk, delmaf = deltamaf)
  Xms <- simulate_corHap(Xj = Xmj, Xk = Xmk, lcor = hcor)
  Xfs <- simulate_corHap(Xj = Xfj, Xk = Xfk, lcor = hcor)
  Xmks <- Xms$Xkstar
  Xfks <- Xfs$Xkstar
  Xmjs <- Xms$Xjstar
  Xfjs <- Xfs$Xjstar
  
  return(data.frame(Xmj = Xmjs, Xmk = Xmks, Xfj = Xfjs, Xfk = Xfks))
}


simulate_twoSNP_vec <- function(mmatriz, idx, jvec, kvec, corvec, delta_MAF){
  l2snp <- simulate_twoSNP(Xmj = mmatriz[,jvec[idx]], Xmk = mmatriz[,kvec[idx]], hcor = corvec[[idx]], deltamaf = delta_MAF)
  ## mmatriz[,jvec[idx]] <- l2snp$Xmjs
  ## fmatriz[,jvec[idx]] <- l2snp$Xfjs   
  ## fmatriz[,kvec[idx]] <- l2snp$Xfk
  return(l2snp)
}

simulate_zztR <- function(p, i, Xf ,Xm, localcor) {
  ### Consider computing this as the inner product
  Xmi <- Xm[,i]
  Xfi <- Xf[,i]
  localdeltacor <- localcor#[[i]]
  N <- length(Xfi)
  ### ZmZf' section
  
  nchng <- floor((N - 1) * abs(localdeltacor) * sd(Xfi) * sd(Xmi) )
  
  # Sanity check 1: Correlation can't be adjusted to be above 1 or below -1
  assert_that( abs(cor(Xfi, Xmi) + localdeltacor) < 1)
  
  Xfistar <- Xfi
  
  if (nchng > 0) {
    # Which are 1s in both, which are 0s in F, 1s in M, which are 1s in F, 0s in M
    # Xfi == Xmi
    agg1 <- which(Xfi * Xmi > 0) ## agreement as 1, can be changed to decrease correlation
    changeUp <- which(Xfi < Xmi) ## 1s in M, 0s in F, can be changed to increase correlation,
    changeDown <- which(Xfi > Xmi) ## 0s in M, 1s in F, change F-1s to maintain constant pi if correlation increased
    agg0 <- which((Xmi == 0) & (Xfi == 0)) ## 0s in both, change F-0s to maintain constant pi if correlation decreased
    
    if (localdeltacor > 0) {
      #     Increase correlation by changing some of the changeUp values to 1, and maintain pi by changing changeDowns  to 0
      #     # Sanity check 2: no need to change or nchange > number of M-1s?
      assert_that(nchng <= min(length(changeUp), length(changeDown))) #  cant change more than could be changed
      Xfistar[sample(changeUp, nchng)] <- 1
      Xfistar[sample(changeDown, nchng)] <- 0
      
    } else {
      # Decrease correlation by changing some agg1 values to 0, maintain pi by changing agg0s to 1s
      assert_that(nchng <= min(length(agg0), length(agg1))) #  cant change more than could be changed
      
      Xfistar[sample(agg0, nchng)] <- 1
      Xfistar[sample(agg1, nchng)] <- 0
    }
  }
  return(Xfistar)
}

simulate_ztz <- function(N, i, Xf ,Xm, localcor) {
  ### Consider computing this as the inner product
  Xmi <- Xm[i,]
  Xfi <- Xf[i,]
  localdeltacor <- localcor#[[i]]
  p <- length(Xfi)
  ### ZmZf' section
  
  nchng <- floor((p - 1) * abs(localdeltacor) * sd(Xmi)*sd(Xfi) ) # thechance of success must be the same in both cases
  
  # Sanity check 1: Correlation can't be adjusted to be above 1 or below -1
  assert_that( abs(cor(Xfi, Xmi) + localdeltacor) < 1)
  
  Xfistar <- Xfi
  
  if (nchng > 0) {
    # Which are 1s in both, which are 0s in F, 1s in M, which are 1s in F, 0s in M
    # Xfi == Xmi
    agg1 <- which(Xfi * Xmi > 0) ## agreement as 1, can be changed to decrease correlation
    changeUp <- which(Xfi < Xmi) ## 1s in M, 0s in F, can be changed to increase correlation,
    changeDown <- which(Xfi > Xmi) ## 0s in M, 1s in F, change F-1s to maintain constant pi if correlation increased
    agg0 <- which((Xmi == 0) & (Xfi == 0)) ## 0s in both, change F-0s to maintain constant pi if correlation decreased
    
    if (localdeltacor > 0) {
      #     Increase correlation by changing some of the changeUp values to 1, and maintain pi by changing changeDowns  to 0
      #     # Sanity check 2: no need to change or nchange > number of M-1s?
      assert_that(nchng <= min(length(changeUp), length(changeDown))) #  cant change more than could be changed
      Xfistar[sample(changeUp, nchng)] <- 1
      Xfistar[sample(changeDown, nchng)] <- 0
      
    } else {
      # Decrease correlation by changing some agg1 values to 0, maintain pi by changing agg0s to 1s
      assert_that(nchng <= min(length(agg0), length(agg1))) #  cant change more than could be changed
      
      Xfistar[sample(agg0, nchng)] <- 1
      Xfistar[sample(agg1, nchng)] <- 0
    }
  }
  return(Xfistar)
}



sim_zztAll <- function(Xm = Xm, Xf = Xfi, zztcor = zztcor){
  lcm <- zztcor/prod(dim(Xm))  #+ sqrt(0.05)*scale(rnorm(ncol(Xm),mean = 0  ))
  sapply(X = 1:ncol(Xm), FUN = simulate_zztR, p = ncol(Xm), localcor = lcm, Xm = Xm, Xf = Xf)
}

sim_ztzAll <- function(Xm = Xm, Xf = Xfi, zztcor = zztcor){
  lcm <-  zztcor/prod(dim(Xm)) #+ sqrt(0.05)*scale(rnorm(nrow(Xm),mean = 0  ))
  t(sapply(X = 1:nrow(Xm), FUN = simulate_ztz,  N = nrow(Xm), localcor = lcm, Xm = Xm, Xf = Xf))
}



sim_chble <- function(N, p, hsqA, hsqD, hsqP, delta_Maf, new_Cor, avg_LD = 0, pair_Count = 0){
  true_EAFs <- runif(p, min = 1e-2, max = 0.5)
  if (N < 250) { 
    true_EAFs <- rep(0.5, p)
  }
  mat_ped_snps <- sapply(true_EAFs, rbinom, n = N, size = 1)
  pat_ped_snps <- apply(mat_ped_snps,2, simulate_oneSNP, delmaf = delta_Maf)
  
  if (avg_LD != 0) {
    
    npairs <- max(pair_Count, floor(0.10 * 0.5 * (p - p^2)))    
    pair1 <- sample(x = p, size = npairs)
    possible_pairs <- expand.grid(pair1, setdiff(1:p, pair1))
    possible_pairs <- possible_pairs[sample(nrow(possible_pairs), size = npairs),]
    lds <- rnorm(nrow(possible_pairs), mean = avg_LD, sd = 0.05)
    
    for (ti in 1:length(lds)) {
      local2snp <- simulate_twoSNP_vec(mmatriz = mat_ped_snps, idx = ti, jvec = possible_pairs[,1], kvec = possible_pairs[,2], corvec = lds, delta_MAF = delta_Maf)           
      mat_ped_snps[,possible_pairs[ti,1]] <- local2snp$Xmj
      mat_ped_snps[,possible_pairs[ti,2]] <- local2snp$Xmk
      pat_ped_snps[,possible_pairs[ti,1]] <- local2snp$Xfj
      pat_ped_snps[,possible_pairs[ti,2]] <- local2snp$Xfk
    }
  }
  
  ##browser()
  
  betaA <- rnorm(n = p, mean = 0, sd = sqrt(hsqA / p))
  betaD <- rnorm(n = p, mean = 0, sd = sqrt(hsqD / p))
  betaP <- rnorm(n = p, mean = 0, sd = sqrt(hsqP / p))
  
  if ( p >= N ) {
    ## use ztz for speed
    basecor <- sum(diag( t(scale(mat_ped_snps)) %*% scale(pat_ped_snps)))/(N*p)
    zzt <- (new_Cor - basecor)*(N*p)
    pat_ped_snps2 <- sim_ztzAll(Xm = mat_ped_snps, Xf = pat_ped_snps,zztcor = zzt)
    ach_Cor <-  sum(diag(t(scale(mat_ped_snps)) %*% scale(pat_ped_snps2)))/(N*p) - basecor
  } else {
    ## use zzt for speed
    basecor <- sum(diag( scale(mat_ped_snps) %*% t( scale(pat_ped_snps))))/(N*p)
    zzt <- (new_Cor - basecor)*(N*p)
    pat_ped_snps2 <- sim_zztAll(Xm = mat_ped_snps, Xf = pat_ped_snps,zztcor = zzt)
    ach_Cor <-  sum(diag(scale(mat_ped_snps) %*% t(scale(pat_ped_snps2))))/(N*p) - basecor
  }
  
  return(list(true_EAFs = true_EAFs, mat_ped_snps = mat_ped_snps, pat_ped_snps = pat_ped_snps2, base_cor = basecor, ach_Cor = ach_Cor,  betaA = betaA,  betaD = betaD, betaP = betaP))
}


make_mldose_and_ped <- function(nsnps, nobs, mat_ped_snps, pat_ped_snps, runid) {

    sim_ids <- paste0("id", 1:nobs)
    sim_snps <- paste0("rs", 1:nsnps)
    
    ## MLDOSE
    ped_snps <- mat_ped_snps + pat_ped_snps
    dosage <- mat_ped_snps - pat_ped_snps


    ## MLINFO
    mlmeans <- colMeans(ped_snps)/2
    mlA1 <- rep("A", nsnps)
    mlA2 <- rep("G", nsnps)
    mlQ <- runif(n = nsnps, min = 0.5, max = 0.999) 
    mlRsq <- runif(n = nsnps, min = 0.5, max = 0.999)

    mlinfo <- data.frame(sim_snps, mlA1, mlA2, mlmeans, mlQ, mlRsq)

    write.table(data.frame(paste(sim_ids, sim_ids, sep = "->"),"ML_DOSE", dosage), file = gzfile(paste0("simd.", runid, ".mldose.gz")), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(data.frame(sim_snps, mlA1, mlA2, sapply(1:nsnps, FUN = function(i) {min(mlmeans[i], 1- mlmeans[i])}), mlmeans, mlQ, mlRsq), file = gzfile(paste0("simd.", runid, ".mlinfo.gz")), row.names = FALSE, col.names = c("SNP", "Al1", "Al2", "Freq1", "MAF", "Quality", "Rsq"), sep = "\t", quote = FALSE)


    ## MAP
    sim_map <- cbind(1, sim_snps, 0, 10000 + 1:nsnps)
    write.table(sim_map, file = paste0("simd.",runid,".map"), row.names = FALSE, col.names = FALSE, quote = FALSE)


    ## PED
    obs_sexes <- c(rep(1, nobs/2), rep(2, nobs/2))
    ped_leading_cols <- cbind(sim_ids, sim_ids, 0, 0, obs_sexes, -9)
    ped_snps_towrite <- apply(ped_snps, 2, gsub, pattern = 0, replacement = "A A")
    ped_snps_towrite <- apply(ped_snps_towrite, 2, gsub, pattern = 1, replacement = "A G")
    ped_snps_towrite <- apply(ped_snps_towrite, 2, gsub, pattern = 2, replacement = "G G")
    write.table(cbind(ped_leading_cols, ped_snps_towrite), file = paste0("simd.",runid,".ped"), row.names = FALSE, col.names = FALSE, quote = FALSE)

    ## Convert to binary and remove big-n-sloppy plaintext
    system(paste0("plink --file simd.", runid," --recode --make-bed --out simd.", runid))
    system(paste0("rm *",runid,".{ped,map}"))
    
    
}



make_adp_grms <- function(runid){

    ## Additive GRM
    system(paste0("./gcta64 --bfile simd.", runid, " --make-grm-bin --out simA.", runid))
    ## Dominance GRM
    system(paste0("./gcta64 --bfile simd.", runid, " --make-grm-d-bin --out simD.", runid))
    ## Parent-of-origin GRM
    system(paste0("./gcta64 --dosage-mach-gz simd.", runid, ".mldose.gz simd.", runid, ".mlinfo.gz --make-grm-bin --out simP.",runid))


}

