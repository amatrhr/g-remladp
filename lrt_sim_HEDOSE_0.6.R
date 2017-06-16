## R script to simulate parent-of-origin effects for analysis with Haseman-Elston regression
source("./simulate_funcs_HEDOSE_0.6.R")

## Called from the commmand line with "Rscript"
args <- as.numeric(commandArgs(TRUE))

#### Test values of the arguments
## args <- c(91433, 0.1, 0.01, 0.05, 0.00, 0.00, 0.25, 10, 1800, 2601, 5); set.seed(seed = args[[1]])

heSim <- function(args){
    ## Wrapper function to generate simulated data sets by calling several other funcitons OB
    runid <- (args[[1]] + floor(1e9*runif(1)) %% as.numeric(format(Sys.time(),"%s")))
    set.seed(seed = runid) 
    hsqA <- args[[2]] ## proportion of variance attributable to additive effects
    hsqD <- args[[3]] ## proportion of variance attributable to dominance effects
    hsqP <- args[[4]] ## proportion of variance attributable to parent-of-origin effects

    delta_Maf <- args[[5]] ## difference in minor allele frequencies between mothers and fathers
    rel_Cor <- args[[6]] ## average correlation between parents' genotypes
    avg_LD <- args[[7]] ## (average LD, for pairs of SNPs in LD )
    npairvar <- args[[8]] ## number of SNP pairs for which to simulate LD > 0
    nobs <- args[[9]] ## Number of simulated probands
    nsnps <- args[[10]] ## Number of simulated SNPsOB
    


    simulate_struct <- sim_chble(N = nobs, p = nsnps, hsqA = hsqA, hsqD = hsqD, hsqP = hsqP, delta_Maf = delta_Maf, new_Cor = rel_Cor, avg_LD = avg_LD, pair_Count = npairvar) ## object holding sim data

    ## Generate proband's genotype from simulated transmissionsOB
    Aeps <- (simulate_struct$mat_ped_snps + simulate_struct$pat_ped_snps)
    Deps <- Aeps * matrix(1, nrow = nrow(Aeps)) %*% matrix(colMeans(Aeps), nrow = 1) - 2*(Aeps > 1)
    Xeps <- (simulate_struct$mat_ped_snps - simulate_struct$pat_ped_snps)

### Standardize the simulated genotype codingsOB
    aSD <-  function(pv, qv){
        sqrt(2 * pv * qv)
    }

    Dmean <- function(pv, qv){
        2 * pv^2
    }

    DSD <- function(pv, qv){
        sqrt(4 * pv^2 * qv^2)
    }
    
    p <- 0.5 * colMeans(Aeps)
    q <- 1 - p

    ## summary(data.frame(Aeps, Deps, Xeps))
    ## crossprod(cbind(Aeps, Deps, Xeps))
    ## cov(data.frame(Aeps, Deps, Xeps))

    
    ## Scale using understood ways
    Aeps <- sweep( sweep(Aeps, 2, 2 * p), 2, aSD(p,q), FUN = "/") 
    Deps <- sweep( sweep(Deps, 2, Dmean(p,q)), 2, DSD(p,q), FUN = "/") 
    Xeps <- sweep( Xeps, 2, aSD(p,q), FUN = "/") 

    ## summary(data.frame(Aeps, Deps, Xeps))
    ## crossprod(cbind(Aeps, Deps, Xeps))
    ## cov(data.frame(Aeps, Deps, Xeps))

    
    aGRM <- (1 / nsnps) * tcrossprod(Aeps)
    dGRM <- (1 / nsnps) * tcrossprod(Deps)

    xGRM <- (1 / nsnps) * tcrossprod(Xeps)

    ## gendata <- data.frame(add = as.vector(aGRM), dom = as.vector(dGRM), xeps = as.vector(xGRM))
    ## summary(gendata)
    ## XtX <- crossprod(data.matrix(gendata))
    ## signif((XtX - diag(nrow(gendata), 3))/nrow(gendata) * 100, 3)
    ## cov(gendata) - XtX/(nrow(gendata))
    
### two ways of phenotype generation  
    pheno <- Aeps %*% matrix(simulate_struct$betaA, ncol = 1) +  Deps %*% matrix(simulate_struct$betaD, ncol = 1)  + Xeps %*% matrix(simulate_struct$betaP, ncol = 1) + rnorm(nobs, mean = 0, sd = sqrt( 1 - (hsqA + hsqP + hsqD)))
    write.table(cbind(paste0("id", 1:nobs),  paste0("id", 1:nobs),  pheno), file = paste0("simd.",runid,".pheno"), row.names = FALSE, col.names = FALSE, quote = FALSE)

### R computations: making GRM

    mGRM <- (1 / nsnps) * tcrossprod(scale(simulate_struct$mat_ped_snps))
    pGRM <- (1 / nsnps) * tcrossprod(scale(simulate_struct$pat_ped_snps))

    phGRM <- tcrossprod(scale(pheno, scale = FALSE)) ## Merely center the phenotype
    
                                        # XGRMalt <- (1 / nsnps) * tcrossprod(scale(Xeps)*(nrow(Xeps) - 1) /nrow(Xeps)) 
                                        # plot(XGRM[lower.tri(XGRM)], XGRMalt[lower.tri(XGRMalt)])
    ecoef <- sum(simulate_struct$betaP^2) ## expected beta valuesOB

    ## Data frame for Haseman-Elston regression
    hedata <- data.frame(pheno = phGRM[lower.tri(phGRM)], mgrm = mGRM[lower.tri(mGRM)], pgrm = pGRM[lower.tri(pGRM)], add = aGRM[lower.tri(aGRM)], dom = dGRM[lower.tri(dGRM)], xeps = xGRM[lower.tri(xGRM)])

    ## signif((crossprod(data.matrix(hedata)) - diag(nrow(hedata), 6))/nrow(hedata) * 100, 3) 

    ### Three Haseman-Elston Models: 
    #### Using parentally transmitted GRMS
    HE1.2 <- lm(pheno ~ mgrm + pgrm, data = hedata)
    #### All three (proband) coded genotype GRMs
    HE2.3 <- lm(pheno ~ add + dom + xeps, data = hedata)
    #### Additive and dominance coded GRMs only
    HE2.4 <- lm(pheno ~ add + dom, data = hedata)

    with(simulate_struct, make_mldose_and_ped(mat_ped_snps = mat_ped_snps, pat_ped_snps = pat_ped_snps, runid = runid, nobs = nobs, nsnps = nsnps ))
    make_adp_grms(runid)

    ## Store results from the simulation in the file system
    system(paste0("touch simgrmlist.", runid, ".txt"))
    system(paste0("echo simA.", runid, ">> simgrmlist.", runid, ".txt"))
    system(paste0("echo simD.", runid, ".d >> simgrmlist.", runid, ".txt"))
    system(paste0("echo simP.", runid, ">> simgrmlist.", runid, ".txt"))
    
    varGCTA <- system(paste0("./gcta64 --reml --reml-no-constrain --mgrm-bin simgrmlist.", runid,".txt --pheno simd.", runid,".pheno --out simd.", runid), intern = TRUE)

    gcta_result <- vector(mode = "numeric", length = 18)
    gcta_result[1:2] <- as.numeric(unlist(strsplit(grep(pattern = "V.G1./Vp", x = varGCTA, value=TRUE), split = "\t"))[2:3])
    gcta_result[3:4] <- as.numeric(unlist(strsplit(grep(pattern = "V.G2./Vp", x = varGCTA, value=TRUE), split = "\t"))[2:3])
    gcta_result[5:6] <- as.numeric(unlist(strsplit(grep(pattern = "V.G3./Vp", x = varGCTA, value=TRUE), split = "\t"))[2:3])
    gcta_result[[7]] <- as.numeric(unlist(strsplit(system(paste0('fgrep -w LRT simd.',runid, '.hsq'), intern = TRUE), "\t"))[[2]])
    gcta_result[[8]] <- as.numeric(unlist(strsplit(system(paste0('fgrep -w Pval simd.',runid, '.hsq'), intern = TRUE), "\t"))[[2]])
    gcta_result[9:18] <- matrix(as.numeric(unlist(strsplit(x = varGCTA[91:94], split = "\t"))), nrow = 4, ncol = 4)[lower.tri(diag(4), diag = TRUE)]

    names(gcta_result) <- c("gvA","segA", "gvD", "segD", "gvP", "segP", "LRT", "gpval", "vgvA", "gcAD", "gcAP", "gcAE", "vgvD", "gcDP", "gcDE", "vgvP", "gcPE", "vgvE")
    
    result <- c(nobs, nsnps, hsqA, hsqD, hsqP, delta_Maf, rel_Cor, avg_LD, npairvar, summary(HE2.3)$fstatistic, anova(HE2.4, HE2.3)$F[[2]], summary(HE2.3)$coefficients['xeps',"Pr(>|t|)"], coef(HE2.3)[-1],sqrt(diag(vcov(HE2.3)))[-1], coef(HE1.2)[-1], sqrt(diag(vcov(HE1.2)))[-1], cov(hedata)[lower.tri(diag(6), diag = TRUE)])

    names(result) <- c("nobs", "nsnps", "hsqA", "hsqD", "hsqP", "delta_Maf", "rel_Cor", "avg_LD", "npairvar", "threeF", "threeFndf", "threeFddf", "twoF","Pp", "heAdd", "heDom", "hePoo", "seheAdd", "seheDom", "sehePoo",  "heMo", "heFa", "seheMo", "seheFa", "phPiV","PMC", "PFC", "PAC","PDC", "PEC", "moPiV","MPC","MAC","MDC","MEC","faPiv","FAC","FDC","FEC","adPiV","ADC","AEC","doPiV","DEC","xePiV")
    
    result <- c(result, colMeans(hedata))
    system(paste0("rm *", runid, "*"))
    
    return(c(result, gcta_result))

}

## Perform n replications of the simulation under conditions given in args vector at command line
outcome <- t(replicate(heSim(args), n = args[[11]]))

write.table(outcome, file = paste0("res", paste0(args, collapse = "V"), "X.txt"), row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
