#!/bin/bash

### res10V0.2V0.1V0.05V0.05V0V0.1V24V3000V4000V10X.txt
awk '$0 !~/nobs/{print $0}' res*.txt > combined_results_faster_HEDOSE_0.txt
awk 'BEGIN { print "nobs", "nsnps", "hsqA", "hsqD", "hsqP", "delta_Maf", "rel_Cor", "avg_LD", "npairvar", "threeF", "threeFndf", "threeFddf", "twoF","Pp", "heAdd", "heDom", "hePoo", "seheAdd", "seheDom", "sehePoo",  "heMo", "heFa", "seheMo", "seheFa", "phPiV","PMC", "PFC", "PAC","PDC", "PEC", "moPiV","MPC","MAC","MDC","MEC","faPiv","FAC","FDC","FEC","adPiV","ADC","AEC","doPiV","DEC","xePiV", "pheno", "mgrm", "pgrm", "add", "dom", "xeps", "gvA","segA", "gvD", "segD", "gvP", "segP", "LRT", "gpval", "vgvA", "gcAD", "gcAP", "gcAE", "vgvD", "gcDP", "gcDE", "vgvP", "gcPE", "vgvE"}; {print $0}' combined_results_faster_HEDOSE_0.txt > combined_results_faster_HEDOSE.txt

