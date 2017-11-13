#PATHFINDER
#Inputs:
#1) Associations file, no header (columns: (1) SNP ids, (2) mark ids, (3) SNP-mark associations, (4) mark-exp associations)
#2) SNP correlations (LD) file (SNPs x SNPs)
#3) Mark correlations file (marks x marks)
#4) Output filename
#5) Prior parameter - variance explained by marks on expression
#Output:
#Posterior probabilities (columns: (1) SNP ids, (2) mark ids, (3) posterior probability for path containing given SNP and mark)

#Process inputs
library(matrixStats, lib.loc=".")
args = commandArgs(trailingOnly=TRUE)
filename = args[1]
filename_corr_snps = args[2]
filename_corr_peaks = args[3]
results_file = args[4]
variance_h = as.numeric(args[5])
Z_table = read.table(filename, header=FALSE)
ld_snps = as.matrix(read.table(filename_corr_snps, header=FALSE))
ld_peaks = as.matrix(read.table(filename_corr_peaks), header=FALSE)
snps = unique(Z_table[,1])
peaks = unique(Z_table[,2])
n.snps = nrow(ld_snps)
n.peaks = nrow(ld_peaks)
if (n.snps<2) {next}
if (n.peaks<2) {next}


#Construct Z_g and Z_h vectors, containing SNP-mark and mark-exp associations
Z_h = c()
Z_g = matrix(0, n.peaks, n.snps)
for (p in 1:n.peaks) {
        for (s in 1:n.snps) {
                snp = Z_table[s+n.snps*(p-1),1]
                peak = Z_table[s+n.snps*(p-1),2]
                z_snp = Z_table[s+n.snps*(p-1),3]
                z_peak = Z_table[s+n.snps*(p-1),4]
                Z_g[p,s]= z_snp
        }
        Z_h = c(Z_h, z_peak)
}

#Compute P(Z_g,Z_h|C) for all possible SNP-mark combinations
X = t(Z_g)
M = matrix(0, nrow=n.snps, ncol=n.peaks)
U = ld_snps+0.01*diag(n.snps)
V = ld_peaks+0.01*diag(n.peaks)
U_inv = solve(U)
V_inv = solve(V)
numerators = c()
for (p in 1:n.peaks) {
        c_h = rep(0, n.peaks)
        c_h[p]=1
        epsilon = rep(0.00001, n.peaks)
        sigma_c_h = (variance_h)*diag(c_h) + diag(epsilon)
        mat = solve(ld_peaks+ld_peaks%*%sigma_c_h%*%ld_peaks+0.0001*diag(n.peaks))
        log_numerator_a = log(exp((-1/2)*t(Z_h)%*%mat%*%Z_h)) #log P(Z_h|C)
        for (s in 1:n.snps) {
                M = matrix(0, nrow=n.snps, ncol=n.peaks)
                lambda_g_c = rep(0, n.snps)
                lambda_g_c[s] = Z_g[p,s]
                for (q in 1:n.peaks) {
                r = ld_peaks[p,q]
                M[,q] = r*ld_snps%*%lambda_g_c
                }
                term = V_inv%*%t(X-M)%*%U_inv%*%(X-M)
                log_numerator_b = -0.5*sum(diag(term)) #log P(Z_g|C)
                log_numerator = log_numerator_a+log_numerator_b #log P(Z_h|C)*P(Z_g|C)
                numerators = c(numerators, log_numerator)
        }
}
#Compute P(Z_g,Z_h)
denominator = logSumExp(numerators)
#Compute P(C|Z_g,Z_h) for all possible SNP-mark combinations
posteriors = exp(numerators - denominator)


#Output posterior probilities
snps = Z_table[,1]
peaks = Z_table[,2]
mat = cbind(snps, peaks)
mat = cbind(mat, posteriors)
write.table(mat, file=results_file, row.names=FALSE, quote=FALSE)


