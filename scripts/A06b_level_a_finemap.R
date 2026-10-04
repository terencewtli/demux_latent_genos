# A06b: Level A smoke test (docs/A06_FINEMAP_PLAN.md), one gene region.
# Genotype arms are built by injecting INDEPENDENT per-entry noise into the true dosages, calibrated to the measured
# GLIMPSE2 / naive accuracy for each donor x SNP site class (A06a; A04g summary_depth / summary_run):
#   *_imp   posterior-mean (Berkson-type) noise: Ghat = mu + r2 (G - mu) + eta, var(eta) = r2 (1 - r2) var(G),
#           so corr(G, Ghat)^2 = r2 and E[G | Ghat] = Ghat (what an imputed dosage is)
#   gex_naive  hard calls at observed sites only (class >= 2): classical noise at the naive r2 of the depth bin, rounded
#           to 0/1/2; unobserved entries mean-imputed per SNP (what using souporcell's calls without imputation does)
#   multiome_imp  per entry, the better of the GEX and ATAC r2 (approximates summing the two modalities' counts)
# Real imputation error is correlated along LD blocks; Level B (GLIMPSE2 on emulated GLs) tests that.
# eQTL: donor-level pseudobulk y = G_true b + e, cis-h2 in {0.05, 0.10, 0.20}; 1 causal (70%) or 2 (30%); causal
# class scenario: any SNP / untyped (GEX-observed in < 5% of donors) / observed (GEX-observed in >= 50% of donors).
# Fine-mapping: susieR (L = 10). Coloc (1-causal phenotypes): GWAS z ~ MVN(lambda R[, c], R) from the true-genotype
# LD, either at the eQTL causal (shared) or at a different SNP with r2 0.3-0.8 to it (distinct); coloc.abf.
# Output: <outdir>/<region>.tsv.gz, one row per phenotype x arm.
# Usage: Rscript A06b_level_a_finemap.R <indir> <region> <outdir> [n_seeds=3] [level=A|B]
suppressPackageStartupMessages({
    library(susieR)
    library(coloc)
    library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]; region <- args[2]; outdir <- args[3]
n_seeds <- if (length(args) >= 4) as.integer(args[4]) else 3
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(outdir, paste0(region, '.tsv.gz'))
if (file.exists(out_file)) { cat(region, 'done already\n'); quit(save = 'no') }
t0 <- Sys.time()

rd <- file.path(indir, region)
snps <- fread(file.path(rd, 'snps.tsv'))
rd_mat <- function(f) { m <- as.matrix(fread(file.path(rd, f), header = FALSE)); dimnames(m) <- NULL; m }
G <- rd_mat('truth.tsv.gz'); storage.mode(G) <- 'double'
CG <- rd_mat('class_gex.tsv.gz'); CA <- rd_mat('class_atac.tsv.gz')
n <- nrow(G); p <- ncol(G)
cat(region, n, 'donors', p, 'SNPs\n')

# measured accuracy by class 0..5 (untyped, listed depth 0, 1-2, 3-5, 6-10, >10)
r2_imp <- list(gex = c(0.556, 0.916, 0.937, 0.948, 0.946, 0.918), atac = c(0.744, 0.852, 0.887, 0.923, 0.949, 0.937))
r2_naive_gex <- c(NA, NA, 0.604, 0.699, 0.765, 0.849)

mu <- colMeans(G); vg <- pmax(apply(G, 2, var), 1e-8)
mu_m <- matrix(mu, n, p, byrow = TRUE); sd_m <- matrix(sqrt(vg), n, p, byrow = TRUE)
posterior_arm <- function(R2) {
    eta <- matrix(rnorm(n * p), n, p) * sqrt(R2 * (1 - R2)) * sd_m
    pmin(pmax(mu_m + R2 * (G - mu_m) + eta, 0), 2)
}
set.seed(1000 + as.integer(sub('r', '', region)))
R2g <- matrix(r2_imp$gex[CG + 1], n, p); R2a <- matrix(r2_imp$atac[CA + 1], n, p)
naive <- function() {
    obs <- CG >= 2
    R2 <- matrix(r2_naive_gex[pmax(CG, 2) + 1], n, p)
    X <- round(G + matrix(rnorm(n * p), n, p) * sd_m * sqrt((1 - R2) / R2))
    X <- pmin(pmax(X, 0), 2); X[!obs] <- NA
    cm <- colMeans(X, na.rm = TRUE); cm[is.na(cm)] <- 0
    idx <- which(is.na(X), arr.ind = TRUE); X[idx] <- cm[idx[, 2]]
    X
}
arms <- list(truth = G, gex_imp = posterior_arm(R2g), atac_imp = posterior_arm(R2a),
             multiome_imp = posterior_arm(pmax(R2g, R2a)), gex_naive = naive())
# Level B (5th argument 'B'): imputed arms are the GLIMPSE2 dosages from emulated GLs (A06d / A06e), read from
# <region>/imp_<arm>.tsv.gz; truth and gex_naive as in Level A. Phenotypes, GWAS and seeds are identical to Level A,
# because they are reseeded per phenotype below.
level <- if (length(args) >= 5) args[5] else 'A'
if (level == 'B') {
    for (a in c('gex', 'atac', 'multiome')) {
        X <- rd_mat(paste0('imp_', a, '.tsv.gz')); storage.mode(X) <- 'double'
        stopifnot(nrow(X) == n, ncol(X) == p)
        arms[[paste0(a, '_imp')]] <- X
    }
}
# per-SNP r2 over SNPs polymorphic in the truth; an arm column that is constant there (no information, e.g. naive
# mean-imputed sites or Level B '.none' regions) counts as r2 = 0 instead of being dropped as NaN
poly <- which(apply(G, 2, function(g) var(g) > 0))
realized_r2 <- sapply(arms, function(X) mean(sapply(poly, function(j) {
    r <- suppressWarnings(cor(X[, j], G[, j]))
    if (is.na(r)) 0 else r^2
})))
cat('realized mean per-SNP r2:', paste(names(realized_r2), round(realized_r2, 3), collapse = ' '), '\n')

frac_obs_gex <- colMeans(CG >= 2); frac_obs_atac <- colMeans(CA >= 2)
pools <- list(any = seq_len(p), untyped = which(frac_obs_gex < 0.05), observed = which(frac_obs_gex >= 0.5))
cat('causal pools:', paste(names(pools), lengths(pools), collapse = ' '), '\n')

Gs <- scale(G); Gs[is.na(Gs)] <- 0
LD <- crossprod(Gs) / (n - 1)
Lchol <- chol(LD + diag(1e-3, p))
maf <- pmin(snps$af, 1 - snps$af)
marginal <- function(X, y) {
    Xc <- scale(X, scale = FALSE); yc <- y - mean(y)
    sxx <- colSums(Xc^2); b <- drop(crossprod(Xc, yc)) / sxx
    res <- sum(yc^2) - b^2 * sxx
    se <- sqrt(pmax(res, 1e-12) / (n - 2) / sxx)
    list(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

rows <- list()
for (scen in names(pools)) {
    pool <- pools[[scen]]
    if (length(pool) < 5) next
    for (h2 in c(0.05, 0.10, 0.20)) for (seed in seq_len(n_seeds)) {
        set.seed(seed * 7919 + round(h2 * 100) * 31 + match(scen, names(pools)))
        k <- if (runif(1) < 0.7) 1 else 2
        causal <- sample(pool, k)
        b <- rnorm(k)
        g <- drop(Gs[, causal, drop = FALSE] %*% b)
        g <- g / sd(g) * sqrt(h2)
        y <- g + rnorm(n, sd = sqrt(1 - h2))
        # GWAS for coloc (1-causal phenotypes): shared or distinct causal, lambda = 6
        gw <- NULL
        if (k == 1) {
            c1 <- causal[1]; r2c <- LD[c1, ]^2
            cand <- setdiff(which(r2c >= 0.3 & r2c <= 0.8), c1)
            for (type in c('shared', 'distinct')) {
                cg <- if (type == 'shared') c1 else if (length(cand)) cand[sample.int(length(cand), 1)] else NA
                if (is.na(cg)) next
                z <- 6 * sign(b[1]) * LD[, cg] + drop(crossprod(Lchol, rnorm(p)))
                gw[[type]] <- list(z = z, causal = cg)
            }
        }
        for (arm in names(arms)) {
            X <- arms[[arm]]
            keep <- apply(X, 2, var) > 1e-8
            Xk <- X[, keep, drop = FALSE]; idx <- which(keep)
            fit <- tryCatch(susie(Xk, y, L = 10, max_iter = 200, verbose = FALSE), error = function(e) NULL)
            pip <- rep(0, p); cs_list <- list()
            if (!is.null(fit)) {
                pip[idx] <- fit$pip
                if (!is.null(fit$sets$cs)) cs_list <- lapply(fit$sets$cs, function(s) idx[s])
            }
            in_cs <- sapply(causal, function(c) any(sapply(cs_list, function(s) c %in% s)))
            lead <- which.max(pip)
            mg <- marginal(X, y); mg$p[!keep] <- 1
            row <- data.table(region = region, scenario = scen, h2 = h2, seed = seed, n_causal = k, arm = arm,
                              n_snps = p, susie_ok = !is.null(fit), n_cs = length(cs_list),
                              cs_size_median = if (length(cs_list)) median(lengths(cs_list)) else NA_real_,
                              causal_in_cs = mean(in_cs), all_causal_in_cs = all(in_cs),
                              causal_pip_mean = mean(pip[causal]), causal_pip_gt05 = mean(pip[causal] > 0.5),
                              lead_is_causal = lead %in% causal,
                              lead_ld_causal = max(LD[lead, causal]^2),
                              lead_frac_obs_gex = frac_obs_gex[lead], causal_frac_obs_gex = mean(frac_obs_gex[causal]),
                              lead_observed_gex = frac_obs_gex[lead] >= 0.5,
                              min_p = min(mg$p), egene_bonf = min(mg$p) < 0.05 / p,
                              causal_beta_ratio = if (k == 1) mg$b[causal] / (drop(lm.fit(cbind(1, G[, causal]), y)$coefficients[2])) else NA_real_,
                              coloc_type = NA_character_, pp_h3 = NA_real_, pp_h4 = NA_real_)
            if (length(gw)) {
                for (type in names(gw)) {
                    d1 <- list(beta = mg$b[keep], varbeta = mg$se[keep]^2, type = 'quant', N = n, sdY = sd(y),
                               snp = as.character(idx))
                    zz <- gw[[type]]$z[keep]; Ng <- 50000; vv <- 2 * maf[keep] * (1 - maf[keep])
                    d2 <- list(beta = zz / sqrt(Ng * vv), varbeta = 1 / (Ng * vv), type = 'quant', N = Ng, sdY = 1,
                               snp = as.character(idx))
                    cc <- tryCatch(suppressMessages(coloc.abf(d1, d2))$summary, error = function(e) NULL)
                    r2 <- copy(row); r2$coloc_type <- type
                    if (!is.null(cc)) { r2$pp_h3 <- cc[['PP.H3.abf']]; r2$pp_h4 <- cc[['PP.H4.abf']] }
                    rows[[length(rows) + 1]] <- r2
                }
            } else rows[[length(rows) + 1]] <- row
        }
    }
}
res <- rbindlist(rows)
res[, realized_r2 := realized_r2[arm]]
fwrite(res, out_file, sep = '\t')
cat(region, nrow(res), 'rows;', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), 'min\n')
