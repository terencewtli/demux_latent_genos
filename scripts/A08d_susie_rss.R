# A08d: SuSiE-RSS on published OneK1K z-scores with LD from four genotype sources (A08d_prep_loci.py output).
#
# Per locus (results/A08_eqtl_anchor/finemap/<ct>/<gene>/): z.tsv + geno_<arm>.tsv.gz (samples x variants dosage).
# Arms: onek1k_all (in-sample LD, ~1,098; the "gold"), onek1k_pilot and latent_pilot (same pilot donors; the
# difference is imputation), kg_eur (1000G EUR founders; out-of-sample reference LD).
# LD = Pearson correlation of dosages; monomorphic / missing -> 0 off the diagonal; then R <- (1 - lambda) R + lambda I
# (default lambda 0.1 for every arm; the pilot arms have n ~ 60, so their R is low-rank without it).
# susie_rss(z, R, n = --n, L = 10, coverage 0.95, min_abs_corr 0.5); estimate_s_rss() records LD-z mismatch per arm.
# Outputs per locus: susie_pip.tsv.gz (arm, variant, pip, cs) and susie_summary.tsv (one row per arm).
#
# usage: Rscript scripts/A08d_susie_rss.R --ct cd4nc [--lambda 0.1] [--n 982] [--L 10]

suppressPackageStartupMessages({
    library(susieR)
    library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
ct <- opt$ct %||% stop("--ct required")
lambda <- as.numeric(opt$lambda %||% "0.1")
n_gwas <- as.numeric(opt$n %||% "982")
L <- as.integer(opt$L %||% "10")
arms <- c("onek1k_all", "onek1k_pilot", "latent_pilot", "kg_eur")

proj <- "/u/project/cluo/terencew/claude/project_ideas/latent_genos"
base <- file.path(proj, "results/A08_eqtl_anchor/finemap", ct)
loci <- list.dirs(base, recursive = FALSE)
cat(sprintf("%s: %d loci, lambda %.3g, n %d, L %d\n", ct, length(loci), lambda, n_gwas, L))

ld_matrix <- function(path, variants) {
    g <- fread(path)
    X <- as.matrix(g[, -1, with = FALSE])
    stopifnot(identical(colnames(X), variants))
    R <- suppressWarnings(cor(X, use = "pairwise.complete.obs"))
    R[!is.finite(R)] <- 0
    diag(R) <- 1
    list(R = (1 - lambda) * R + lambda * diag(nrow(R)), n_samples = nrow(X),
         n_mono = sum(apply(X, 2, function(x) var(x, na.rm = TRUE)) < 1e-8, na.rm = TRUE))
}

t0 <- Sys.time()
for (d in loci) {
    out_sum <- file.path(d, "susie_summary.tsv")
    if (file.exists(out_sum)) next
    zt <- fread(file.path(d, "z.tsv"))
    lead <- zt$variant[zt$is_lead == 1][1]
    pips <- list(); sums <- list()
    for (arm in arms) {
        f <- file.path(d, sprintf("geno_%s.tsv.gz", arm))
        if (!file.exists(f)) next
        ld <- ld_matrix(f, zt$variant)
        s_rss <- tryCatch(estimate_s_rss(zt$z, ld$R, n = n_gwas), error = function(e) NA_real_)
        fit <- tryCatch(susie_rss(z = zt$z, R = ld$R, n = n_gwas, L = L, coverage = 0.95, min_abs_corr = 0.5,
                                  max_iter = 200),
                        error = function(e) { message(sprintf("  %s %s: %s", basename(d), arm, conditionMessage(e))); NULL })
        if (is.null(fit)) {
            sums[[arm]] <- data.table(arm = arm, ok = 0L, n_samples = ld$n_samples, n_mono = ld$n_mono, s_rss = s_rss)
            next
        }
        cs <- fit$sets$cs %||% list()
        cs_id <- integer(nrow(zt))
        for (k in seq_along(cs)) cs_id[cs[[k]]] <- k
        pips[[arm]] <- data.table(arm = arm, variant = zt$variant, pip = fit$pip, cs = cs_id)
        lead_i <- which(zt$variant == lead)
        sums[[arm]] <- data.table(arm = arm, ok = 1L, converged = as.integer(isTRUE(fit$converged)),
                                  n_samples = ld$n_samples, n_mono = ld$n_mono, s_rss = s_rss, n_var = nrow(zt),
                                  n_cs = length(cs), cs_sizes = paste(lengths(cs), collapse = ","),
                                  lead_pip = fit$pip[lead_i], lead_in_cs = as.integer(cs_id[lead_i] > 0),
                                  top_pip_variant = zt$variant[which.max(fit$pip)], max_pip = max(fit$pip))
    }
    if (length(pips)) fwrite(rbindlist(pips), file.path(d, "susie_pip.tsv.gz"), sep = "\t")
    fwrite(rbindlist(sums, fill = TRUE)[, `:=`(ct = ct, gene = basename(d), lead = lead)], out_sum, sep = "\t")
}
cat(sprintf("done %s in %.1f min\n", ct, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
