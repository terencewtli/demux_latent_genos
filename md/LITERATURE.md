# Literature review: latent genotypes from genotype-free demultiplexers, and imputing them

Compiled 2026-09-29. Searched PubMed (E-utilities), Europe PMC (full-text search, which includes PMC and indexed preprints), the bioRxiv API (abstracts only; bioRxiv full-text pages were rate-limited/blocked, so preprint methods were not always readable), general web search, and GitHub issue search (souporcell, vireo, popscle, Demuxafy docs, scSplit, demuxalot, Monopogen repos).

Legend: **[V]** = the paper exists and the claim was checked against the abstract or full text. **[V-abs]** = exists, but only the abstract (or a secondary summary) was read, so method details are not confirmed. **[UNVERIFIED]** = could not confirm.

---

## Bottom line: has this been done?

**(b) Imputing latent genotypes with a reference panel and then re-demultiplexing: I found no evidence it has been done.** No paper, preprint, GitHub issue or forum post I found feeds souporcell `cluster_genotypes.vcf`, vireo `GT_donors.vireo.vcf`, freemuxlet or scSplit cluster genotypes into Minimac4/TOPMed, Beagle, GLIMPSE or IMPUTE5 and then runs demuxlet, vireo-with-genotypes or ambimux on the result.

Searches that came back empty:
- Europe PMC full text for `"cluster_genotypes.vcf"` (0 hits) and `"GT_donors" AND imputation` (0 hits).
- `souporcell AND (minimac4 OR Beagle OR GLIMPSE OR "imputation server") AND ("cluster genotypes" OR "inferred genotypes")` (0 hits).
- `freemuxlet AND imputation AND "reference panel"`. The hits were only papers that impute *array* genotypes.
- GitHub issues for "imputation/impute/imputed" in souporcell, vireo, popscle and Demuxafy. The only hits were about using imputed *array* genotypes as `--known_genotypes`, or about matching clusters to them.

Two caveats limit how confident this can be. bioRxiv full text could not be read directly, and GitHub *code* search was not available without authentication. A lab could have done this privately or in a supplement without saying so in any indexed text.

**(a) Checking latent-genotype accuracy against truth: done, but only briefly, as a side analysis inside tool papers.**
- scSplit reported cluster-genotype concordance of about 98.4%, against 95.9% for genotypes built from demuxlet clusters.
- souporcell compared cluster genotypes with WGS in one supplementary panel (Fig. S1i–j).
- vireo reported that genotypes are accurate only where a donor has more than about 10 reads (Fig. S4).

No benchmark I found treats latent-genotype quality (concordance or dosage R² as a function of MAF, depth, donor number or modality) as a main outcome. The large benchmarks (Demuxafy; Li et al. 2025 eLife; Fu et al. 2025 BiB; Weber 2021; Souporcell3) score droplet assignment, not genotypes. Li et al. 2025 say explicitly that for genotype-free methods "latent genotype inference is the driving force" of errors, but they did not measure it.

**The closest prior work is next door, not on this question.** Several groups impute genotypes from single-cell reads of **one donor per library**:
- Monopogen (Beagle-HMM LD refinement).
- Mu et al. 2026, Cell Genomics: GLIMPSE on aggregated scATAC reads, validated against array+Minimac4.
- Zhang et al. 2026, pig sc-eQTL: cellSNP-lite + SHAPEIT2/Beagle, CR > 0.90, R² 0.80–0.90 against WGS.

The only step our idea adds is running this on *per-cluster pseudo-bulk from a genotype-free pool*, then closing the loop by re-demultiplexing. That step is small conceptually but clearly unoccupied. It is also exactly what both Kulhankova et al. 2023 ("Genetic imputation could further augment the data…") and Alvarez et al. 2025 (ambimux: "an extension of a genotype-free or missing-genotype framework would further expand the utility") name as future work.

---

## 1. Genotype-free demultiplexers and what they said about latent genotype accuracy

**Xu J, Falconer C, Nguyen Q, … Coin LJM (2019). Genotype-free demultiplexing of pooled single-cell RNA-seq (scSplit). *Genome Biology* 20:290. doi:10.1186/s13059-019-1852-7** [V]
- HMM/EM on allele fractions, tested on pools of 2–32 samples. It builds presence/absence genotype matrices at "distinguishing variants" from its clusters and compares them with known genotypes. Concordance was about 98.4% for scSplit clusters vs 95.9% for demuxlet clusters (Table S1).
- *Relevance:* this is the earliest direct latent-vs-truth comparison, but only at a small set of distinguishing sites and only as presence/absence. No imputation.

**Heaton H, Talman AM, Knights A, … Lawniczak MKN (2020). Souporcell: robust clustering of single-cell RNA-seq data by genotype without reference genotypes. *Nature Methods* 17:615–620. doi:10.1038/s41592-020-0820-1** [V]
- Sparse mixture clustering plus joint inference of ambient RNA and cluster genotypes. It compares the cluster genotypes of souporcell, vireo and scSplit against WGS variant calls on a synthetic mixture: "genotype accuracy is significantly lower than one would attain with genome sequencing". The main error mode for vireo and scSplit is hom-ref sites called as het, which is attributed to ambient RNA (Fig. S1i–j).
- *Relevance:* this directly motivates imputation, since LD could fix exactly this kind of error. It is only one supplementary panel, with no stratification by MAF or depth and no dosage R².

**Huang Y, McCarthy DJ, Stegle O (2019). Vireo: Bayesian demultiplexing of pooled single-cell RNA-seq data without genotype reference. *Genome Biology* 20:273. doi:10.1186/s13059-019-1865-2** [V]
- Variational Bayes with a genotype prior. It reports that estimated donor genotypes are accurate only where a donor has more than about 10 reads (Additional file 1, Fig. S4). It supports a **two-step run**: first genotype-free, then replace the inferred genotype probabilities with known values for the matched donors and use that mixed matrix as the prior. It also provides functions to match donors across pools or 'omics by genotype.
- *Relevance:* the two-step mixed-prior design is the natural slot for plugging in imputed genotypes. The vireo docs recommend imputation only for array or WES genotypes of real donors, not for GT_donors ([docs](https://vireosnp.readthedocs.io/en/latest/genotype.html)) [V].

**Kang HM, Subramaniam M, Targ S, … Ye CJ (2018). Multiplexed droplet single-cell RNA-sequencing using natural genetic variation (demuxlet). *Nature Biotechnology* 36:89–94. doi:10.1038/nbt.4042** [V-abs]
- The original genotype-based method. It used imputed BeadChip genotypes (seen in search summaries; the full text was not open to me). freemuxlet is in the popscle suite and has no standalone paper.
- *Relevance:* this is the target method for the re-demultiplexing step.

**Weerakoon M, Vu H, Behboudi R, Heaton H (2026). Souporcell3: robust demultiplexing for high-donor single-cell RNA-seq datasets. *Bioinformatics* 42(3):btag117. doi:10.1093/bioinformatics/btag117** (bioRxiv 10.1101/2025.07.10.664218) [V]
- Scales souporcell to 64 donors using K-harmonic means and re-initializes low-quality clusters. Evaluated by ARI against demuxlet and ground-truth labels, with no genotype-accuracy analysis and no imputation.
- *Relevance:* the high-donor regime is where latent genotypes get sparsest per donor, so imputation should help most there.

**Rogozhnikov A, Ramkumar P, Shah K, Bedi R, Kato S, Escola GS (2021). Demuxalot: scaled up genetic demultiplexing for single-cell sequencing. bioRxiv 10.1101/2021.05.22.443646** [V-abs]
- Genotype-based, with a `learn_genotypes` step that refines and extends known (array) genotypes from the pooled reads. The motivation is that each genotype has about 10× more SNVs than the array captures.
- *Note:* one web summary claimed it was published in *Nat Commun* 12:5692. It is **not in PubMed**, so treat that claim as false; it is a bioRxiv preprint.
- *Relevance:* the closest existing "refine genotypes, then re-demultiplex" loop, but it works in the other direction (from known genotypes outward, using no LD or reference panel). In Demuxafy, "Demuxalot (refined)" was the best single method. That is a threat to our novelty framing and also a baseline to beat.

**Ranjbaran A, Luca F, Pique-Regi R (2026). fastdemux: robust SNP-based demultiplexing of single-cell population genomics data. bioRxiv 10.64898/2026.02.10.705082** [V]
- A DLDA-based genotype demultiplexer that uses DS, then GP, then GT from the VCF. Missing genotypes are set to the mean dosage. It does not evaluate how genotype uncertainty affects performance and does not take genotype-free output as input.
- *Relevance:* it natively consumes imputed **dosages**, so it would be a good downstream consumer of imputed latent genotypes.

**Schaefer NK, Pavlovic BJ, Schmitz MT, Pollen AA (2026). CellBouncer, a unified toolkit for single-cell demultiplexing and ambient RNA analysis. *Cell Genomics*. doi:10.1016/j.xgen.2026.101275** [V]
- A demultiplexing and ambient toolkit (the great-ape use case). Its full text has nothing on imputing inferred genotypes.

**Curion F, Wu X, Heumos L, … Theis FJ (2024). hadge: a comprehensive pipeline for donor deconvolution in single-cell studies. *Genome Biology*. doi:10.1186/s13059-024-03249-z** [V]
- Runs hashing and genotype methods in parallel, matches donors with a Phi score, and can pull donor-informative variants from **vireo-estimated donor genotypes** (sites with read depth > 10). Imputation (TOPMed-r2, Minimac4) appears only for the real array genotypes.
- *Relevance:* reuses vireo latent genotypes downstream but does not impute them.

**Zoodsma M, Zhan Q, Kumar S, … Li Y (2024). CellDemux: coherent genetic demultiplexing in single-cell and single-nuclei experiments. bioRxiv 10.1101/2024.01.18.576186** [V-abs]
- A consensus framework over 187 libraries and 800 samples (RNA, ATAC, multiome) that assigns 88% of donors. Methods not read.

**Hartoularos GC, Si Y, Zhang F, … Kang HM, Ye CJ (2023). Reference-free multiplexed single-cell sequencing identifies genetic modifiers of the human immune response ("clue"). bioRxiv 10.1101/2023.05.29.542756** [V-abs]
- Uses freemuxlet to encode 64 individuals across pools without reference genotypes or hashing, then "integrat[es] genotyping data" for response eQTLs.
- *Unverified detail:* whether freemuxlet genotypes were matched to or imputed from array data. Full text blocked. **Worth reading; it is a plausible threat.**

## 2. Tools that match latent clusters to external genotypes (common practice, never with imputation of the latent side)

- **Demuxafy docs, `Assign_Indiv_by_Geno.R`** [V]: Pearson correlation between souporcell cluster genotypes (GT/DS/GP) and reference genotypes, with a heatmap and key file. No imputation. https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Souporcell.html
- **souporcell GitHub issue #77 (2020)** [V]:
  - A user merged two samples genotyped on a UK Biobank Axiom array. souporcell used about 3,000 array SNPs directly, about 60,000 after Michigan imputation, and about 170,000 with `--common_variants`.
  - Imputed genotypes and common variants gave "much better results".
  - Heaton points to `shared_samples.py` for matching clusters to genotypes or across experiments.
  - *Relevance:* anecdotal evidence that imputation density matters for demultiplexing, but this is imputation of the **true** genotypes. https://github.com/wheaton5/souporcell/issues/77
- **souporcell GitHub issue #145 (2022)** [V]: a user merged `cluster_genotypes.vcf` with TOPMed-imputed array genotypes via bcftools and plink and got ambiguous cluster-to-donor matches. A commenter with 56 donors attributes some of this to sample mislabelling. https://github.com/wheaton5/souporcell/issues/145
- **souporcell issues #181 and #155** [V]: failures when passing imputed VCFs with `./.` entries or all-ref sites as `--known_genotypes`.
- **Demuxafy docs issue #17** [V]: demuxlet `--geno-error-coeff` needs an R2 INFO field from imputed VCFs. This is evidence that users already feed imputation-quality metrics into demuxlet.

## 3. Genotype calling / imputation from single-cell or RNA reads (adjacent)

**Dou J, Tan Y, Kock KH, … Chen K (2024). Single-nucleotide variant calling in single-cell sequencing data with Monopogen. *Nature Biotechnology* 42:803–812. doi:10.1038/s41587-023-01873-x** [V]
- Pools reads across cells of **one sample**, computes GLs, and refines them with a Beagle HMM on the 1KG3 panel. It gets more than 95% genotype accuracy against WGS for panel SNVs, and "the high accuracy is largely due to the LD-based genotyping refinement". It treats heterozygous calls that LD imputes to homozygous as sequencing errors. Works on scRNA, snRNA, snATAC and scDNA.
- *Relevance:* the closest methodological analogue. Running Monopogen on per-cluster BAMs from souporcell or vireo would be one concrete way to implement our idea. Europe PMC found no paper that combines Monopogen with souporcell, vireo or demultiplexing (3 hits, none relevant). Monopogen GitHub issue #34 is about missing genotypes, not pools.

**Buralkin I, Chen H, Park J, Liu Z (2026). scDeepVariant: a population-informed deep learning framework for germline variant calling in scRNA-seq. bioRxiv 10.64898/2025.12.31.696877** [V-abs]
- A DeepVariant model trained on paired WGS and snRNA, with gnomAD/1KG allele-frequency channels. It beats Monopogen above 10× depth and on rare variants, where LD refinement is weakest.
- *Relevance:* an alternative to LD imputation for "fixing" latent calls. It warns that LD refinement is limited for rare variants.

**Mu Z, Randolph HE, Aguirre-Gamboa R, … Barreiro LB (2026). Impact of disease-associated chromatin accessibility QTLs across immune cell types and contexts. *Cell Genomics*. doi:10.1016/j.xgen.2025.101061** [V]
- Per individual: bcftools GLs from **aggregated scATAC reads** at 1KG sites (DP ≤ 15), then joint GLIMPSE imputation across three studies, then Eagle phasing. They validated against array+Minimac4 in 13 Benaglio et al. individuals (INFO by MAF bin, high concordance) and used the result for caQTL mapping in 48 individuals.
- *Relevance:* the strongest evidence that GLIMPSE on pseudo-bulk single-cell reads gives QTL-grade genotypes. The difference from our idea is that the samples were one donor per library, not latent clusters from a pool.

**Zhang Q, Bao Q, Zeng L, He Z, Wang Z, Li C, Yi G (2026). A scalable framework for single-cell eQTL mapping … in pigs. *J Anim Sci Biotechnol* 17:132. doi:10.1186/s40104-026-01452-5** [V]
- Benchmarked SNP calling and imputation from 42 sc/snRNA-seq datasets. They chose cellSNP-lite, then SHAPEIT2 + Beagle 5.4, against a WGS truth set: CR > 0.90 and R² 0.80–0.90, with DR² > 0.95 and MAF ≥ 0.05 filters. After imputation, the variant distribution shifts from 3′UTR/intronic toward the genomic background.
- *Relevance:* a direct precedent for the "sc reads → cellSNP-lite → imputation" pipeline and its accuracy numbers. Again, not pooled.

**Guo K, Zhong Z, Zeng H, … Zhang Z (2025). Comparative analysis of genotype imputation strategies for SNPs calling from RNA-seq. *BMC Genomics* 26:245. doi:10.1186/s12864-025-11411-5** [V]
- Masks WGS down to RNA-SNP sites using 6,567 pig RNA-seq samples and compares Beagle, Minimac4 and IMPUTE5: CR 0.906–0.917 and r² 0.78–0.79, with no tool clearly better.
- *Relevance:* a site-pattern-only simulation (no sequencing error). It gives an upper bound for imputation from an RNA-like SNP footprint.

**Deelen P, Zhernakova DV, de Haan M, … Franke L (2015). Calling genotypes from public RNA-sequencing data enables identification of genetic variants that affect gene-expression levels. *Genome Medicine* 7:30. doi:10.1186/s13073-015-0152-4** [V-abs]
- Classic bulk RNA-seq genotype calling followed by imputation for eQTL mapping. Cited as the bulk precedent (my recollection is that they used Beagle; that detail was not re-checked).

**Quinones-Valdez G, Fu T, Chan TW, Xiao X (2022). scAllele. *Science Advances* 8:eabn6398. doi:10.1126/sciadv.abn6398** [V]: scRNA SNV and indel calling plus allele-linked splicing. No imputation.
**Huang X, Huang Y (2021). Cellsnp-lite: an efficient tool for genotyping single cells. *Bioinformatics* 37:4569–4571. doi:10.1093/bioinformatics/btab358** [V]: the pileup engine used before vireo.
**Wiens M, Farahani H, Scott RW, Underhill TM, Bashashati A (2024). Benchmarking bulk and single-cell variant-calling approaches on Chromium scRNA-seq and scATAC-seq libraries. *Genome Research* 34:1196–1210. doi:10.1101/gr.277066.122** [V]: bulk callers on pooled reads beat per-cell callers, measured against matched WGS. No imputation.
**Yao J, Gazal S (2026). Evaluating genetic ancestry inference from single-cell transcriptomic datasets. *HGG Advances* 7:100564. doi:10.1016/j.xhgg.2026.100564** [V]: ADMIXTURE on read-derived SNPs is accurate, applied to 401 HCA donors. Full text has no imputation.
**Hong SC, Muyas F, Cortés-Ciriano I, Hormoz S (2025). scAI-SNP. *BMC Methods* 2:10. doi:10.1186/s44330-025-00029-4** [V-abs]: ancestry inference from sc reads.
**Tomofuji Y, Edahiro R, … Okada Y (2024). Quantification of escape from X chromosome inactivation with single-cell omics data (scLinaX). *Cell Genomics*. doi:10.1016/j.xgen.2024.100625** [V]
- Uses imputed array genotypes where available (Michigan, Minimac4, relaxed R² 0.3 because "the genotype could be also confirmed by the allele information of the scRNA-seq reads"). Otherwise it uses genotypes called by cellSNP-lite from scRNA.
- *Relevance:* informal combining of read-level and imputed evidence (also relevant to the XCI project).

## 4. Low-pass WGS + imputation as a cheap genotype source

**Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O (2021). Efficient phasing and imputation of low-coverage sequencing data using large reference panels (GLIMPSE). *Nature Genetics* 53:120–126. doi:10.1038/s41588-020-00756-0** [V]
**Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O (2023). Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes (GLIMPSE2). *Nature Genetics* 55:1088–1090. doi:10.1038/s41588-023-01438-3** [V]
**Li JH, Mazur CA, Berisa T, Pickrell JK (2021). Low-pass sequencing increases the power of GWAS and decreases measurement error of polygenic risk scores compared to genotyping arrays. *Genome Research* 31:529–537. doi:10.1101/gr.266486.120** [V]
**Davies RW, Kucka M, Su D, … Myers S (2021). Rapid genotype imputation from sequence with reference panels (QUILT). *Nature Genetics* 53:1104–1111. doi:10.1038/s41588-021-00877-0** [V]; STITCH (panel-free): Davies et al. 2016, *Nat Genet*, doi:10.1038/ng.3594 [V]
**Kumar KH, Rubinacci S, Zöllner S (2026). MetaGLIMPSE: meta-imputation of low-coverage sequencing data. *AJHG*. doi:10.1016/j.ajhg.2026.02.004** [V-abs]
- *Relevance of this group:* a latent cluster from a pool is effectively a "low-pass, heavily biased (exonic/peak), ambient-contaminated" genome. GLIMPSE and QUILT are the right engines because they take GLs rather than hard calls.
- Nobody I found uses low-pass WGS **for demultiplexing**. The one sc paper that used low-pass WGS (Schott et al. 2022 *Cell Genomics*, doi:10.1016/j.xgen.2022.100207 [V]: 4× BGI plus Gencove imputation) used it for GWAS of challenge-study volunteers, while demultiplexing used 1000G LCL genotypes.

## 5. Re-identification / privacy / forensic matching of single-cell clusters to external genotypes

**Kulhankova L, Montiel González D, Bindels E, Kling D, Kayser M, Mulugeta E (2023). Single-cell transcriptome sequencing allows genetic separation, characterization and identification of individuals in multi-person biological mixtures ("De-goulash"). *Communications Biology* 6:201. doi:10.1038/s42003-023-04557-z** [V]
- Mixed blood from 2–5 people, separated by souporcell-style SNP clustering. They then compare the per-cluster SNP profiles against a reference genotype database to assign cluster identity, and report high match rates with as few as tens of cells per cluster.
- *Relevance:* the only paper I found that explicitly raises imputation of cluster-level genotypes, in its discussion, as a way to densify sparse per-cluster SNP sets. It is stated as an idea for future work, not implemented. This is the clearest "someone has thought of this" precedent for part (b).

**Walker CR, Li X, Chakravarthy M, Lounsbery-Scaife W, … Gürsoy G (2024). Private information leakage from single-cell count matrices. *Cell* 187:6537. doi:10.1016/j.cell.2024.09.012** [V-abs]
- Builds pseudo-bulk profiles per cell type and predicts genotypes through eQTL associations, reporting high linking accuracy to genotype panels (>80% with GTEx eQTLs; >70% with one cell type). Commentary: **Cho H (2024). Privacy of single-cell gene expression data. *Patterns* 5:101096. doi:10.1016/j.patter.2024.101096** [V], which notes that "accurate reconstruction of genotypes at individual genomic positions is not necessary for data linkage".
- *Relevance:* genotype prediction here comes from **expression**, not reads, and no imputation is used. It shows that cluster-to-panel matching is tractable even with poor per-site accuracy, which is encouraging for cluster-identity recovery but does not address the accuracy question we care about.

**Gürsoy G, Emani P, … Gerstein M (2020). Data sanitization to reduce private information leakage from functional genomics. *Cell* 183:905–917. doi:10.1016/j.cell.2020.09.036** [V-abs] and **Harmanci A, Gerstein M (2016). Quantification of private information leakage from phenotype-genotype data: linking attacks. *Nature Methods* 13:251–256. doi:10.1038/nmeth.3746** [V-abs]
- Foundational quantification of variant leakage from functional-genomics reads, including 10x scRNA.
- *Relevance:* background for why per-cluster genotype quality matters ethically, and the framework for scoring cluster-to-individual matches.

## 6. Benchmarks that scale donor number (and whether they measured latent genotype quality)

**Neavin D, Senabouth A, Arora H, … Powell JE (2024). Demuxafy. *Genome Biology* 25:94. doi:10.1186/s13059-024-03224-8** [V]
- Benchmarks 7 demultiplexers and several doublet detectors on pools from 2 to at least 16 individuals; 0.4–78.8% of droplets were misclassified depending on method and pool size. **Demuxalot (refined)** plus DoubletDetection had the fewest errors (about 1% in 2-donor pools, about 3% at 16+). The reference genotypes were arrays imputed on the Michigan server (Minimac3, HRC r1.1 and 1000G p3).
- *Relevance:* the field-standard benchmark. It scores droplets only, never genotypes. Its "refined" winner shows that improving genotypes improves assignment, which is our premise, with no reference panel involved.

**Li T, Alvarez M, Liu C, … Zaitlen N (2025). The impact of ambient contamination on demultiplexing methods for single-nucleus multiome experiments. *eLife* (reviewed preprint) 106769; bioRxiv 10.1101/2025.02.06.636969. doi:10.7554/eLife.106769** [V]
- Introduces **ambisim**, a genotype-aware read-level joint snRNA/snATAC simulator, and varies ambient fraction, doublet rate, donor number and coverage. Key line for us: for freemuxlet and scSplit, ambient fraction did **not** separate accurate from inaccurate droplets, "suggesting that latent genotype inference is the driving force". Genotype-based methods beat genotype-free ones on average; ATAC genotype-free methods are the most coverage-sensitive. It also introduces a "variant consistency" metric. Reference VCFs were TOPMed-imputed, filtered at R² > 0.90.
- *Relevance:* the most direct statement in the literature that latent genotype quality is the limiting factor, and the paper stops exactly short of measuring or fixing it. ambisim is also the natural simulator for our evaluation, since ground-truth genotypes are known by construction.

**Alvarez M, Li T, Lee SHT, … (2025). Integrated ambient modeling and genetic demultiplexing of single-cell RNA+ATAC multiome experiments with Ambimux. bioRxiv 10.1101/2025.08.21.671671; PMC12407767** [V]
- A genotype-based multiome demultiplexer that jointly models ambient molecules. Genotypes were Eagle v2.4 phased and **Minimac4 imputed**, filtered at R² > 0.99. Limitation stated by the authors: ambimux "was designed for multiplexed experiments in which all donors have genotype data available. In practice, it may be challenging to obtain this for every donor", and "an extension of a genotype-free or missing-genotype framework would further expand the utility of ambimux".
- *Relevance:* this is the stated gap our proposal fills, and ambimux is the obvious re-demultiplexing engine for the second pass.

**Fu et al. (2025). Benchmarking of computational demultiplexing methods for single-cell RNA sequencing. *Briefings in Bioinformatics* 26(4):bbaf371. doi:10.1093/bib/bbaf371** [V]
- Vireo, souporcell, freemuxlet and scSplit on simulated pools of 2, 4 and 6 from six donors, with either SNP-array or matched bulk RNA-seq variants as reference. It compares variant callers (bcftools, cellSNP, FreeBayes) on the bulk RNA, excludes scSplit for low accuracy, and validates with sex-linked genes. **No imputation and no latent-genotype evaluation.**
- *Relevance:* it shows bulk-RNA-derived genotypes being used as a reference for demultiplexing, one step short of imputing them.

**Weber LM, Hippen AA, Hickey PF, … Hicks SC (2021). Genetic demultiplexing of pooled single-cell RNA-sequencing samples in cancer facilitates effective experimental design. *GigaScience* 10:giab062. doi:10.1093/gigascience/giab062** [V-abs]
- In-silico pooling in ovarian and lung cancer; genetic demultiplexing works despite somatic variation, but high ambient RNA hurts. Droplet-level metrics only.

**Nassiri I, Kwok AJ, … (2024). Demultiplexing of single-cell RNA-sequencing data using interindividual variation in gene expression. *Bioinformatics Advances* 4:vbae085. doi:10.1093/bioadv/vbae085** [V-abs]; **Tang Z, et al. (2024). MitoSort. *Genomics Proteomics Bioinformatics* 22:qzae073. doi:10.1093/gpbjnl/qzae073** [V-abs]; **Sen E, et al. (2026). Disagreement between demultiplexing methods reveals structured cell quality gradients (Split-flow). bioRxiv 10.64898/2026.05.10.724135** [V-abs]; **Mears J, Orchard P, Varshney A, … (2026). Genetic demultiplexing and TSS identification from nanopore sequencing of 10x multiome libraries. bioRxiv 10.64898/2026.03.31.715454** [V-abs]
- Orthogonal or alternative demultiplexing signals (expression, mtDNA, hashing concordance, long reads). None involve reference-panel imputation of inferred genotypes.

**Neavin DR, Steinmann AM, Farbehi N, … Powell JE (2023). A village in a dish model system for population-scale hiPSC studies. *Nature Communications* 14:3240. doi:10.1038/s41467-023-38704-1** [V-abs]
- Pools many hiPSC lines in one culture ("village"), which is the design regime where genotype-free scaling matters most.

## 7. Imputation as genotype error correction, in general

- **Monopogen** (above) is the clearest published example of using an LD reference panel to *correct* per-site calls, and it explicitly reclassifies het calls that impute to hom as sequencing errors.
- **Browning BL, Zhou Y, Browning SR (2018). A one-penny imputed genome from next-generation reference panels. *AJHG* 103:338–348. doi:10.1016/j.ajhg.2018.07.015** [V] and **Browning BL, Tian X, Zhou Y, Browning SR (2021). Fast two-stage phasing of large-scale sequence data. *AJHG* 108:1880–1890. doi:10.1016/j.ajhg.2021.08.005** [V]: Beagle 5.x, the engine most likely to accept GLs from cluster pseudo-bulk.
- **Das S, Forer L, Schönherr S, … Fuchsberger C (2016). Next-generation genotype imputation service and methods (Minimac3/Michigan Imputation Server). *Nature Genetics* 48:1284–1287. doi:10.1038/ng.3656** [V]; TOPMed panel: **Taliun D, et al. (2021). *Nature* 590:290–299. doi:10.1038/s41586-021-03205-y** [V-abs].
- **Lau W, Ali A, Maude H, … (2024). The hazards of genotype imputation when mapping disease susceptibility variants. *Genome Biology* 25:17. doi:10.1186/s13059-023-03140-3** [V-abs]: a caution that imputation introduces its own errors, especially for rare variants and mismatched panels.
- **Furuta T, Yamamoto T, Ashikari M (2023). GBScleanR: robust genotyping error correction using an HMM with error pattern recognition. *Genetics* 224:iyad055. doi:10.1093/genetics/iyad055** [V-abs]: explicit error-correction-by-HMM framing, but in biparental plant populations, not human panels.
- **Day TRC, et al. (2025). Adjustment for genotype imputation uncertainty corrects for inflated type I error in family-based association testing. *Genetic Epidemiology*. doi:10.1002/gepi.70021** [V-abs]: a reminder to propagate imputation uncertainty (dosages or GP) rather than hard-calling, which matters if imputed latent genotypes feed a likelihood-based demultiplexer.

---

## Gaps / novelty

1. **No one has published latent-genotype accuracy as a primary result.** There is no curve of concordance or dosage R² against MAF, per-cluster depth, donor number, ambient fraction, or modality (RNA vs ATAC vs multiome). The three existing datapoints (scSplit ~98.4% at distinguishing sites; souporcell Fig. S1i; vireo >10 reads) are side analyses on small pools.
2. **No one has imputed latent genotypes.** The souporcell `cluster_genotypes.vcf` (GT/GL and ref/alt counts) and vireo `GT_donors.vireo.vcf` are exactly the GL-bearing inputs GLIMPSE, Beagle or Minimac4 want, and nothing in the indexed literature or in the tools' issue trackers does this.
3. **No one has closed the loop.** Genotype-free pass → per-cluster GLs → panel imputation → re-demultiplex with demuxlet, vireo-with-GT, ambimux or fastdemux, then check whether droplet assignment, doublet detection and unassigned rates improve. Demuxalot's `learn_genotypes` is the closest, and it refines from known genotypes without a panel.
4. **The ATAC/multiome case is wide open.** Latent genotypes in ATAC are sparser per site but spread over far more of the genome, which is the regime where LD imputation helps most, and Li et al. 2025 show ATAC genotype-free methods are the most coverage-sensitive. Mu et al. 2026 already proved GLIMPSE works on aggregated scATAC reads for single donors.
5. **A useful negative-control question nobody has answered:** how much of the apparent gain would come simply from denser common-variant sites rather than from genuine error correction? souporcell issue #77 hints that density alone matters a lot.
6. **Ambient contamination interacts with imputation in an untested way.** Ambient molecules create the specific hom-ref-to-het error that souporcell documents, and LD refinement is exactly what should undo it, but nobody has tested whether panel imputation removes or amplifies ambient-induced errors.

## Threats (closest prior work)

| Work | Why it is a threat | Why it does not scoop us |
|---|---|---|
| **Monopogen** (Dou 2024 *Nat Biotech*) | Same core mechanic: pooled single-cell reads → GLs → LD/panel refinement → >95% accuracy. Someone could apply it per cluster in a weekend. | One sample per library; never demultiplexing; no re-demultiplexing loop; no evaluation with pooled or ambient-contaminated data. |
| **Mu et al. 2026** (*Cell Genomics*, GLIMPSE on aggregated scATAC) | Proves the imputation step works on single-cell pseudo-bulk and validates against array+Minimac4. | Individually processed donors, not latent clusters; genotypes used for caQTL mapping, not fed back into demultiplexing. |
| **Demuxalot (refined)** (Rogozhnikov 2021; best method in Demuxafy) | An existing, well-benchmarked "refine genotypes then re-demultiplex" loop. | Requires known genotypes to start; refines from data only, with no reference panel, no LD and no imputation of unobserved sites; cannot start from a genotype-free pool. |
| **Kulhankova 2023** (*Commun Biol*) | Explicitly floats imputing cluster-level genotypes to densify them. | Stated as future work in the discussion, not implemented, and aimed at identity matching rather than demultiplexing accuracy. |
| **Hartoularos 2023** ("clue", freemuxlet at 64 donors) **[methods unverified]** | Large-scale reference-free design with genotype integration; could contain an unnoticed matching or refinement step. | Abstract describes integrating genotype data for eQTL mapping, not imputing freemuxlet genotypes. **Read the full text before writing anything up.** |
| **Ambimux** (Alvarez 2025) | Names the genotype-free/missing-genotype extension as its own next step, and the group overlaps with the Li 2025 benchmark. | Not implemented there. Also an argument that the idea is timely rather than that it is taken. |
| **Souporcell3** (2026) | Active development on high-donor pools by the original author; could add genotype refinement next. | Purely a clustering improvement; ARI-based evaluation; no genotype output analysis. |
| **scDeepVariant** (2026 preprint) | A learned alternative that might outcompete LD imputation for fixing sc variant calls. | Single-sample germline calling; no demultiplexing; needs paired WGS training data. |

### Verification caveats
- bioRxiv full-text HTML and PDF were rate-limited throughout, so **Hartoularos 2023, CellDemux, Souporcell3 preprint, Split-flow, Mears 2026 and demuxalot** are recorded from abstracts or secondary summaries.
- GitHub *code* search (as opposed to issue search) needs authentication and was not run, so a script in someone's pipeline repo that pipes `cluster_genotypes.vcf` into Beagle would not have been found.
- Claims marked [V-abs] should be re-checked before they go into a paper.
- One search summary claimed demuxalot was published in *Nature Communications* 12:5692. PubMed does not index it, so that citation appears to be wrong; cite the bioRxiv DOI.
