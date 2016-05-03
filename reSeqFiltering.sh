# Annotation
java -jar /bin/GenomeAnalysisTK.jar -R ~/Work/Genomes/ucsc.hg19.fasta -T VariantAnnotator -V ~/Work/Genomes/hapmap_3.3.hg19.vcf -o ~/Work/Genomes/hapmap_3.3.hg19_VariantAnnotator.vcf -nt 3 -all
java -jar /bin/GenomeAnalysisTK.jar -R ~/Work/Genomes/ucsc.hg19.fasta -T VariantAnnotator -V ~/Work/Genomes/1000G_omni2.5.hg19.vcf -o ~/Work/Genomes/1000G_omni2.5.hg19_VariantAnnotator.vcf -nt 3 -all
java -jar /bin/GenomeAnalysisTK.jar -R ~/Work/Genomes/ucsc.hg19.fasta -T VariantAnnotator -V ~/Work/Genomes/1000G_phase1.snps.high_confidence.hg19.vcf -o ~/Work/Genomes/1000G_phase1.snps.high_confidence.hg19_VariantAnnotator.vcf -nt 3 -all
java -jar /bin/GenomeAnalysisTK.jar -R ~/Work/Genomes/ucsc.hg19.fasta -T VariantAnnotator -V ~/Work/Genomes/dbsnp_138.hg19.vcf -o ~/Work/Genomes/dbsnp_138.hg19_VariantAnnotator.vcf -nt 3 -all
java -jar /bin/GenomeAnalysisTK.jar -R ~/Work/Genomes/ucsc.hg19.fasta -T VariantAnnotator -V round1_mpileup_multisampleUCSCannotated.vcf.gz -o round1_mpileup_multisampleUCSCannotated_VariantAnnotator.vcf -nt 3 -all

# Site-specific filtering
# 1. Variant Quality Score Recalibration (VQSR)
java -jar /bin/GenomeAnalysisTK.jar -T VariantRecalibrator -R /Users/robin/Work/Genomes/ucsc.hg19.fasta -input round1_mpileup_multisampleUCSCannotated_VariantAnnotator.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 Work/Genomes/hapmap_3.3.hg19_VariantAnnotator.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 Work/Genomes/1000G_omni2.5.hg19_VariantAnnotator.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 Work/Genomes/1000G_phase1.snps.high_confidence.hg19_VariantAnnotator.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Work/Genomes/dbsnp_138.hg19_VariantAnnotator.vcf -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -nt 3 -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R
java -jar GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R  /Users/robin/Work/Genomes/ucsc.hg19.fasta \
    -input round1_mpileup_multisample.vcf.gz \
    -mode SNP \
    --ts_filter_level 99.9 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -o round1_mpileup_multisample.recalibrated.vcf
# 2. Genotype-level checks: Exclude sites that have a) depth < 8 OR b) GQ < 20
# 3. Genotype call rate for sites must be >90%, if not site should be excluded
vcftools --max-missing 0.1 --gzvcf round1_mpileup_multisample.vcf.gz --minDP 8 --minGQ 20 --recode-INFO-all --recode --out step_a2_a3
# 4. Call rate between cases and controls per site must be similar - examine the call rate difference per site via a chi-squared test, exclude extreme values
! need caseControl.txt

# Sample-level QC
# 1. Ts/Tv ratio: For exomes this is normally between 2.7 and 3. Check in exons only, generate a distribution of these ratios, Ts/Tv vs frequency, log-normalise and remove samples 3 s.d. from the mean
# 2. Number of variants (SNPs and indels) across samples must be similar
# 3. Ancestry via PCA (Note: We don't have enough information in 63 loci to do this, but ask Leeds to confirm these are all European. I know from previous projects that there are some that are not)
! check with Leeds
# 4. Remove related and duplicate samples. Normally this is done by IBD but again we do not have enough information, so again ask Leeds to confirm.
! check with Leeds
