from itertools import product

METHODS = ["JTI", "PrediXcan"]
GWAS_TISSUE_PAIRS =  [
  ("30740_irnt", "Adipose_Visceral_Omentum"),    # glucose (quantile)
  ("30740_irnt", "Liver"),                       # glucose (quantile)
  ("30740_irnt", "Muscle_Skeletal"),             # glucose (quantile)
  ("30740_irnt", "Pancreas"),                    # glucose (quantile)
  ("30760_irnt", "Liver"),                       # HDL (quantile)
  ("30780_irnt", "Liver"),                       # LDL (quantile)
  ("30890_irnt", "Skin_Sun_Exposed_Lower_leg"),  # Vitamin D (quantile)
  ("30710_irnt", "Whole_Blood"),                 # C-reactive protein (quantile)
  ("30700_irnt", "Kidney_Cortex"),               # Creatinine (quantile)
]
LOOKUP = {
  # ref: https://nealelab.github.io/UKBB_ldsc/downloads.html
  "30740_irnt": "https://dropbox.com/s/d16kjpukj1qmg6e/30740_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  "30760_irnt": "https://dropbox.com/s/65jisgxwbbdrkaw/30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  "30780_irnt": "https://dropbox.com/s/4rnjzczwjgs5pgl/30780_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  "30890_irnt": "https://dropbox.com/s/bbbnj239w1kjkyi/30890_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  "30710_irnt": "https://dropbox.com/s/4o9qs1s6jkrzrq7/30710_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  "30700_irnt": "https://dropbox.com/s/n6xdh0149u0yoy0/30700_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
}

SPRED_OUT = [f"{g}-{m}_{t}" for (g, t), m in product(GWAS_TISSUE_PAIRS, METHODS)]

rule all:
  input: expand("results/{out}.csv", out=SPRED_OUT)

rule clean:
    shell: "rm -rf data results hakyimlab-MetaXcan-e211445"
  
rule install_MetaXcan:
  output: "hakyimlab-MetaXcan-e211445/software/SPrediXcan.py"
  shell:
    """
    wget https://github.com/hakyimlab/MetaXcan/zipball/e211445 -O tmp.zip
    unzip tmp.zip
    rm tmp.zip
    """

rule download_db:
  output: "data/{sample}.db"
  shell: "wget https://zenodo.org/record/3842289/files/{wildcards.sample}.db -O {output}"

rule download_cov:
  output: "data/{sample}.txt.gz"
  shell: "wget https://zenodo.org/record/3842289/files/{wildcards.sample}.txt.gz -O {output}"

rule download_and_merge_gwas:
  input: "data/GWAS/variants.sorted.tsv"
  output: "data/GWAS/{ukbb_id}.gwas.imputed_v3.both_sexes.tsv"
  params: lambda wildcards: LOOKUP[wildcards.ukbb_id]
  shell: "paste <(cut -f 2- {input}) <(wget -qO - {params} | gunzip) > {output}"

rule download_variants:
  output: "data/GWAS/variants.tsv.bgz"
  shell: "wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz -O {output}"

rule sort_variants:
  input: "data/GWAS/variants.tsv.bgz"
  output: "data/GWAS/variants.sorted.tsv"
  shell:
    """
    gunzip -c {input} | head -n 1 | cut -f -6 > {output} || true
    sort <(gunzip -c {input} | cut -f -6 | sed 1d) >> {output} || true
    """

rule run_SPrediXcan:
  input:
    covariance="data/{sample}.txt.gz",
    model_db="data/{sample}.db",
    gwas_file="data/GWAS/{ukbb_id}.gwas.imputed_v3.both_sexes.tsv",
    SPrediXcan="hakyimlab-MetaXcan-e211445/software/SPrediXcan.py"
  output: "results/{ukbb_id}-{sample}.csv"
  shell:
    """
    python {input.SPrediXcan} --model_db_path {input.model_db} \
      --covariance {input.covariance} \
      --gwas_file {input.gwas_file} \
      --snp_column rsid \
      --effect_allele_column alt \
      --non_effect_allele_column ref \
      --chromosome_column chr \
      --position_column pos \
      --beta_column beta \
      --se_column se \
      --pvalue_column pval \
      --freq_column minor_AF \
      --output_file {output} \
    """

rule render_rmd:
  input: "results/30780_irnt-PrediXcan_Liver.csv", "results/30780_irnt-JTI_Liver.csv"
  output: "plots.html"
  shell: 
    """
    R -e "rmarkdown::render('plots.Rmd')"
    """