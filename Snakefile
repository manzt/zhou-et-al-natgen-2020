from itertools import product

from snakemake.remote.HTTP import RemoteProvider

METHODS = ["JTI", "PrediXcan"]
GWAS_TISSUE_PAIRS =  [
  ("30740_irnt", "Adipose_Visceral_Omentum"),
  ("30740_irnt", "Liver"),
  ("30740_irnt", "Muscle_Skeletal"),
  ("30740_irnt", "Pancreas"),
  ("30760_irnt", "Liver"),
  ("30780_irnt", "Liver"),
  ("30890_irnt", "Skin_Sun_Exposed_Lower_leg"),
  ("30710_irnt", "Whole_Blood"),
  ("30700_irnt", "Kidney_Cortex"),
]
LOOKUP = { # ref: https://nealelab.github.io/UKBB_ldsc/downloads.html
  "30740_irnt": ("https://dropbox.com/s/d16kjpukj1qmg6e/30740_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "https://dropbox.com/s/9zuen4u9ikity1q/30740_irnt.imputed_v3.ldsc.both_sexes.tsv.gz"),
  "30760_irnt": ("https://dropbox.com/s/65jisgxwbbdrkaw/30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "https://dropbox.com/s/qu1qp9pvket42m6/30760_irnt.imputed_v3.ldsc.both_sexes.tsv.gz"),
  "30780_irnt": ("https://dropbox.com/s/4rnjzczwjgs5pgl/30780_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "https://dropbox.com/s/fk0th6tarp8q3e5/30780_irnt.imputed_v3.ldsc.both_sexes.tsv.gz"),
  "30890_irnt": ("https://dropbox.com/s/bbbnj239w1kjkyi/30890_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "https://dropbox.com/s/wwbxxjnm3k9wiv0/30890_irnt.imputed_v3.ldsc.both_sexes.tsv.gz"),
  "30710_irnt": ("https://dropbox.com/s/4o9qs1s6jkrzrq7/30710_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "https://dropbox.com/s/tie2ft1fie5p9ey/30710_irnt.imputed_v3.ldsc.both_sexes.tsv.gz"),
  "30700_irnt": ("https://dropbox.com/s/n6xdh0149u0yoy0/30700_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "https://dropbox.com/s/rlyuzamxja4j6lg/30700_irnt.imputed_v3.ldsc.both_sexes.tsv.gz"),
}

SPRED_OUT = [f"{g}-{m}_{t}" for (g, t), m in product(GWAS_TISSUE_PAIRS, METHODS)]
http = RemoteProvider()

rule all:
  input: expand("results/{out}.csv", out=SPRED_OUT)

rule clean:
    shell: "rm -rf data results hakyimlab-MetaXcan-e211445"
  
rule install_metaXcan:
  input: http.remote("https://github.com/hakyimlab/MetaXcan/zipball/e211445")
  output: "hakyimlab-MetaXcan-e211445/software/SPrediXcan.py"
  shell: "unzip {input} && rm {input}"

rule download_db:
  input: http.remote("https://zenodo.org/record/3842289/files/{sample}.db", keep_local=True)
  output: "data/{sample}.db"
  shell: "mv {input} {output}"

rule download_cov:
  input: http.remote("https://zenodo.org/record/3842289/files/{sample}.txt.gz", keep_local=True)
  output: "data/{sample}.txt.gz"
  shell: "mv {input} {output}"

rule download_gwas:
  output: "data/GWAS/{ukbb_id}.gwas.imputed_v3.both_sexes.tsv.bgz"
  params: lambda wildcards: LOOKUP[wildcards.ukbb_id][0]
  shell: "wget {params} -O {output}"

rule download_ldsc:
  output: "data/GWAS/{ukbb_id}.imputed_v3.ldsc.both_sexes.tsv.gz"
  params: lambda wildcards: LOOKUP[wildcards.ukbb_id][1]
  shell: "wget {params} -O {output}"

rule run_SPrediXcan:
  input:
    covariance="data/{sample}.txt.gz",
    model_db="data/{sample}.db",
    gwas_file="data/GWAS/{ukbb_id}.imputed_v3.ldsc.both_sexes.tsv.gz",
    SPrediXcan="hakyimlab-MetaXcan-e211445/software/SPrediXcan.py"
  output: "results/{ukbb_id}-{sample}.csv"
  shell:
    """
    python {input.SPrediXcan} --model_db_path {input.model_db} \
      --covariance {input.covariance} \
      --gwas_file {input.gwas_file} \
      --zscore_column Z \
      --output_file {output} \
    """