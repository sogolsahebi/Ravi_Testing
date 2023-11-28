from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

# Configuration
configfile: "config.yaml"

# Define S3 Remote Provider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]

#TODO: change the github link later
data_source = "https://github.com/sogolsahebi/Ravi_Testing"

# Rule to download data
rule download_data:
    output:
        S3.remote(prefix + "download/Table_S1_Clinical_Annotations.csv"),
        S3.remote(prefix + "download/SU2C-MARK_Harmonized_Clinical_Annotations_Supplement_v1.txt"),
        S3.remote(prefix + "download/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct"),
        S3.remote(prefix + "download/Source_Data.zip")
    resources:
        mem_mb=2000
    shell:
        """
        wget -O {output[3]} {data_source}/files/Source_Data.zip
        unzip {output[3]} -d {prefix}download/
        cp {prefix}download/Clinical/Table_S1_Clinical_Annotations.csv {output[0]}
        cp {prefix}download/Clinical/SU2C-MARK_Harmonized_Clinical_Annotations_Supplement_v1.txt {output[1]}
        cp {prefix}download/RNA/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct {output[2]}
        """

# Rule to format downloaded data
rule format_download_data:
    input:
        S3.remote(prefix + "download/Table_S1_Clinical_Annotations.csv"),
        S3.remote(prefix + "download/SU2C-MARK_Harmonized_Clinical_Annotations_Supplement_v1.txt"),
        S3.remote(prefix + "download/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct")
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt.gz"),
        S3.remote(prefix + "download/SNV.txt.gz")
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \\
        {input[0]} \\
        {input[1]} \\
        {input[2]} \\
        {output[0]} \\
        {output[1]} \\
        {output[2]}
        """

# Rule to format clinical data
rule format_clin:
    input:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    output:
        S3.remote(prefix + "processed/CLIN.csv")
    shell:
        """
        Rscript scripts/Format_CLIN.R \\
        {input[1]} \\
        {input[2]} \\
        {input[3]} \\
        {output}
        """

# Rule to format expression data
rule format_expr:
    input:
        S3.remote(prefix + "download/EXPR.txt.gz"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    output:
        S3.remote(prefix + "processed/EXPR.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_EXPR.R \\
        {input[0]} \\
        {input[1]} \\
        {output}
        """

# Rule to format cased sequenced data
rule format_cased_sequenced:
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt.gz"),
        S3.remote(prefix + "download/SNV.txt.gz")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv")
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \\
        {input[0]} \\
        {input[1]} \\
        {input[2]} \\
        {output}
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v19.annotation.RData?raw=true -O {output[0]}
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_drug.csv -O {output[1]}
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_tissue.csv -O {output[2]}
        """

# Rule to get MultiAssayExp
rule get_MultiAssayExp:
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData")
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb=3000,
        disk_mb=3000
    shell:
        """
        Rscript -e " \\
        'load(\"{input[4]}\"); \\
        source(\"https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R\"); \\
        saveRDS( \\
            get_MultiAssayExp(study = \"Ravi\", input_dir = \"{prefix}processed\"), \\
            \"{output}\" \\
        );' \\
        """
