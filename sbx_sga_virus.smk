ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
try:
    SBX_SGA_VERSION = get_ext_version("sbx_sga")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_SGA_VERSION = "0.0.0"


localrules:
    all_sga_virus,


rule all_sga_virus:
    input:
        expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_virus_summary.tsv",
            sample=Samples,
        ),


rule sga_genomad_download_db:
    """Download Genomad database"""
    output:
        version=Path(Cfg["sbx_sga"]["genomad_ref"]) / "genomad_db" / "version.txt",
    log:
        LOG_FP / "sga_genomad_download_db.log",
    benchmark:
        BENCHMARK_FP / "sga_genomad_download_db.tsv"
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        GENOMAD_DB_DIR=$(dirname {output.version})
        genomad download-database $(dirname "$GENOMAD_DB_DIR") > {log} 2>&1 || true
        """


rule sga_genomad_end_to_end:
    """Run Genomad end-to-end pipeline"""
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
        db_version=Path(Cfg["sbx_sga"]["genomad_ref"]) / "genomad_db" / "version.txt",
    output:
        assembly_summary=ISOLATE_FP
        / "genomad"
        / "{sample}"
        / "{sample}_summary"
        / "{sample}_virus_summary.tsv",
    log:
        LOG_FP / "genomad_end_to_end_{sample}.log",
    benchmark:
        BENCHMARK_FP / "genomad_end_to_end_{sample}.tsv"
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        ASSEMBLY_SUMMARY_DIR=$(dirname {output.assembly_summary})
        DB_DIR=$(dirname {input.db_version})
        
        if [ ! -s {input.contigs} ]; then
            touch {output.assembly_summary}
        else
            genomad end-to-end --cleanup --splits 8 {input.contigs} $(dirname "$ASSEMBLY_SUMMARY_DIR") $DB_DIR > {log} 2>&1
        fi
        """
