def get_sga_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_sga":
            return Path(fp)
    raise Error(
        "Filepath for sbx_sga not found, are you sure it's installed under extensions/sbx_sga?"
    )


ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
SBX_SBA_VERSION = open(get_sga_path() / "VERSION").read().strip()
TOOLS = ["shovill", "mlst", "checkm", "amr", "bakta", "mash"]


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


localrules:
    all_sga,


rule all_sga:
    input:
        # QC
        expand(ISOLATE_FP / "mash" / "{sample}_sorted_winning.tab", sample=Samples),
        # Assembly QC
        expand(
            ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv", sample=Samples
        ),
        expand(ISOLATE_FP / "quast" / "{sample}" / "report.tsv", sample=Samples),
        # Typing
        expand(ISOLATE_FP / "mlst" / "{sample}.mlst", sample=Samples),
        # Annotation
        expand(ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt", sample=Samples),
        # AMR Profiling
        expand(ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out", sample=Samples),
        f"{ISOLATE_FP}/reports/shovill.report",
        f"{ISOLATE_FP}/reports/mlst.report",
        f"{ISOLATE_FP}/reports/checkm.report",
        f"{ISOLATE_FP}/reports/amr.report",
        f"{ISOLATE_FP}/reports/bakta.report",
        f"{ISOLATE_FP}/reports/mash.report",
        f"{ISOLATE_FP}/final_summary.tsv",


rule sga_mash:
    input:
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        agg=temp(ISOLATE_FP / "mash" / "{sample}.fastq"),
        win=temp(ISOLATE_FP / "mash" / "{sample}_winning.tab"),
        sort=ISOLATE_FP / "mash" / "{sample}_sorted_winning.tab",
    params:
        ref=Cfg["sbx_sga"]["mash_ref"],
    log:
        LOG_FP / "sga_mash_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mash_{sample}.tsv"
    conda:
        "envs/mash.yml"
    shell:
        """
        zcat  {input.reads} > {output.agg}
        mash screen -w -p 8 {params.ref} {output.agg} > {output.win} 2> {log}
        sort -gr {output.win} > {output.sort} 2>> {log}
    """


rule sga_shovill:
    input:
        rp1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        rp2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    log:
        LOG_FP / "sga_shovill_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_shovill_{sample}.tsv"
    conda:
        "envs/shovill.yml"
    shell:
        """
        (shovill --force --assembler skesa --outdir $(dirname {output.contigs}) --R1 {input.rp1}  --R2 {input.rp2} &> {log} &&
        mv $(dirname {output.contigs})/contigs.fa {output.contigs}) ||
        touch {output.contigs}
        """


### Assembly QC
rule sga_checkm:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        quality_report=ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv",
    params:
        ref=Cfg["sbx_sga"]["checkm_ref"],
    log:
        LOG_FP / "sga_checkm_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_checkm_{sample}.tsv"
    conda:
        "envs/checkm2.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            checkm2 predict \\
            --force \\
            -x fa \\
            -i {input.contigs} \\
            -o $(dirname {output.quality_report}) \\
            --database_path {params.ref} \\
            &> {log} || touch {output.quality_report}
        else
            touch {output.quality_report}
        fi
        """


rule sga_quast:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        quast_dir=ISOLATE_FP / "quast" / "{sample}" / "report.tsv",
    log:
        LOG_FP / "sga_quast_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_quast_{sample}.tsv"
    conda:
        "envs/quast.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            quast.py \\
                -o $(dirname {output.quast_dir}) \\
                {input.contigs} \\
                &> {log} || touch {output.quast_dir}
        else
            touch {output.quast_dir}
        fi
        """


### Typing
rule sga_mlst:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        mlst=ISOLATE_FP / "mlst" / "{sample}.mlst",
    log:
        LOG_FP / "sga_mlst_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mlst_{sample}.tsv"
    conda:
        "envs/mlst.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            mlst --nopath {input.contigs} > {output.mlst} 2> {log}
        else
            touch {output.mlst}
        fi
        """


### Annotation
rule sga_bakta:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        bakta=ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt",
    params:
        ref=Cfg["sbx_sga"]["bakta_ref"],
    log:
        LOG_FP / "sga_bakta_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_bakta_{sample}.tsv"
    conda:
        "envs/bakta.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            bakta --db {params.ref} \\
            --output $(dirname {output.bakta}) \\
            --prefix {wildcards.sample} \\
            --skip-plot {input.contigs} \\
            &> {log}
        else
            touch {output.bakta}
        fi
        """


### AMR Profiling
rule sga_abritamr:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        abritamr=ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out",
    log:
        LOG_FP / "sga_abritamr_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_abritamr_{sample}.tsv"
    conda:
        "envs/abritamr.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            abritamr run \\
                --contigs {input.contigs} \\
                --prefix {wildcards.sample} \\
                &> {log}
            mv {wildcards.sample} $(dirname $(dirname {output.abritamr}))
        else
            mkdir -p $(dirname {output.abritamr})
            touch {output.abritamr}
        fi 
    """


### MLST Report
rule mlst_summary:
    input:
        reports=expand(ISOLATE_FP / "mlst" / "{sample}.mlst", sample=Samples),
    output:
        mlst_report=ISOLATE_FP / "reports" / "mlst.report",
    shell:
        """
        for f in {input.reports}; do
            if [ -s "$f" ]; then
                cat "$f"
            else
                sample=$(basename "$f" .mlst)
                echo -e "$sample\t\t\t\t\t\t\t\t\t"
            fi
        done > {output.mlst_report}
        """


### CheckM Report
rule checkm_summary:
    input:
        reports=expand(
            ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv", sample=Samples
        ),
    output:
        checkm_report=ISOLATE_FP / "reports" / "checkm.report",
    script:
        "scripts/summarize_checkm.py"


### AbritAMR Report:
rule abritamr_summary:
    input:
        reports=expand(
            ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out", sample=Samples
        ),
    output:
        amr_report=ISOLATE_FP / "reports" / "amr.report",
    script:
        "scripts/summarize_amr_output.py"


### Bakta Report
rule bakta_summary:
    input:
        reports=expand(
            ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt", sample=Samples
        ),
    output:
        bakta_report=ISOLATE_FP / "reports" / "bakta.report",
    script:
        "scripts/summarize_bakta.py"


## Shovill Report
rule shovill_summary:
    input:
        assemblies=expand(
            ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa", sample=Samples
        ),
    output:
        shovill_report=ISOLATE_FP / "reports" / "shovill.report",
    script:
        "scripts/summarize_shovill.py"


## Mash Report
rule mash_summary:
    input:
        reports=expand(
            ISOLATE_FP / "mash" / "{sample}_sorted_winning.tab", sample=Samples
        ),
    output:
        mash_report=ISOLATE_FP / "reports" / "mash.report",
    script:
        "scripts/summarize_mash.py"


## Summary Report
rule all_summary:
    input:
        reports=expand(ISOLATE_FP / "reports" / "{tool}.report", tool=TOOLS),
    output:
        final_report=ISOLATE_FP / "final_summary.tsv",
    script:
        "scripts/summarize_all.py"
