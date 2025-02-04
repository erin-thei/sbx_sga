def get_sga_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_sga":
            return Path(fp)
    raise Error(
        "Filepath for sbx_sga not found, are you sure it's installed under extensions/sbx_sga?"
    )


ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
SBX_SBA_VERSION = open(get_sga_path() / "VERSION").read().strip()

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
        ISOLATE_FP / "checkm" / "quality_report.tsv",
        expand(ISOLATE_FP / "quast" / "{sample}" / "report.tsv", sample=Samples),
        # Typing
        ISOLATE_FP / "mlst" / "mlst_report.tsv",
        # Annotation
        expand(ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt", sample=Samples),
        # AMR Profiling
        expand(ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out", sample=Samples),


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
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "contigs.fa",
    log:
        LOG_FP / "sga_shovill_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_shovill_{sample}.tsv"
    conda:
        "envs/shovill.yml"
    shell:
        """
        shovill --assembler skesa --outdir $(dirname {output.contigs}) --R1 {input.rp1}  --R2 {input.rp2} &> {log}
        """


### Assembly QC
rule sga_checkm:
    input:
        contigs=expand(
            ISOLATE_FP / "shovill" / "{sample}" / "contigs.fa", sample=Samples
        ),
    output:
        quality_report=ISOLATE_FP / "checkm" / "quality_report.tsv",
    params:
        ref=Cfg["sbx_sga"]["checkm_ref"],
    log:
        LOG_FP / "sga_checkm.log",
    benchmark:
        BENCHMARK_FP / "sga_checkm.tsv"
    conda:
        "envs/checkm2.yml"
    shell:
        """
        checkm2 predict \\
        -x fa \\
        -i $(dirname {input.contigs[0]}) \\
        -o $(dirname {output.quality_report}) \\
        --database_path {params.ref} \\
        &> {log}
        """


rule sga_quast:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "contigs.fa",
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
        quast.py \\
        -o $(dirname {output.quast_dir}) \\
        {input.contigs} \\
        &> {log}
        """


### Typing
rule sga_mlst:
    input:
        contigs=expand(
            ISOLATE_FP / "shovill" / "{sample}" / "contigs.fa", sample=Samples
        ),
    output:
        mlst=ISOLATE_FP / "mlst" / "mlst_report.tsv",
    log:
        LOG_FP / "sga_mlst.log",
    benchmark:
        BENCHMARK_FP / "sga_mlst.tsv"
    conda:
        "envs/mlst.yml"
    shell:
        """
        mlst {input.contigs} > {output.mlst} 2> {log}
        """


### Annotation
rule sga_bakta:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "contigs.fa",
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
        bakta --db {params.ref} \\
        --output $(dirname {output.bakta}) \\
        --prefix {wildcards.sample} \\
        --skip-plot {input.contigs} \\
        &> {log}
        """


### AMR Profiling
rule sga_abritamr:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "contigs.fa",
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
        abritamr run \\
        --contigs {input.contigs} \\
        --prefix {output.abritamr} \\
        &> {log}
        """
