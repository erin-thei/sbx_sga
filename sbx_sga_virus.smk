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
            VIRUS_FP / "virsorter" / "{sample}" / "final-viral-combined.fa",
            sample=Samples,
        ),


rule sga_install_virsorter:
    output:
        VIRUS_FP / "virsorter" / ".installed",
    benchmark:
        BENCHMARK_FP / "install_virsorter.tsv"
    log:
        LOG_FP / "install_virsorter.log",
    params:
        db_fp=Cfg["sbx_sga"]["virsorter_ref"],
    resources:
        runtime=2400,
    threads: 4
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        # First check if directory exists and has files
        if [ -d {params.db_fp} ] && [ "$(ls -A {params.db_fp})" ]; then
            echo "VirSorter database already installed" > {log}
            touch {output}
            exit 0
        fi

        echo "Installing VirSorter database" > {log}
        virsorter setup -d {params.db_fp} -j {threads} > {log} 2>&1
        touch {output}
        """


rule sga_virsorter:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
        install=VIRUS_FP / "virsorter" / ".installed",
    output:
        combined_viral=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-combined.fa",
        scores=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-score.tsv",
        boundaries=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-boundary.tsv",
    benchmark:
        BENCHMARK_FP / "virsorter_{sample}.tsv"
    log:
        LOG_FP / "virsorter_{sample}.log",
    params:
        out_dir=str(VIRUS_FP / "virsorter" / "{sample}"),
        db_fp=Cfg["sbx_sga"]["virsorter_ref"],
    resources:
        mem_mb=24000,
        runtime=720,
    threads: 4
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        virsorter \
        run \
        -w {params.out_dir} \
        -i {input.contigs} \
        --min-length 1000 \
        -j {threads} \
        --db-dir {params.db_fp} \
        all \
        > {log} 2>&1
        """


rule sga_install_checkv:
    output:
        VIRUS_FP / "checkv" / ".installed",
    benchmark:
        BENCHMARK_FP / "install_checkv.tsv"
    log:
        LOG_FP / "install_checkv.log",
    params:
        db_fp=Cfg["sbx_sga"]["checkv_ref"],
    resources:
        runtime=2400,
    threads: 4
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        # First check if directory exists and has files
        if [ -d {params.db_fp} ] && [ "$(ls -A {params.db_fp})" ]; then
            echo "CheckV database already installed" > {log}
            touch {output}
            exit 0
        fi

        echo "Installing CheckV database" > {log}
        checkv download_database {params.db_fp} >> {log} 2>&1
        touch {output}
        """


rule sga_checkv:
    input:
        contigs=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-combined.fa",
        install=VIRUS_FP / "checkv" / ".installed",
    output:
        proviruses=VIRUS_FP / "checkv" / "{sample}" / "proviruses.fna",
        viruses=VIRUS_FP / "checkv" / "{sample}" / "viruses.fna",
    benchmark:
        BENCHMARK_FP / "checkv_{sample}.tsv"
    log:
        LOG_FP / "checkv_{sample}.log",
    params:
        db_fp=Cfg["sbx_sga"]["checkv_ref"],
    resources:
        mem_mb=24000,
        runtime=720,
    threads: 4
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        checkv \
        end_to_end \
        {input.contigs} \
        $(dirname (output.proviruses)) \
        -t {threads} \
        -d {params.db_fp} \
        > {log} 2>&1
        """


rule sga_combine_checkv:
    input:
        proviruses=VIRUS_FP / "checkv" / "{sample}" / "proviruses.fna",
        viruses=VIRUS_FP / "checkv" / "{sample}" / "viruses.fna",
    output:
        virus=VIRUS_FP / "checkv" / "{sample}" / "combined.fna",
    shell:
        """
        cat {input.proviruses} {input.viruses} > {output.virus}
        """


rule sga_pharokka_install:
    output:
        VIRUS_FP / "pharokka" / ".installed",
    benchmark:
        BENCHMARK_FP / "install_pharokka.tsv"
    log:
        LOG_FP / "install_pharokka.log",
    params:
        db_fp=Cfg["sbx_sga"]["pharokka_ref"],
    resources:
        runtime=2400,
    threads: 4
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        # First check if directory exists and has files
        if [ -d {params.db_fp} ] && [ "$(ls -A {params.db_fp})" ]; then
            echo "Pharokka database already installed" > {log}
            touch {output}
            exit 0
        fi

        echo "Installing Pharokka database" > {log}
        install_databases.py -o {params.df_fp} >> {log} 2>&1
        touch {output}
        """


rule sga_pharokka:
    input:
        contigs=VIRUS_FP / "checkv" / "{sample}" / "combined.fna",
        install=VIRUS_FP / "pharokka" / ".installed",
    output:
        pharokka=VIRUS_FP / "pharokka" / "{sample}" / "pharokka.tsv",
    benchmark:
        BENCHMARK_FP / "pharokka_{sample}.tsv"
    log:
        LOG_FP / "pharokka_{sample}.log",
    params:
        db_fp=Cfg["sbx_sga"]["pharokka_ref"],
    resources:
        mem_mb=24000,
        runtime=720,
    threads: 4
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        pharokka.py -i <phage fasta file> -o <output directory> -d <path/to/database_dir> -t <threads>
        """