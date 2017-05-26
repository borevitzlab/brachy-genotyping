import json
configfile: "config.yml"

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

SAMP2LANE = json.load(open("metadata/samp2lane.json"))


rule all:
    input:
        expand("data/reads/qc/{sample}.fastq.gz", sample=SAMP2LANE),


rule axe_demux:
    input:
        r1="rawdata/gbs/{lane}/{lane}_R1.fastq.gz",
        r2="rawdata/gbs/{lane}/{lane}_R2.fastq.gz",
        keyfile="metadata/keyfiles/{lane}.axe"
    output:
        files=dynamic("data/reads/raw/{lane}/{sample}_il.fastq.gz"),
    threads:
        1
    log:
        "data/logs/axe/{lane}.log"
    params:
        statsfile="data/stats/axe/{lane}.tsv",
        mode=lambda wc: "-c" if config['lanes'][wc.lane].get("combo", True) else "",
    shell:
        "mkdir -p $(dirname {params.statsfile}); "
        "axe-demux "
        "   {params.mode}"
        "   -m 0"  # No mismatches allowable with borevitz lab barcodes
        "   -t {params.statsfile}"
        "   -b {input.keyfile}"
        "   -f <(cat {input.r1})"
        "   -r <(cat {input.r2})"
        "   -I data/reads/raw/{wildcards.lane}/"
        "   -z 1" # zip them
        " >{log} 2>&1"

rule qcreads:
    input:
        reads=lambda wc: "data/reads/raw/{lane}/{sample}_il.fastq.gz".format(lane=SAMP2LANE[wc.sample], sample=wc.sample),
    output:
        reads="data/reads/qc/{sample}.fastq.gz",
        settings="data/log/adapterremoval/{sample}_settings.txt",
    log:
        "data/log/adapterremoval/{sample}.log",
    threads:
        2
    params:
        adp1=config["qc"]["adapter1"],
        adp2=config["qc"]["adapter2"],
    shell:
        "AdapterRemoval"
        "   --file1 <(cat {input.r1})"
        "   --file2 <(cat {input.r2})"
        "   --adapter1 {params.adp1}"
        "   --adapter2 {params.adp2}"
        "   --combined-output"
        "   --interleaved"
        "   --gzip"
        "   --collapse"
        "   --trimns"
        "   --trimqualities"
        "   --threads {threads}"
        "   --settings {output.settings}"
        "   --output1 {output.reads}"
        " >{log} 2>&1"
