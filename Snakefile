import snkmk
import json
configfile: "config.yml"

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

SAMP2LANE = json.load(open("metadata/samp2lane.json"))
LANE2SAMP = json.load(open("metadata/lane2samp.json"))
REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])

rule all:
    input:
        #expand("data/reads/qc/{lane}/{sample}.fastq.gz", sample=SAMP2LANE),
        expand("data/alignments/{aligner}/{ref}/{sample}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sample=SAMP2LANE),
        expand("data/alignments/{aligner}/{ref}.bam",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"]),
#        expand("data/variants/{caller}/{aligner}/{ref}.bcf",
#               caller=config["varcall"]["callers"],
#               aligner=config["varcall"]["aligners"],
#               ref=config["varcall"]["refs"]),

rule qcreads:
    input:
        reads="data/reads/raw/{lane}/{sample}_il.fastq"
    output:
        reads="data/reads/qc/{lane}/{sample}.fastq.gz",
        settings="data/stats/adapterremoval/{lane}/{sample}_settings.txt",
    log:
        "data/log/adapterremoval/{lane}/{sample}.log",
    threads:
        4
    params:
        adp1=config["qc"]["adapter1"],
        adp2=config["qc"]["adapter2"],
    shell:
        "AdapterRemoval"
        "   --file1 {input.reads}"
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

rule identify_adaptors:
    input:
        idlog=expand("data/qcstats/adatpter-id/{lane}.txt", lane=LANE2SAMP)

rule ar_id_adaptors:
    input:
        r1="rawdata/gbs/{lane}/{lane}_R1.fastq.gz",
        r2="rawdata/gbs/{lane}/{lane}_R2.fastq.gz",
    output:
        idlog="data/qcstats/adatpter-id/{lane}.txt",
    log:
        "data/log/adapterremoval-idadapt/{lane}.log",
    threads:
        4
    params:
        adp1=config["qc"]["adapter1"],
        adp2=config["qc"]["adapter2"],
    shell:
        "AdapterRemoval"
        "   --file1 {input.r1}"
        "   --file2 {input.r2}"
        "   --adapter1 {params.adp1}"
        "   --adapter2 {params.adp2}"
        "   --identify-adapters"
        "   --threads {threads}"
        " >{output.idlog} 2>{log}"

rule map:
    input:
        expand("data/alignments/{aligner}/{ref}/{sample}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sample=SAMP2LANE),
        expand("data/alignments/{aligner}/{ref}.bam",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"]),

rule ngmap:
    input:
        reads=lambda wc: "data/reads/qc/{lane}/{sample}.fastq.gz".format(
                                lane=SAMP2LANE[wc.sample], sample=wc.sample),
        ref=lambda wc: config['refs'][wc.ref]
    output:
        bam="data/alignments/ngm/{ref}/{sample}.bam",
        bai="data/alignments/ngm/{ref}/{sample}.bam.bai",
    log:
        "data/log/ngm/{ref}_{sample}.log"
    threads:
        8
    shell:
        "( ngm"
        "   -q {input.reads}"
        "   -p" # paired input
        "   -r {input.ref}"
        "   -t {threads}"
        "   --rg-id {wildcards.sample}"
        "   --rg-sm {wildcards.sample}"
        "   --very-sensitive"
        " | samtools view -Suh -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.sample}"
        "   -@ {threads}"
        "   -m 1G"
        "   -o {output.bam}"
        "   -" # stdin
        " && samtools index {output.bam}"
        " ) >{log} 2>&1"

rule bwamem:
    input:
        reads=lambda wc: "data/reads/qc/{lane}/{sample}.fastq.gz".format(
                                lane=SAMP2LANE[wc.sample], sample=wc.sample),
        ref=lambda wc: config['refs'][wc.ref]
    output:
        bam="data/alignments/bwa/{ref}/{sample}.bam",
        bai="data/alignments/bwa/{ref}/{sample}.bam.bai",
    log:
        "data/log/bwa/{ref}/{sample}.log"
    threads:
        8
    shell:
        "( bwa mem"
        "   -p" # paired input
        "   -t {threads}"
        "   -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}'"
        "   {input.ref}"
        "   {input.reads}"
        " | samtools view -Suh -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.sample}"
        "   -@ {threads}"
        "   -m 1G"
        "   -o {output.bam}"
        "   -" # stdin
        " && samtools index {output.bam}"
        " ) >{log} 2>&1"

rule mergebam:
    input:
        expand("data/alignments/{{aligner}}/{{ref}}/{sample}.bam", sample=SAMP2LANE),
    output:
        bam="data/alignments/{aligner}/{ref}.bam",
        bai="data/alignments/{aligner}/{ref}.bam.bai",
    log:
        "data/log/mergebam/{ref}.log"
    threads: 3
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   {output.bam}"
        "   {input}"
        " && samtools index {output.bam}"
        " ) >{log} 2>&1"



rule varcall:
    input:
        expand("data/variants/{caller}/{aligner}/{ref}.bcf",
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"]),


rule freebayes:
    input:
        bam="data/alignments/{aligner}/{ref}.bam",
        bai="data/alignments/{aligner}/{ref}.bam.bai",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/freebayes/{aligner}/{ref}/split/{region}.bcf",
        idx="data/variants/freebayes/{aligner}/{ref}/split/{region}.bcf.csi",
    log:
        "data/log/freebayes/{aligner}/{ref}/{region}.log"
    threads: 1
    params:
        region=lambda wc: "' --region '".join(REGIONS[wc.ref][wc.region])
    shell:
        "( freebayes"
        "   --theta 0.02" # higher prior on mutation rate
        "   --use-reference-allele"
        "   --min-mapping-quality 10"
        "   --min-base-quality 10"
        "   --min-alternate-fraction 0.1"
        "   --min-alternate-count 1"
        "   --min-alternate-total 4"
        "   --use-mapping-quality"
        "   --genotype-qualities"
        "   --region '{params.region}'"
        "   -f {input.ref}"
        "   {input.bam}"
        " | bcftools view"
        "   -O b"
        "   -o {output.bcf}"
        " && bcftools index -f {output.bcf}"
        " ) >{log} 2>&1"

rule mpileup:
    input:
        bam="data/alignments/{aligner}/{ref}.bam",
        bai="data/alignments/{aligner}/{ref}.bam.bai",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/mpileup/{aligner}/{ref}/split/{region}.bcf",
        idx="data/variants/mpileup/{aligner}/{ref}/split/{region}.bcf.csi",
    log:
        "data/log/mpileup/{aligner}/{ref}/{region}.log"
    threads: 1
    params:
        region=lambda wc: "' --region '".join(REGIONS[wc.ref][wc.region]),
        targets=lambda wc: "' --targets '".join(REGIONS[wc.ref][wc.region]) # for bcftools
    shell:
        "( samtools mpileup"
        "   --output-tags DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR" #output everything
        "   --region '{params.region}'"
        "   --fasta-ref {input.ref}"
        "   --redo-BAQ"
        "   --BCF --uncompressed"
        "   {input.bam}"
        " | bcftools call"
        "   --targets '{params.targets}'" # might not be needed
        "   --multiallelic-caller"
        "   --prior 0.01" # increase mutation rate prior
        "   -O b"
        "   -o {output.bcf}"
        " && bcftools index -f {output.bcf}"
        " ) >{log} 2>&1"

rule bcfmerge:
    input:
        bcf=lambda wc: expand("data/variants/{caller}/{aligner}/{ref}/split/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref,
                              region=sorted(REGIONS[wc.ref])),
        idx=lambda wc: expand("data/variants/{caller}/{aligner}/{ref}/split/{region}.bcf.csi",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref,
                              region=sorted(REGIONS[wc.ref])),
    output:
        bcf="data/variants/{caller}/{aligner}/{ref}.bcf",
    log:
        "data/log/merge/{caller}/{aligner}/{ref}.log"
    threads: 2
    shell:
        "( bcftools concat"
        "   --allow-overlaps"
        "   --remove-duplicates"
        "   --threads {threads}"
        "   -O b"
        "   -o {output.bcf}"
        "   {input.bcf}"
        " ) >{log} 2>&1"
