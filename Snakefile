import snkmk
import json
configfile: "config.yml"

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])
SAMPLESETS = snkmk.make_samplesets()


localrules: all, qc, map, varcall, qcreads

rule varcall:
    input:
        expand("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.{ext}",
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"],
               ext=["bcf", "bcf.csi", "vcf.gz", "vcf.gz.csi"]),

rule map:
    input:
        expand("data/alignments/{aligner}/{ref}/samples/{sample}.bam",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sample=SAMPLESETS['all-samples']),
        expand("data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sampleset=SAMPLESETS),

rule qc:
    input:
        expand("data/reads/qc/{sample}.fastq.gz", sample=SAMPLESETS['all-samples'])

rule all:
    input:
        rules.qc.input, rules.map.input, rules.varcall.input

#######################################################################
#                           Read processing                           #
#######################################################################


rule qcreads:
    input:
        reads="rawdata/samples/{sample}.fastq.gz"
    output:
        reads="data/reads/qc/{sample}.fastq.gz",
        settings="data/stats/adapterremoval/{sample}_settings.txt",
    log:
        "data/log/adapterremoval/{sample}.log",
    threads:
        2
    params:
        adp1=config["qc"]["adapter1"],
        adp2=config["qc"]["adapter2"],
        minqual=config["qc"]["minqual"],
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
        "   --minquality {params.minqual}"
        "   --threads {threads}"
        "   --settings {output.settings}"
        "   --output1 {output.reads}"
        " >{log} 2>&1"


#######################################################################
#                       Alignment to reference                        #
#######################################################################

rule ngmap:
    input:
        reads="data/reads/qc/{sample}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref]
    output:
        bam="data/alignments/ngm/{ref}/samples/{sample}.bam",
    log:
        "data/log/ngm/{ref}/{sample}.log"
    threads:
        4
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
        " ) >{log} 2>&1"

rule bwamem:
    input:
        reads="data/reads/qc/{sample}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref]
    output:
        bam="data/alignments/bwa/{ref}/samples/{sample}.bam",
    log:
        "data/log/bwa/{ref}/{sample}.log"
    threads:
        4
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
        " ) >{log} 2>&1"

rule mergebam:
    input:
        lambda wc: expand("data/alignments/{aligner}/{ref}/samples/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),
    output:
        bam="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
    log:
        "data/log/mergebam/{aligner}/{ref}/{sampleset}.log"
    threads: 8
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   {output.bam}"
        "   {input}"
        " ) >{log} 2>&1"

rule bamidx:
    input:
        "{path}.bam"
    output:
        "{path}.bam.csi"
    shell:
        "samtools index -c {input}"

#######################################################################
#                           Variant Calling                           #
#######################################################################

rule freebayes:
    input:
        bam="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
        bamidx="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam.csi",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/freebayes~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/freebayes/{aligner}~{ref}~{sampleset}/{region}.log"
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
    shell:
        "( samtools view"  # some versions of freebayes don't seek to region
        "   -u"            # also region format is zero-based in freebayes
        "    {input.bam}"  # so we extract the correct region from the BAM
        "   '{wildcards.region}'"
        " | freebayes"
        "   --theta {params.theta}"
        "   --use-best-n-alleles 4"
        "   --min-alternate-fraction 0"
        "   --min-alternate-count 1" # per sample
        "   --min-alternate-total 3" # across all samples
        "   --min-coverage 5" # across all samples
        "   --strict-vcf"
        "   --stdin"
        "   --fasta-reference {input.ref}"
        " | bcftools view"
        "   -O u  -o {output.bcf}"
        " ) >{log} 2>&1"

rule mpileup:
    input:
        bam="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
        bamidx="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam.csi",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/mpileup~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/mpileup/{aligner}~{ref}~{sampleset}/{region}.log"
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
    shell:
        "( samtools mpileup"
        "   --output-tags DP,AD,SP,INFO/AD" #output everything
        "   --region '{wildcards.region}'"
        "   --fasta-ref {input.ref}"
        "   --redo-BAQ"
        "   --BCF --uncompressed"
        "   {input.bam}"
        " | bcftools call"
        "   --targets '{wildcards.region}'" # might not be needed
        "   --multiallelic-caller"
        "   --prior {params.theta}"
        "   -O u"
        "   -o {output.bcf}"
        " ) >{log} 2>&1"


rule bcfnorm:
    input:
        bcf="data/variants/raw_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        # Not a pipe! can't run multiple filters if a pipe
        bcf=temp("data/variants/norm_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf"),
    log:
        "data/log/bcfnormalise/{caller}~{aligner}~{ref}~{sampleset}/{region}.log"
    shell:
        "( bcftools norm"
        "   --fasta-ref {input.ref}"
        "   -O u"
        "   {input.bcf}"
        " | vt decompose_blocksub + -o -" # decompose MNP to multipe SNPs
        " | bcftools norm" # Split multi-alleics
        "   --do-not-normalize"
        "   --multiallelics -snps"
        "   -O u  -o {output.bcf}"
        " ) >{log} 2>&1"

rule bcffilter:
    input:
        bcf="data/variants/norm_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        # Not a pipe! can't run all regions separately if this is a pipe into merge
        bcf=temp("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf"),
    log:
        "data/log/bcffilter/{caller}~{aligner}~{ref}~{sampleset}/{filter}/{region}.log"
    params:
        filtarg=lambda wc: config["varcall"]["filters"][wc.filter].replace('\n', ' ')
    shell:
        "( bcftools view"
        "   {params.filtarg}"
        "   -O u"
        "   {input.bcf}"
        " | bcftools norm" # We normalise here to re-join multi-allelic sites, after filtering with multi-allelics split
        "   --fasta-ref {input.ref}"
        "   --do-not-normalize"
        "   --multiallelics +snps" # Split multi-alleic sites
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"

localrules: bcfmerge_fofn
rule bcfmerge_fofn:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.ref])),
    output:
        fofn=temp("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN"),
    run:
        with open(output[0], "w") as fh:
            for s in input:
                print(s, file=fh)

rule bcfmerge:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.ref])),
        fofn="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN",
    output:
        bcf="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf",
    log:
        "data/log/mergebcf/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}.log"
    threads: 8
    shell:
        "( bcftools concat"
        "   --threads {threads}"
        "   -O b"
        "   -o {output.bcf}"
        "   --file-list {input.fofn}"
        " ) >{log} 2>&1"


rule bcf2vcf:
    input:
        bcf="{path}.bcf",
    output:
        vcf="{path}.vcf.gz",
    log:
        "data/log/bcf2vcf/{path}.log"
    threads: 8
    shell:
        "( bcftools view"
        "   {input.bcf}"
        "   -O z"
        "   --threads {threads}"
        "   -o {output.vcf}"
        " ) >{log} 2>&1"

rule variantidx:
    input:
        "{path}.{ext}"
    output:
        "{path}.{ext}.csi"
    wildcard_constraints:
        ext="bcf|vcf.gz"
    priority: 2
    shell:
        "bcftools index -f {input}"


