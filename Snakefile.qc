import gbshacks

config["all_samples"] = gbshacks.keyfile_to_samplenames(config['keyfile'])

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")


rule all:
    input:
        dynamic("data/reads/raw/{samp}_il.fastq.gz")

rule axe_demux:
    input:
        r1=config['qc']['inputs_r1'],
        r2=config['qc']['inputs_r2'],
        keyfile=config['keyfile'],
    output:
        files=expand("data/tmp/demuxed/{samp}_il.fastq", samp=config['all_samples']),
        statsfile="data/{}_demux-stats.tab".format(config['analysis_name']),
    threads:
        1
    log:
        "data/logs/axe/demux.log"
    params:
        mode=lambda wc: "-c" if config['lanes'][wc.lane].get("combo", True) else "",
    shell:
        "axe-demux "
        "   {params.mode}"
        "   -m 0"  # No mismatches allowable with borevitz lab barcodes
        "   -t {output.statsfile}"
        "   -b {input.keyfile}"
        "   -f <(cat {input.r1})"
        "   -r <(cat {input.r2})"
        "   -I data/tmp/demuxed/"
        " >{log} 2>&1"
