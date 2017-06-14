#!/usr/bin/env python3
import yaml
import subprocess as sp
from multiprocessing.pool import ThreadPool
from sys import stdin, stdout, stderr, argv
import sys
import os


def parse_cfg():
    with open("./config.yml") as fh:
        d = yaml.load(fh)
    return d["lanes"]


def runaxe(lane):
    lanedata = parse_cfg()[lane]
    keyfile = "metadata/keyfiles/" + lane + ".axe"
    statsfile = "data/stats/demux/" + lane + ".tsv"
    r1fq = "rawdata/gbs/" + lane + "/" + lane + "_R1.fastq.gz"
    r2fq = "rawdata/gbs/" + lane + "/" + lane + "_R2.fastq.gz"
    outprefix = "data/reads/raw/" + lane + "/"
    logfile = "data/log/demux/" + lane + ".log"

    os.makedirs("data/stats/demux", exist_ok=True)
    os.makedirs("data/log/demux", exist_ok=True)
    os.makedirs(outprefix.rstrip("/"), exist_ok=True)

    cmd = ["axe-demux",
           "-m", "0",
           "-t", statsfile,
           "-b", keyfile,
           "-f", r1fq,
           "-r", r2fq,
           "-I", outprefix]

    if lanedata["combo"]:
        cmd.append("-c")

    print("Running lane", lane, file=stderr)
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)

    with open(logfile, "wb", buffering=1000) as fh:
        fh.write((" ".join(cmd) + "\n").encode("utf8"))
        while True:
            fh.write(proc.stdout.read(100))
            if proc.poll() is not None:
                break
    print("Finished lane", lane, "(returned", proc.returncode, ")", file=stderr)
    return proc.returncode


def main():
    lanes = parse_cfg()
    if len(argv) > 1:
        # Filter on whitelist of plates
        lanes = {k: v for k, v in lanes.items() if k in argv[1:]}

    pool = ThreadPool(16)
    rets = pool.map(runaxe, lanes)
    sys.exit(max(rets))

if __name__ == "__main__":
    main()
