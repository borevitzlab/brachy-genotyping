from __future__ import print_function
from collections import defaultdict
import csv
from glob import glob
from os.path import (
    basename,
    dirname,
    splitext,
)


def parsefai(fai):
    with open(fai) as fh:
        for l in fh:
            cname, clen, _, _, _ = l.split()
            clen = int(clen)
            yield cname, clen


def make_regions(rdict, window=1e6, base=1):
    window = int(window)
    ret = {}
    for refname, refpath in rdict.items():
        fai = refpath+".fai"
        windows = []
        curwin = []
        curwinlen = 0
        for cname, clen in parsefai(fai):
            for start in range(0, clen, window):
                wlen = min(clen - start, window)
                windows.append("{}:{:09d}-{:09d}".format(cname, start + base, start+wlen))
        ret[refname] = windows
    return ret


def make_chroms(rdict):
    ret = {}
    for refname, refpath in rdict.items():
        fai = refpath+".fai"
        ref = dict()
        scafs = []
        for cname, clen in parsefai(fai):
            if cname.lower().startswith("chr"):
                ref[cname] = [cname]
            else:
                scafs.append(cname)
        ref["scaffolds"] = scafs
        ret[refname] = ref
        print(refname, "has", len(ref), "chromosome sets")
    return ret


def make_readcountdict(lanes):
    readcounts = dict()
    for lane in lanes:
        try:
            with open("data/stats/demux/{}.tsv".format(lane)) as fh:
                next(fh)  # skip header
                for line in fh:
                    line = line.rstrip().split("\t")
                    samp = line[-2]
                    if samp == "No Barcode":
                        continue
                    rcount = int(line[-1])
                    readcounts[samp] = rcount
        except FileNotFoundError as exc:
            pass
    return readcounts


def s2l2s(metadatafile):
    samples = csv.DictReader(open(metadatafile))
    s2l = {}
    l2s = defaultdict(list)
    for sample in samples:
        s2l[sample["anon.name"]] = sample["lane"]
        l2s[sample["lane"]].append(sample["anon.name"])
    return s2l, l2s

def readlines(file):
    with open(file) as fh:
        return [l.rstrip("\n\r") for l in fh]

def make_samplesets():
    sets = dict()
    all_samples = list()
    for sset in glob("metadata/samplesets/*.txt"):
        setname = splitext(basename(sset))[0]
        samples =  readlines(sset)
        sets[setname] = samples
        all_samples.extend(samples)
    sets["all-samples"] = list(set(all_samples))
    return sets
