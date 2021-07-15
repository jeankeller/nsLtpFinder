# -*- coding: utf-8 -*-
# @Time    : 13/07/2021 13:50
# @Author  : Jean Keller
# @Email   : kellerjeanphd@gmail.com
# @File    : generic_functions.py
# @Software: PyCharm

import pandas as pds
from glob import glob
from os import path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def get_sequences_id_tp_extract(species_code, path_hmmer_res, path_ms_res):
    df_res_summary = pds.DataFrame(columns=["seqName"], dtype=str)
    hmmer_res = pds.read_csv(glob(f"{path_hmmer_res}{path.sep}{species_code}*.domtblout")[0], delim_whitespace=True,
                             comment="#", names=["seqName", "accession", "seqLength", "queryName", "pfamCode",
                                                 "pfamLength", "fullSeqEvalue", "fullSeqScore", "fullSeqBias",
                                                 "domNum", "of", "domCEvalue", "domIEvalue", "domScore", "domBias",
                                                 "hmmStart", "hmmEnd", "aliStart", "aliEnd", "envStart", "envEnd",
                                                 "acc", "targetDesc"])
    ms_res = pds.read_csv(glob(f"{path_ms_res}{path.sep}{species_code}_motifs_count.txt")[0], sep="\t",
                          names=["seqIndex", "seqName", "motifRegex", "NumberOfMotifs"])
    non_null_motifsearch_res = ms_res[ms_res["NumberOfMotifs"] != 0].copy()
    seq_to_keep = list(set(hmmer_res["seqName"].tolist()) | set(non_null_motifsearch_res["seqName"].tolist()))
    df_res_summary["seqName"] = seq_to_keep
    df_res_summary["hmmsearch"] = df_res_summary["seqName"].\
        apply(lambda x: 1 if x in hmmer_res["seqName"].tolist() else 0)
    df_res_summary["motifsearch"] = df_res_summary["seqName"].\
        apply(lambda x: 1 if x in non_null_motifsearch_res["seqName"].tolist() else 0)
    return seq_to_keep, df_res_summary


def extract_seq(species_code, fasta_file, path_out, list_ids, prefix_file):
    dbseq = SeqIO.index(fasta_file, "fasta")
    with open(path_out + path.sep + f"{species_code}_{prefix_file}_retained_seq.fa", "w") as fasta_out:
        for seq_id in list_ids:
            rec = dbseq.get(seq_id, default=None)
            if rec:
                SeqIO.write(rec, fasta_out, "fasta-2line")
            else:
                raise ValueError(f"{seq_id} not found in fasta file: {fasta_file}")
    fasta_out.close()


def get_prot_properties(s, retained_seq):
    seq = ProteinAnalysis(clean_seq(str(retained_seq[s["seqName"]].seq)))
    s["isoelectric"] = seq.isoelectric_point()
    s["gravy"] = seq.gravy()
    s["molecular_weight"] = seq.molecular_weight()
    return s


def clean_seq(seq):
    return seq.replace("*", "").replace("X", "").replace(".", "")