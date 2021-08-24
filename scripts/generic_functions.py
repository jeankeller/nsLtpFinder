#!/usr/bin/python3
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
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def clean_fasta(fasta, path_clean_fasta, species_code):
    with open(path_clean_fasta + path.sep + species_code + "_clean.fasta", "w") as fasta_clean:
        for rec in SeqIO.parse(fasta, "fasta"):
            clean_rec = SeqRecord(seq=rec.seq.upper(), id=rec.id, description="")
            SeqIO.write(clean_rec, fasta_clean, "fasta-2line")


def get_sequences_id_tp_extract(species_code, path_hmmer_res, path_ms_res):
    """
    Function to extract IDs of proteins from HMMSEARCH & motifSearch results for extracting their sequences later
    :param species_code: six-letter prefix species-specific
    :param path_hmmer_res: path to HMMSEARCH results
    :param path_ms_res: path to motifSearch results
    :return: list of protein IDs to keep and dataframe summarizing results (will be used downstream and written at
    the end)
    """
    # Create empty dataframe and read HMMSEARCH/motifSearch results
    df_res_summary = pds.DataFrame(columns=["seqName"], dtype=str)
    hmmer_res = pds.read_csv(glob(f"{path_hmmer_res}{path.sep}{species_code}*.domtblout")[0], delim_whitespace=True,
                             comment="#", names=["seqName", "accession", "seqLength", "queryName", "pfamCode",
                                                 "pfamLength", "fullSeqEvalue", "fullSeqScore", "fullSeqBias",
                                                 "domNum", "of", "domCEvalue", "domIEvalue", "domScore", "domBias",
                                                 "hmmStart", "hmmEnd", "aliStart", "aliEnd", "envStart", "envEnd",
                                                 "acc", "targetDesc"])
    ms_res = pds.read_csv(glob(f"{path_ms_res}{path.sep}{species_code}_motifs_count.txt")[0], sep="\t",
                          names=["seqIndex", "seqName", "motifRegex", "NumberOfMotifs"])

    # Get protein with motif identified
    non_null_motifsearch_res = ms_res[ms_res["NumberOfMotifs"] != 0].copy()
    # Get intersection of motifSearch & HMMSEARCH results
    seq_to_keep = list(set(hmmer_res["seqName"].tolist()) | set(non_null_motifsearch_res["seqName"].tolist()))
    df_res_summary["seqName"] = seq_to_keep
    # Add, for each protein if it has been detected by motifSearch (1 in motifsearch column), hmmsearch (1 in
    # hmmsearch column) or both (1 in both column). 0 indicates no result for the considered analysis
    df_res_summary["hmmsearch"] = df_res_summary["seqName"].\
        apply(lambda x: 1 if x in hmmer_res["seqName"].tolist() else 0)
    df_res_summary["motifsearch"] = df_res_summary["seqName"].\
        apply(lambda x: 1 if x in non_null_motifsearch_res["seqName"].tolist() else 0)
    return seq_to_keep, df_res_summary


def extract_seq(species_code, fasta_file, path_out, list_ids, prefix_file):
    """
    Function to extract sequences from list of IDs
    :param species_code: six-letter prefix species-specific
    :param fasta_file: path to protein file (FASTA format)
    :param path_out: path to output directory
    :param list_ids: list of IDs to extract
    :param prefix_file: prefix to determine if out file is for all candidates or only top candidates
    :return: write a FASTA file of sequences to extract
    """
    dbseq = SeqIO.index(fasta_file, "fasta")
    with open(path_out + path.sep + f"{species_code}_{prefix_file}_retained_seq.fa", "w") as fasta_out:
        for seq_id in list_ids:
            rec = dbseq.get(seq_id, default=None)
            if rec:
                SeqIO.write(rec, fasta_out, "fasta-2line")
            else:
                raise ValueError(f"{seq_id} not found in fasta file: {fasta_file}")
    fasta_out.close()


def get_prot_properties(s, retained_seq, cat):
    """
    Function to add column to the summarizing dataframe with protein properties (isoelectric point, gravy and
    molecular weight so far, other can be added)
    :param s: dataframe
    :param retained_seq: FASTA file of sequences to analyze
    :param cat: suffix to append at the end of column name if protein properties is made on different protein versions
    (e.g mature and non mature)
    :return: dataframe with new columns containing protein properties
    """
    seq = ProteinAnalysis(clean_seq(str(retained_seq[s["seqName"]].seq)))
    s["isoelectric_" + cat] = seq.isoelectric_point()
    s["gravy_" + cat] = seq.gravy()
    s["molecular_weight_" + cat] = seq.molecular_weight()
    return s


def clean_seq(seq):
    """
    Function to clean protein sequences prior propeties analysis. Trim invalid characters such as *, . and X
    :param seq: sequence to clean
    :return: cleaned sequence
    """
    return seq.replace("*", "").replace("X", "").replace(".", "")


def count_cysteines(seq):
    return seq.upper().count("C")


def merge_fasta(path_fasta_dir, species_code):
    """
    Function to merge fasta files
    :param path_fasta_dir: path to directory containing FASTA files to merge
    :param species_code: six-letter prefix species-specific
    :return: number of sequences after merging
    """
    list_fasta_to_merge = [SeqIO.index(x, "fasta") for x in
                           glob(f"{path_fasta_dir}{path.sep}{species_code}*_matureProt_retained_seq.fa")]
    fasta_merged = open(f"{path_fasta_dir}{path.sep}{species_code}_top_and_low_conf_candidates_matureProt.fa", "w")
    tmp_index = {k: v for d in list_fasta_to_merge for k, v in d.items()}
    for rec in tmp_index:
        SeqIO.write(tmp_index[rec], fasta_merged, "fasta-2line")
    return len(tmp_index)

