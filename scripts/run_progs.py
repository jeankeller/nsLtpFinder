#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 13/07/2021 10:29
# @Author  : Jean Keller
# @Email   : kellerjeanphd@gmail.com
# @File    : run_progs.py
# @Software: PyCharm

import subprocess
from os import path
from scripts import motif_finder


def run_hmmsearch(fasta_query, path_rep_out, global_evalue, domain_evalue, hmm_models, hmm_cpus):
    """
    Function to perform HMMSEARCH on a protein FASTA file
    Search is performed with a default e-value of 1e-04 for BOTH hit and domain, can be changed
    :param hmm_cpus: number of cpus to use for the HMMSEARCH algorithm
    :param hmm_models: path to a text file containing all HMMs model to find
    :param domain_evalue: e-value threshold for domain alignment
    :param global_evalue: e-value threshold for the whole alignment
    :param path_rep_out: path to output directory where result files will be written
    :param fasta_query: path to query proteins in FASTA format
    :return: write three results files : (1) .res: HMMSEARCH results in human-readable format, (2) .domtblout: HMMSEARCH
    results in tabular format (use for computational parsing), (3) .ali: alignments between proteins and models
    """
    filename = path.splitext(path.basename(fasta_query))[0].split("_")[0]
    domtblout_out = path_rep_out + path.sep + "{0}_hmmsearch.domtblout".format(filename)
    hmmsearch_res = path_rep_out + path.sep + "{0}_hmmsearch.res".format(filename)
    alignment_out = path_rep_out + path.sep + "{0}_hmmsearch.ali".format(filename)
    hmmsearch_args = ["hmmsearch", "--domtblout", domtblout_out, "-A", alignment_out, "-o", hmmsearch_res, "-E",
                      str(global_evalue), "--domE", str(domain_evalue), "--cpu", str(hmm_cpus), hmm_models, fasta_query]
    subprocess.run(hmmsearch_args)


def run_motifsearch(fasta_query, output_dir, species_code):
    """
    Function to run motifSearch package without installing it
    :param fasta_query: path to query proteins in FASTA format
    :param output_dir: path to output directory where result files will be written
    :return: write all results in output directory
    """
    fcount = open(output_dir + path.sep + f"{species_code}_motifs_count.txt", "w")
    fpos = open(output_dir + path.sep + f"{species_code}_motifs_positions.txt", "w")
    fseq = open(output_dir + path.sep + f"{species_code}_motifs_sequences.fa", "w")
    search_motif = motif_finder.GetMotifs(motif="C.{6,14}C.{8,19}CC.{9,20}C.C.{12,34}C.{5,14}C", fasta_file=fasta_query,
                                          path_count=fcount, fseq=fseq, fpos=fpos)
    search_motif.checkinput()
    search_motif.writemotif()


def run_signalp(fasta_query, path_out, species_code):
    """
    Function to run signalp program
    :param fasta_query: path to query proteins in FASTA format
    :param path_out: path to output directory
    :param species_code: six-letter prefix species-specific
    :return: write results (mature proteins, signalp log and signalp table results) in out directory
    """
    signalp_out = open(path_out + path.sep + species_code + ".signalp.log", "w")
    sp_args = ["signalp5", "-fasta", fasta_query, "-format", "short", "-mature", "-prefix", species_code, "-org", "euk"]
    subprocess.run(sp_args, stdout=signalp_out, cwd=path_out)


def run_meme(fasta_query, meme_res_out, species_code, analysis_type, cpus):
    """
    Function to run MEME program
    :param fasta_query: path to query proteins in FASTA format
    :param meme_res_out: path to output directoru
    :param species_code: six-letter prefix species-specific
    :param analysis_type: string to precise if the analysis concerns all candidate or top candidate
    :param cpus: number of threads to use. ONLY SET UP IF PARALLEL VERSION OF MEME AVAILABLE!!! (openMPI)
    :return: write all results (logo pictures in eps/png format, html and txt formatted results) in out directory
    """
    meme_out = open(meme_res_out + path.sep + species_code + "_" + analysis_type + ".log", "w")
    meme_args = ["meme", fasta_query, "-oc", meme_res_out + path.sep + species_code + "_" + analysis_type, "-protein",
                 "-mod", "anr", "-nmotifs", "10", "-minw", "10", "-maxw", "30", "-evt", "0.001"]
    subprocess.run(meme_args, stdout=meme_out, stderr=meme_out)
