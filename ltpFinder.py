#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 12/07/2021 17:17
# @Author  : Jean Keller; University of Toulouse III
# @Email   : kellerjeanphd@gmail.com
# @File    : ltpFinder.py
# @Software: PyCharm

import sys
from os import path, mkdir, getcwd
from glob import glob
from scripts.check_dependencies import check_installed_programs, check_installed_modules
from scripts.run_progs import run_hmmsearch, run_motifsearch, run_signalp, run_meme
from scripts.generic_functions import *


def main():
    # Check if all modules are installed
    check_installed_modules()

    # import modules
    import argparse
    import shutil
    import pandas as pds
    from Bio import SeqIO

    # check if all programs are installed
    check_installed_programs()

    args_parser = argparse.ArgumentParser(usage="ltpFinder.py [options...]", description="""This script allows 
    identifying LTP in a given proteome""")
    args_parser.add_argument("-f", "--files", required=True, type=str, help="Path to directory containing fasta "
                                                                            "queries")
    args_parser.add_argument("-t", "--threads", required=False, type=int, default=1, help="Number of CPUs to use")
    args_parser.add_argument("--hmm", required=True, type=str, help="Path to file containing HMM models")
    args_parser.add_argument("-o", "--output", required=True, help="Path to output directory")
    args_parser.add_argument("--hmm_eval_glob", required=False, default=1e-04, help="E-value threshold for HMM whole "
                                                                                    "alignment")
    args_parser.add_argument("--hmm_eval_dom", required=False, default=1e-04, help="E-value threshold for HMM domain "
                                                                                   "alignment")
    args = args_parser.parse_args()

    out_dir = args.output
    if not path.isdir(out_dir):
        mkdir(out_dir)

    # Create working directory
    workdir = out_dir + path.sep + "workdir"
    if not path.isdir(workdir):
        mkdir(workdir)

    # Create HMM output directory
    hmmsearch_outdir = out_dir + path.sep + "hmm_res"
    if not path.isdir(hmmsearch_outdir):
        mkdir(hmmsearch_outdir)

    # Create motifSearch output directory
    ms_outdir = out_dir + path.sep + "motifSearch_res"
    if not path.isdir(ms_outdir):
        mkdir(ms_outdir)

    # Create signalP output directory
    path_signalp_out = out_dir + path.sep + "signalp_res"
    if not path.isdir(path_signalp_out):
        mkdir(path_signalp_out)

    # Create retained sequences output directory
    path_fasta_seq = out_dir + path.sep + "retained_sequences"
    if not path.isdir(path_fasta_seq):
        mkdir(path_fasta_seq)

    # Create MEME output directory
    path_meme_out = out_dir + path.sep + "meme_res"
    if not path.isdir(path_meme_out):
        mkdir(path_meme_out)

    # Create results output directory
    path_res_out = out_dir + path.sep + "final_results"
    if not path.isdir(path_res_out):
        mkdir(path_res_out)

    # Get fasta files to analyze CHANGE EXTENSION ACCORDINGLY!
    list_fasta_queries = glob(f"{args.files}{path.sep}*.fa")
    sys.stdout.write(f"{len(list_fasta_queries)} FASTA files to process\n")
    sys.stdout.flush()

    for fasta_file in list_fasta_queries:
        sys.stdout.write(f"\n#########\nProcessing file {path.basename(fasta_file)}...  \n")
        species_code = path.splitext(path.basename(fasta_file))[0].split("_")[0]

        sys.stdout.write("Cleaning input FASTA...  ")
        sys.stdout.flush()
        clean_fasta(fasta_file, workdir, species_code)
        fasta_file_cleaned = glob(f"{workdir}{path.sep}{species_code}_clean.fasta")[0]
        sys.stdout.write("done\n")
        sys.stdout.flush()

        #  Run HMMSEARCH
        sys.stdout.write("Performing HMSEARCH...  ")
        sys.stdout.flush()
        # path_hmm_models = f"{path.dirname(sys.argv[0])}{path.sep}data{path.sep}ltp_domains.hmm"
        # run_hmmsearch(fasta_query=fasta_file_cleaned, path_rep_out=hmmsearch_outdir, global_evalue=args.hmm_eval_glob,
        #               domain_evalue=args.hmm_eval_dom, hmm_models=path_hmm_models, hmm_cpus=args.threads)

        run_hmmsearch(fasta_query=fasta_file_cleaned, path_rep_out=hmmsearch_outdir, global_evalue=args.hmm_eval_glob,
                      domain_evalue=args.hmm_eval_dom, hmm_models=args.hmm, hmm_cpus=args.threads)
        sys.stdout.write("done\n")
        sys.stdout.flush()

        # Run motifSearch
        sys.stdout.write("Performing motifSearch...  ")
        sys.stdout.flush()
        run_motifsearch(fasta_query=fasta_file_cleaned, output_dir=ms_outdir, species_code=species_code)
        sys.stdout.write("done\n")
        sys.stdout.flush()
        sequences_to_get, df_summary_res = get_sequences_id_tp_extract(species_code=species_code,
                                                                       path_hmmer_res=hmmsearch_outdir,
                                                                       path_ms_res=ms_outdir)
        sys.stdout.write(f"{len(sequences_to_get)} sequences kept from HMMSEARCH & motifSearch analysis!\n")
        sys.stdout.flush()

        # Extract sequences from protein FASTA file
        sys.stdout.write("Extracting retained sequences for downstream analysis...  ")
        sys.stdout.flush()
        extract_seq(species_code=species_code, fasta_file=fasta_file_cleaned, path_out=path_fasta_seq,
                    list_ids=sequences_to_get, prefix_file="all_candidates")
        sys.stdout.write("done\n")
        sys.stdout.flush()
        path_retained_seq = glob(f"{path_fasta_seq}{path.sep}{species_code}*")[0]
        retained_seq = SeqIO.index(path_retained_seq, "fasta")

        sys.stdout.write("Running signalp...  ")
        sys.stdout.flush()

        # signalp5
        run_signalp(fasta_query=path_retained_seq, species_code=species_code, path_out=path_signalp_out)
        signalp_res = pds.read_csv(glob(f"{path_signalp_out}{path.sep}{species_code}*.signalp5")[0], sep="\t",
                                   comment="#", names=["seqName", "signalp_pred", "sp_score", "other_score", "cs_pos"])
        mature_prot = glob(f"{path_signalp_out}{path.sep}{species_code}_mature.fasta")[0]
        mature_prot_idx = SeqIO.index(mature_prot, "fasta")
        cysteine_count = {}
        for rec in SeqIO.parse(mature_prot, "fasta"):
            cysteine_count[rec.id] = count_cysteines(str(rec.seq))
        df_cysteines = pds.DataFrame.from_dict(cysteine_count, orient="index", columns=["cysteines_count"]).\
            reset_index()
        df_cysteines.rename(columns={"index": "seqName"}, inplace=True)
        df_cysteines["cysteines_count"] = df_cysteines["cysteines_count"].astype("int32")
        sys.stdout.write("done\n")
        sys.stdout.flush()

        # Calculate protein properties
        sys.stdout.write("Calculating protein properties...  ")
        sys.stdout.flush()

        # Prepare output tables
        df_summary_res = df_summary_res.merge(signalp_res[["seqName", "signalp_pred"]], on="seqName", how="left")
        df_summary_res = df_summary_res.merge(df_cysteines[["seqName", "cysteines_count"]], on="seqName", how="left")
        df_summary_res = df_summary_res.apply(get_prot_properties, args=[retained_seq, "wholeProt"], axis=1)
        df_summary_prot_w_sp = df_summary_res[df_summary_res["signalp_pred"] == "SP(Sec/SPI)"].copy()
        df_summary_prot_w_sp = df_summary_prot_w_sp.apply(get_prot_properties, args=[mature_prot_idx, "matureProt"],
                                                          axis=1)
        df_summary_res = pds.merge(df_summary_res,
                                   df_summary_prot_w_sp[["seqName", "isoelectric_matureProt", "gravy_matureProt",
                                                        "molecular_weight_matureProt"]],
                                   on="seqName", how="left")
        df_summary_res.fillna(value={"cysteines_count": "ND", "isoelectric_matureProt": "ND", "gravy_matureProt": "ND",
                                     "molecular_weight_matureProt": "ND"}, inplace=True)
        df_top_candidates = df_summary_res[(df_summary_res["hmmsearch"] == 1) & (df_summary_res["motifsearch"] == 1) &
                                           (df_summary_res["signalp_pred"] == "SP(Sec/SPI)") &
                                           (df_summary_res["cysteines_count"] == 8)]
        seq_ids_top_candidates = df_top_candidates["seqName"].to_list()
        extract_seq(species_code=species_code, fasta_file=fasta_file_cleaned, path_out=path_res_out,
                    list_ids=seq_ids_top_candidates, prefix_file="top_candidates_wholeProt")
        extract_seq(species_code=species_code, fasta_file=mature_prot, path_out=path_res_out,
                    list_ids=seq_ids_top_candidates, prefix_file="top_candidates_matureProt")
        path_seq_top_candidates = glob(f"{path_res_out}{path.sep}{species_code}*.fa")[0]

        sys.stdout.write("done\n")
        sys.stdout.flush()

        sys.stdout.write(f"{len(df_top_candidates)} top candidates identified!!\n")
        sys.stdout.flush()

        df_summary_res.to_csv(f"{path_res_out}{path.sep}{species_code}_all_results.tsv", sep="\t", index=False)
        df_top_candidates.to_csv(f"{path_res_out}{path.sep}{species_code}_top_candidates_results.tsv", sep="\t",
                                 index=False)

        # Run MEME
        sys.stdout.write("Running MEME...  ")
        sys.stdout.flush()
        run_meme(fasta_query=path_retained_seq, meme_res_out=path_meme_out, analysis_type="all_candidates",
                 species_code=species_code, cpus=args.threads)
        run_meme(fasta_query=path_seq_top_candidates, meme_res_out=path_meme_out, analysis_type="top_candidates",
                 species_code=species_code, cpus=args.threads)

        sys.stdout.write("done\n")
        sys.stdout.flush()

        sys.stdout.write(f"done processing {path.basename(fasta_file)}!\n#########\n")
        sys.stdout.flush()

    # Running MEME on merged top candidates from all query species
    sys.stdout.write("Running MEME on all top candidates...  ")
    sys.stdout.flush()
    top_candidates_files = glob(f"{path_res_out}{path.sep}*.fa")
    with open(f"{path_res_out}{path.sep}all_top_candidates.fasta", "wb") as outfile:
        for file in top_candidates_files:
            with open(file, "rb") as readfile:
                shutil.copyfileobj(readfile, outfile)
    outfile.close()
    all_top_candidates = path_res_out + path.sep + "all_top_candidates.fasta"
    run_meme(fasta_query=all_top_candidates, meme_res_out=path_meme_out, analysis_type="top_candidates",
             species_code="all", cpus=args.threads)
    sys.stdout.write("done\n")
    sys.stdout.flush()

    sys.stdout.write(f"Results location:\n=> HMMSEARCH: {hmmsearch_outdir}\n=> motifSearch: {ms_outdir}\n"
                     f"=> signalP: {path_signalp_out}\n=> MEME: {path_meme_out}\n=> Final results: {path_res_out}\n")
    sys.stdout.flush()


if __name__ == '__main__':
    main()
