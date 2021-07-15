# -*- coding: utf-8 -*-
# @Time    : 13/07/2021 11:29
# @Author  : Jean Keller
# @Email   : kellerjeanphd@gmail.com
# @File    : motif_finder.py
# @Software: PyCharm

import re
from Bio import SeqIO


class GetMotifs:
    def __init__(self, motif, fasta_file, path_count, fseq, fpos):
        self.mot = motif
        self.fseq = fasta_file
        self.fcount = path_count
        self.outseq = fseq
        self.outpos = fpos

    def checkinput(self):
        aa_ok = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
                 "M", "F", "P", "S", "T", "W", "Y", "V", "X"]
        GetMotifs.completter(self, aa_ok)

    def completter(self, lref):
        valid_mot = True
        for let in str(self.mot):
            if str(let).upper() not in lref and str(let).isalpha():
                valid_mot = False
                raise ValueError("Invalid residue {0} in motif {1}".format(let, self.mot))

    def replaceN(self):
        return self.mot.replace("X", ".")

    def searchmotif(self, mot2search, seq):
        hits = mot2search.finditer(str(seq))
        n = 0
        mot_spec = []
        for hsp in hits:
            n += 1
            mot_spec.append((hsp.start(), hsp.end(), hsp.group()))
        mot_spec.append(n)
        return mot_spec

    def writemotif(self):
        mot2search = re.compile(GetMotifs.replaceN(self))
        seqidx = 0
        for rec in SeqIO.parse(self.fseq, "fasta"):
            res = GetMotifs.searchmotif(self, mot2search, str(rec.seq))
            n_fwd, c_seqlen = 0, len(str(rec.seq))
            seqidx += 1
            if len(res) > 1:
                n_fwd = res[-1]
                for i in res[:-1]:
                    pstart, pend, seq = i[0] + 1, i[1], i[2]
                    GetMotifs.writepos(self.outpos, rec.id, self.mot, seq, pstart, pend)
                    GetMotifs.writeseq(self.outseq, rec.id, self.mot, pstart, pend, seq)
            GetMotifs.writecount(self.fcount, rec.id, self.mot, str(n_fwd), "seq"+str(seqidx))

    @staticmethod
    def writeseq(fseq, idseq, regexpr, pos_start, pos_end, seq2write):
        fseq.write(f">{idseq}_{pos_start}:{pos_end} {regexpr}\n{seq2write}\n")

    @staticmethod
    def writepos(fpos, idseq, regexpr, found_comb, pos_start, pos_end):
        fpos.write(f"{idseq}\t{regexpr}\t{found_comb}\t{pos_start}\t{pos_end}\n")

    @staticmethod
    def writecount(fcount, idseq, regexpr, count_fwd, seqidx):
        fcount.write(f"{seqidx}\t{idseq}\t{regexpr}\t{count_fwd}\n")
