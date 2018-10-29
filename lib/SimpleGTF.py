#!env python3

import logging
import re


class Feature():
    """A simple genomic feature"""
    def __init__(self, id, chrom, start, end, strand):
        self.id = id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand


class Exon(Feature):
    """A simple exon object"""


class Transcript(Feature):
    """A simple transcript object"""

    def __init__(self, id, chrom, start, end, strand):
        super().__init__(id, chrom, start, end, strand)
        self.exons = []

    def add_exon(self, exon):
        self.exons.append(exon)


class Gene(Feature):
    """A simple gene object"""

    def __init__(self, id, chrom, start, end, strand):
        super().__init__(id, chrom, start, end, strand)
        self.transcripts = []

    def add_transcript(self, transcript):
        self.transcripts.append(transcript)


class SimpleGTF():
    """A simplified GTF parser"""

    @classmethod
    def parse_GTF(self, lines):

        genes = []
        cur_gene, cur_transcript, cur_exon = None, None, None

        for line in lines:
            (chrom, source, type, start, end,
             score, strand, phase, group) = line.split("\t")

            names = SimpleGTF.get_names(group)
            if len(names) == 0:
                continue

            if type == 'gene':
                cur_gene = Gene(names['gene_id'],
                                chrom, start, end, strand)
                genes.append(cur_gene)
            elif type == 'transcript':
                if 'transcript_id' not in names:
                    logging.warning("No transcript_id found in %s" % group)
                cur_transcript = Transcript(names['transcript_id'],
                                            chrom, start, end, strand)
                cur_gene.add_transcript(cur_transcript)
            elif type == 'exon':
                cur_exon = Exon(names['exon_number'],
                                chrom, start, end, strand)
                cur_transcript.add_exon(cur_exon)

        return genes

    @staticmethod
    def get_names(group):
        names = {}
        for field in group.split(';'):
            search = re.search('\s*([\w-]+) "([\w-]+)"\s*', field)
            if search is not None:
                type, name = search.group(1, 2)
                names[type] = name
        return names
