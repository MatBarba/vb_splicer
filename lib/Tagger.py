#!env python3

import argparse
import logging

from Splice import Splice, SpliceCollection, SpliceDB, SpliceGeneCollection, SpliceTranscriptCollection
from SimpleGTF import SimpleGTF
from BCBio import GFF
from bx.intervals.intersection import Interval, IntervalTree

import eHive
import os
import csv

import requests, sys

class Tagger(eHive.BaseRunnable):
    para_cache = {}

    def param_default(self):
        return {
                rest_server: "https://www.vectorbase.org/rest",
                do_donor_acceptor: False
        }

    def run(self):
        logging.basicConfig(level=logging.INFO)

        splice_db = self.param_required('splice_db')
        gtf_file = self.param_required('gtf_file')
        do_not_retag = self.param('do_not_retag')
        species = self.param('species')
        rest_server = self.param('rest_server')
        do_donor_acceptor = self.param('do_donor_acceptor')
        paralogs_dir = self.param('paralogs_dir')

        if do_not_retag:
            return

        if gtf_file is None:
            logging.warn("No GTF file, abort")
            return

        # Run it
        Tagger.tag_splices(splice_db, gtf_file, species, rest_server, do_donor_acceptor, paralogs_dir)

    def tag_splices(input, gtf_path, species, rest_server, do_donor_acceptor, paralogs_dir):
        input_db = SpliceDB(input)
        input_db.detag()

        # Get data from the GTF
        logging.info("Get introns from GTF")
        genes = Tagger.get_genes(gtf_path)
        genes_intervals = Tagger.make_gene_interval_tree(genes)
        introns = Tagger.get_introns(genes)
        sg_collection = SpliceGeneCollection.from_GTF_genes(genes)
        st_collection = SpliceTranscriptCollection.from_GTF_genes(genes)

        stats = {
            'known': 0,
            'inbridge': 0,
            'outbridge': 0,
            'duplicates': 0,
            'left': 0,
            'right': 0,
            'overlap': 0,
            'nocontact': 0,
            'outgene_left': 0,
            'outgene_right': 0,
            'ingene': 0,
        }

        logging.info("Get chroms from db")
        chroms = input_db.get_chroms()
        
        # Prepare sequence cache
        cur_chrom = ''
        cur_start = 0
        cur_end = 0
        cur_seq = ''
        seq_cache_length = 1000000

        nc = 1

        for chrom in chroms.keys():
            collection = input_db.get_collection(chroms=[chrom])

            # Compare all the splices with the introns
            splices = collection.get_splices()
            logging.info("%d/%d Tag %d splices in chrom %s" % (nc, len(chroms), len(splices), chrom))
            nc += 1
            for splice in splices:
                # Update sequence cache
                if do_donor_acceptor:
                    if chrom != cur_chrom or cur_end < splice.end:
                        cur_chrom = chrom
                        cur_start = splice.start
                        cur_end = cur_start + seq_cache_length
                        cur_seq = Tagger.get_rest_sequence(rest_server, species, cur_chrom, cur_start, cur_end)

                    # Get donor/acceptor (if there is a direction)
                    if splice.strand != '.':
                        donor, acceptor = Tagger.get_donor_acceptor(cur_seq, cur_start, splice)
                        splice.set_donor_acceptor(donor, acceptor)

                # Check against known intron
                tag = 'unknown'
                if introns.is_known(splice):
                    same = introns.get_same_splice(splice)
                    
                    tag = 'known'
                    splice.set_gene(same.gene)

                    # Record this intron in the genes and transcripts table
                    sg_collection.add_known_intron(same.gene, splice.coverage)
                    for tr in same.transcripts:
                        st_collection.add_known_intron(tr, splice.coverage)

                # Check overlaps etc. otherwise
                else:
                    # Do start and end overlap with some genes?
                    gene1, gene2 = splice.gene_overlap(genes_intervals)
                    
                    # Start or end overlap
                    if gene1 or gene2:

                        # Identify gene touching an intron from the left
                        left_gene = st_collection.splice_to_the_left(splice)
                        if introns.left_is_known(splice):
                            left_gene = introns.get_left_splices(splice)[0].gene
                        # Identify gene touching an intron from the right
                        right_gene = st_collection.splice_to_the_right(splice)
                        if introns.right_is_known(splice):
                            right_gene = introns.get_right_splices(splice)[0].gene

                        # Both overlap
                        if gene1 and gene2:
                            # Both start and end overlap the same gene
                            if gene1 == gene2:
                                tag = 'ingene'
                                splice.set_gene(gene1)

                                if left_gene is not None and right_gene is not None:
                                    tag = 'inbridge'
                                elif left_gene is not None:
                                    tag = 'left'
                                elif right_gene is not None:
                                    tag = 'right'

                            # The genes overlapped by start and end are different
                            else:
                                splice.set_gene(gene1)
                                splice.set_gene2(gene2)
                                tag = 'overlap'

                                # The overlap starts and ends at the end and starts of exons
                                if left_gene is not None and right_gene is not None:
                                    tag = 'outbridge'

                                # Check duplicates using paralogy
                                if paralogs_dir and Tagger.check_duplicates(paralogs_dir, species, gene1, gene2):
                                    tag = 'duplicates'

                        elif gene1:
                            splice.set_gene(gene1)
                            tag = 'outgene_right'
                        elif gene2:
                            splice.set_gene(gene2)
                            tag = 'outgene_left'
                    # No overlap
                    else:
                        tag = 'nocontact'

                stats[tag] += 1
                splice.set_tag(tag)

            # Store the tags back in the database
            input_db.tag_back(collection)

        input_db.add_genes_collection(sg_collection)
        input_db.add_transcripts_collection(st_collection)

        logging.info("Splices filter stats:")
        for stat, num in stats.items():
            logging.info("%s : %d" % (stat, num))

    def make_gene_interval_tree(genes):
        intervals = {}

        for gene in genes:
            if not gene.chrom in intervals:
                intervals[gene.chrom] = IntervalTree()

            intervals[gene.chrom].insert(gene.start, gene.end, gene)

        return intervals


    def get_rest_sequence(rest_server, species, chrom, start, end):
        web_species = species[0].upper() + species[1:]
        ext = "/sequence/region/%s/%s:%d-%d?" % (web_species, chrom, start, end)
        r = requests.get(rest_server + ext, headers={ "Content-Type" : "text/plain"})
        
        # Too bad if we fail, but not catastrophic
        if not r.ok:
            r.raise_for_status()
        #    sys.exit()

        return(r.text)

    def get_donor_acceptor(seq, seq_start, splice):
        donor_start = splice.start - seq_start + 1
        donor_end = donor_start + 2
        acceptor_start = splice.end - seq_start - 1
        acceptor_end = acceptor_start + 2

        donor = seq[donor_start:donor_end]
        acceptor = seq[acceptor_start:acceptor_end]

        if splice.strand == "-":
            new_donor = Tagger.reverse_strand(acceptor)
            new_acceptor = Tagger.reverse_strand(donor)
            donor = new_donor
            acceptor = new_acceptor

        return donor, acceptor

    def reverse_strand(seq):
        return seq.translate(str.maketrans('CGTA', 'GCAT'))[::-1]

    def check_duplicates(paralogs_dir, species, gene1, gene2):
        duplicates = {}

        # Check para cache
        if len(Tagger.para_cache) == 0:
            # Get gene1 paralogs
            paralogues = Tagger.get_species_paralogues(paralogs_dir, species)
            Tagger.para_cache = paralogues

        if gene1 in Tagger.para_cache and gene2 in Tagger.para_cache[gene1]:
            return True
        else:
            return False

    def get_species_paralogues(paralogs_dir, species):
        para_file = paralogs_dir + "/" + species + ".txt"
        if not os.path.exists(para_file):
            return {}
        
        paras = {}
        with open(para_file, 'r') as para_lines:
            for line in csv.reader(para_lines, delimiter="\t"):
                gene1, gene2 = line[0], line[1]
                if gene1 not in paras:
                    paras[gene1] = {}
                paras[gene1][gene2] = True
        return paras

    def get_introns(genes):

        # Find known introns
        col = SpliceCollection()

        stats = {
            'gene': 0,
            'transcript': 0,
            'exon': 0,
            'intron': 0,
        }

        for gene in genes:
            stats['gene'] += 1
            for tr in gene.transcripts:
                stats['transcript'] += 1
                exons = tr.exons
                if len(exons) > 0:
                    stats['exon'] += 1
                    start = 0
                    end = 0
                    for exon in sorted(exons, key=lambda x: x.start):
                        if start > 0:
                            stats['intron'] += 1
                            end = exon.start
                            intron = Splice(
                                    exon.chrom,
                                    start,
                                    end,
                                    exon.strand,
                                    gene=gene.id,
                                    transcripts=[tr.id]
                                    )
                            col.add_splice(intron)
                        start = exon.end

        logging.info("Genes stats:")
        for stat, num in stats.items():
            logging.info("%s : %d" % (stat, num))

        return col

    def get_genes(gtf_path):
        if gtf_path is None:
            logging.info("No gtf provided: no gene data used")
            return []

        # First import the genes
        logging.info("Import genes")
        with open(gtf_path) as gtf:
            genes = SimpleGTF.parse_GTF(gtf.readlines())
        return genes

    def print_gff(self, collection, output):
        records = []
        for splice in collection.get_splices():
            records.append(splice.get_gff_record())

        logging.info("%d records to write in GFF %s" % (len(records), output))

        with open(output, 'w') as gff:
            GFF.write(records, gff)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Splice db to use')
    parser.add_argument(
        '--gtf', dest='gtf', help='GFT file with genes features')
    parser.add_argument(
            '-d', '--debug',
            help="Print lots of logging.debugging statements",
            action="store_const", dest="loglevel", const=logging.DEBUG,
            default=logging.WARNING,
            )
    parser.add_argument(
            '-v', '--verbose',
            help="Be verbose",
            action="store_const", dest="loglevel", const=logging.INFO,
            )
    args = parser.parse_args()

    logging.basicConfig(level=args.loglevel)

    Tagger.tag_splices(args.input, args.gtf)


if __name__ == "__main__":
    main()
