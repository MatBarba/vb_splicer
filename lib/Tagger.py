#!env python3

import argparse
import logging

from Splice import Splice, SpliceCollection, SpliceDB, SpliceGeneCollection, SpliceTranscriptCollection
from SimpleGTF import SimpleGTF
from BCBio import GFF
from bx.intervals.intersection import Interval, IntervalTree

import eHive
import os

import requests, sys

class Tagger(eHive.BaseRunnable):
    para_cache = {}

    def param_default(self):
        return {
                rest_server: "https://www.vectorbase.org/rest",
                do_donor_acceptor: False,
                do_duplicates: False
        }

    def run(self):
        logging.basicConfig(level=logging.INFO)

        splice_db = self.param_required('splice_db')
        gtf_file = self.param_required('gtf_file')
        do_not_retag = self.param('do_not_retag')
        species = self.param('species')
        rest_server = self.param('rest_server')
        do_donor_acceptor = self.param('do_donor_acceptor')
        do_duplicates = self.param('do_duplicates')

        if do_not_retag:
            return

        if gtf_file is None:
            logging.warn("No GTF file, abort")
            return

        # Run it
        Tagger.tag_splices(splice_db, gtf_file, species, rest_server, do_donor_acceptor, do_duplicates)

    def tag_splices(input, gtf_path, species, rest_server, do_donor_acceptor, do_duplicates):
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
            'nocontact': 0
        }

        logging.info("Get chroms from db")
        chroms = input_db.get_chroms()
        
        # Prepare sequence cache
        cur_chrom = ''
        cur_start = 0
        cur_end = 0
        cur_seq = ''
        seq_cache_length = 1000000

        nc = 0

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
                if introns.is_known(splice):
                    same = introns.get_same_splice(splice)

                    stats['known'] += 1
                    splice.set_tag('known')
                    splice.set_gene(same.gene)
                    sg_collection.add_known_intron(same.gene, splice.coverage)
                    for tr in same.transcripts:
                        st_collection.add_known_intron(tr, splice.coverage)

                # Maybe the start/end of the splice touch a known exon?
                else:
                    # One side known from introns
                    left_ok = introns.left_is_known(splice)
                    right_ok = introns.right_is_known(splice)
                    lefts = introns.get_left_splices(splice)
                    rights = introns.get_right_splices(splice)
                    # Or just to the left/right of a transcript (terminal exons)
                    gene_to_the_right = st_collection.splice_to_the_left(splice)
                    gene_to_the_left = st_collection.splice_to_the_right(splice)
                    
                    # Determine which gene touches the left/right sides
                    left_gene = None
                    right_gene = None
                    if left_ok:
                        left_gene = lefts[0].gene
                    elif gene_to_the_left:
                        left_gene = gene_to_the_left
                    if right_ok:
                        right_gene = rights[0].gene
                    elif gene_to_the_right:
                        right_gene = gene_to_the_right
                    
                    # Categorization of unknown splice junctions
                    if left_gene and right_gene:
                        # Links exons in the same gene
                        if left_gene == right_gene:
                            stats['inbridge'] += 1
                            splice.set_tag('inbridge')
                            splice.set_gene(left_gene)
                            splice.set_gene2(right_gene)

                        # Links exons in different genes
                        else:
                            # Check for duplication!
                            if do_duplicates and Tagger.check_duplicates(rest_server, left_gene, right_gene):
                                stats['duplicates'] += 1
                                splice.set_tag('duplicates')
                            else:
                                stats['outbridge'] += 1
                                splice.set_tag('outbridge')
                            splice.set_gene(left_gene)
                            splice.set_gene2(right_gene)
                    elif left_gene:
                        stats['left'] += 1
                        splice.set_tag('left')
                        splice.set_gene(left_gene)
                    elif right_gene:
                        stats['right'] += 1
                        splice.set_tag('right')
                        splice.set_gene(right_gene)

                    # Just an overlap?
                    else:
                        gene1, gene2 = splice.gene_overlap(genes_intervals)

                        if gene1 and gene2:
                            splice.set_gene(gene1)
                            splice.set_gene2(gene2)
                            if do_duplicates and Tagger.check_duplicates(rest_server, gene1, gene2):
                                stats['duplicates'] += 1
                                splice.set_tag('duplicates')
                            else:
                                splice.set_tag('overlap')
                                stats['overlap'] += 1

                        elif gene1:
                            splice.set_gene(gene1)
                            splice.set_tag('overlap')
                            stats['overlap'] += 1
                        elif gene2:
                            splice.set_gene(gene2)
                            splice.set_tag('overlap')
                            stats['overlap'] += 1
                        else:
                            stats['nocontact'] += 1
                            splice.set_tag('nocontact')

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

    def check_duplicates(rest_server, gene1, gene2):
        # Check para cache
        if gene1 in Tagger.para_cache and gene2 in Tagger.para_cache[gene1]:
            return Tagger.para_cache[gene1][gene2]

        paralogues = Tagger.get_paralogues(rest_server, gene1)
        if paralogues is None:
            return False
        
        # If the both genes are close paralogs (with >50% similarity)
        is_duplicate = False
        for para in paralogues:
            target = para["target"]
            if target["id"] == gene2 and float(target["perc_pos"]) > 50:
                is_duplicate = True
                break

        # Keep in cache
        if gene1 in Tagger.para_cache:
            Tagger.para_cache[gene1][gene2] = is_duplicate
        else:
            Tagger.para_cache[gene1] = { gene2: is_duplicate }
            
        return is_duplicate

    def get_paralogues(rest_server, gene):
        ext = "/homology/id/%s?" % gene
        r = requests.get(rest_server + ext, headers={ "Content-Type" : "application/json"}, params={"type": "paralogues"})
        
        # Too bad if we fail, but not catastrophic
        if not r.ok:
            r.raise_for_status()
        #    sys.exit()
        
        homologues = r.json()
        if "data" in homologues and len(homologues["data"]) == 1:
            return homologues["data"][0]['homologies']
        else:
            return None


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
