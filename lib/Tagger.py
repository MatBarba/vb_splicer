#!env python3

import argparse
import logging

from Splice import Splice, SpliceCollection, SpliceDB
from SimpleGTF import SimpleGTF
from BCBio import GFF

import eHive
import os


class Tagger(eHive.BaseRunnable):
    def run(self):
        splice_db = self.param_required('splice_db')
        gtf_file = self.param_required('gtf_file')

        if gtf_file is None:
            logging.warn("No GTF file, abort")
            return

        # Run it
        tag_splices(splice_db, gtf_file)

    def tag_splices(input, gtf_path):
        input_db = SpliceDB(input)
        collection = input_db.get_collection()

        if collection.size == 0:
            logging.warn("No splices, abort")
            return
        
        # Get data from the GTF
        introns = Tagger.get_introns(gtf_path)

        stats = {
            'known': 0,
            'inbridge': 0,
            'outbridge': 0,
            'left': 0,
            'right': 0,
            'nocontact': 0
        }

        # Compare all the splices with the introns
        for splice in collection.get_splices():
            left_ok = introns.left_is_known(splice)
            right_ok = introns.right_is_known(splice)
            lefts = introns.get_left_splices(splice)
            rights = introns.get_right_splices(splice)

            if introns.is_known(splice):
                same = introns.get_same_splice(splice)

                stats['known'] += 1
                splice.set_tag('known')
                splice.set_gene(same.gene)
            elif left_ok and right_ok:
                if lefts[0].gene == rights[0].gene:
                    stats['inbridge'] += 1
                    splice.set_tag('inbridge')
                    splice.set_gene(lefts[0].gene)
                else:
                    stats['outbridge'] += 1
                    splice.set_tag('outbridge')

            elif left_ok:
                stats['left'] += 1
                splice.set_tag('left')
                splice.set_gene(lefts[0].gene)
            elif right_ok:
                stats['right'] += 1
                splice.set_tag('right')
                splice.set_gene(rights[0].gene)
            else:
                stats['nocontact'] += 1
                splice.set_tag('nocontact')

        logging.info("Splices filter stats:")
        for stat, num in stats.items():
            logging.info("%s : %d" % (stat, num))

        # Store the tags back in the database
        input_db.tag_back(collection)

    def get_introns(gtf_path):
        genes = Tagger.get_genes(gtf_path)

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
                            intron = Splice(exon.chrom, start, end, '.', gene=gene.id)
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
