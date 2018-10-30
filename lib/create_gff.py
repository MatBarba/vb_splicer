#!env python3

import argparse
import logging

from Splice import Splice, SpliceCollection, SpliceDB
from SimpleGTF import SimpleGTF
from BCBio import GFF


def create_gff(input, output, gtf_path, filters=None, coverage=1):
    input_db = SpliceDB(input)
    collection = input_db.get_collection(coverage=coverage)

    if collection.size == 0:
        logging.warn("No splices, abort")
        return

    # FILTER BY GENES
    if filters is not None and gtf_path is not None:
        logging.warn("Filtering now")
        collection = filter_splices(collection, gtf_path, filters)

    logging.info("Final number of splices: %d" % collection.size)
    print_gff(collection, output)


def filter_known_splices(collection, gtf_path):
    genes = get_genes(gtf_path)

    # Find known introns
    known_introns = SpliceCollection()
    stats = {
        'gene': 0,
        'transcript': 0,
        'exon': 0,
        'intron': 0,
        'known': 0
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
                        intron = Splice(exon.chrom, start, end, '.')

                        if collection.is_known(intron):
                            stats['known'] += 1
                            db_splice = collection.get_same_splice(intron)
                            known_introns.add_splice(db_splice)
                    start = exon.end

    logging.info("Genes stats:")
    for stat, num in stats.items():
        logging.info("%s : %d" % (stat, num))

    return known_introns


def filter_splices(collection, gtf_path, filters):

    introns = get_introns(gtf_path)
    filtered = SpliceCollection()

    stats = {
        'known': 0,
        'links': 0,
        'left': 0,
        'right': 0,
        'links': 0,
        'uncontacted': 0
    }

    # Compare all the splices with the introns
    for splice in collection.get_splices():
        left_ok = introns.left_is_known(splice)
        right_ok = introns.right_is_known(splice)

        if introns.is_known(splice):
            stats['known'] += 1
            if filters == 'known':
                filtered.add_splice(splice)
        elif left_ok and right_ok:
            stats['links'] += 1
            if filters == 'links':
                filtered.add_splice(splice)
        elif left_ok:
            stats['left'] += 1
            if filters == 'startends':
                filtered.add_splice(splice)
        elif right_ok:
            stats['right'] += 1
            if filters == 'startends':
                filtered.add_splice(splice)
        else:
            stats['uncontacted'] += 1
            if filters == 'others':
                filtered.add_splice(splice)

    logging.info("Splices filter stats:")
    for stat, num in stats.items():
        logging.info("%s : %d" % (stat, num))

    return filtered


def get_introns(gtf_path):
    genes = get_genes(gtf_path)

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
                        intron = Splice(exon.chrom, start, end, '.')
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


def print_gff(collection, output):
    records = []
    for splice in collection.get_splices():
        records.append(splice.get_gff_record())

    logging.info("%d records to write in GFF" % len(records))

    with open(output, 'w') as gff:
        GFF.write(records, gff)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Path to output the created gff file")
    parser.add_argument('input', help='Splice db to use')
    parser.add_argument(
        '--gtf', dest='gtf', help='GFT file with genes features')
    parser.add_argument(
        '--filter', dest='filter', default='known',
        help='GTF filter: known, startends, links, others')
    parser.add_argument(
        '--coverage', dest='coverage', default=1, help='Minimum coverage')
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

    create_gff(
        args.input, args.output, args.gtf, args.filter, int(args.coverage))

if __name__ == "__main__":
    main()
