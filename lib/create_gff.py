#!env python3

import argparse
import logging

from Splice import SpliceDB
from SimpleGTF import SimpleGTF
from BCBio import GFF


def create_gff(input, output, gtf_path, coverage):
    input_db = SpliceDB(input)
    splices = input_db.get_splices(coverage=coverage)

    if len(splices) == 0:
        logging.warn("No splices, abort")
        return

    # FILTER BY GENES
    # genes = get_genes(gtf_path)

    print_gff(splices, output)


def get_genes(gtf_path):

    # First import the genes
    with open(gtf_path) as gtf:
        genes = SimpleGTF.parse_GTF(gtf.readlines())
    return genes


def print_gff(splices, output):
    records = []
    for splice in splices:
        records.append(splice.get_gff_record())

    logging.info("%d records to write in GFF" % len(records))

    with open(output, 'w') as gff:
        GFF.write(records, gff)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Path to output the created gff file")
    parser.add_argument('input', help='Splice db to use')
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

    create_gff(args.input, args.output, int(args.coverage))


if __name__ == "__main__":
    main()
