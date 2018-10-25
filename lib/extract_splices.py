#!env python3

import argparse
import logging

import HTSeq
from Splice import Splice, SpliceCollection, SpliceDB


def extract_splices(bam_input, sqlite_output):
    """Extract splices from a bam file, and write them in an SQLite db file"""
    logging.info("Read bam file " + bam_input)

    # Read bam file
    bam_reader = HTSeq.BAM_Reader(bam_input)

    collection = SpliceCollection()
    db = SpliceDB(sqlite_output)
    db.init()

    count = 0
    count_splices = 0

    cur_chrom = ''
    for aln in bam_reader:
        count += 1

        new_splices = Splice.from_aln(aln)
        if len(new_splices) > 0:
            count_splices += len(new_splices)
            # New chromosome? Store it!
            splice_chrom = new_splices[0].chrom
            if cur_chrom != splice_chrom:
                if cur_chrom != '':
                    logging.info("%s done: Writing %d splices to %s" % (
                        cur_chrom, collection.size, sqlite_output))
                    db.add_collection(collection)
                collection = SpliceCollection()
                cur_chrom = splice_chrom

            collection.add_splices(new_splices)

    logging.info("%s done: Writing %d splices to %s" % (
        cur_chrom, collection.size, sqlite_output))
    db.add_collection(collection)

    logging.info("Total read alignments: " + str(count))
    logging.info("Total splices: " + str(count_splices))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to input bam file")
    parser.add_argument("output", help="Path to output a gff file")
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

    extract_splices(args.input, args.output)
