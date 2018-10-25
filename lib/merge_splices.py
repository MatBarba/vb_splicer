#!env python3

import argparse
import logging

import sqlite3
from Splice import SpliceCollection, SpliceDB


def merge_dbs(output_db, input_dbs):
    """Merge a list of SpliceDBs into one"""

    db = SpliceDB(output_db)
    db.init()
    logging.info("Created output db %s" % output_db)

    # First, get the list of chroms from all dbs
    chroms = get_chroms(input_dbs)
    logging.info("Total chroms: %d" % len(chroms))

    size = 0
    for chrom in sorted(chroms):
        col = SpliceCollection()
        for input_db in input_dbs:
            logging.debug("Merge db %s" % input_db)
            input_col = SpliceDB(input_db)
            in_col = input_col.get_collection(chrom=chrom)

            # Merge all splices
            col.add_splices(in_col.get_splices())
        db.add_collection(col)
        size += col.size
        logging.info(
            "Completed merging of %s (%d total merged splices)" % (chrom, size))


def get_chroms(dbs):

    chroms = {}
    for db in dbs:
        sql = "SELECT DISTINCT(chrom) FROM splices"
        conn = sqlite3.connect(db)
        c = conn.cursor()
        c.execute(sql)

        for row in c.fetchall():
            chrom = row[0]
            if chrom not in chroms:
                chroms[chrom] = 1

    return list(chroms.keys())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Path to output the merged db")
    parser.add_argument('input_dbs', metavar='input_dbs', nargs='+',
                        help='dbs to merge')
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

    # print(args.input_dbs)
    merge_dbs(args.output, args.input_dbs)
