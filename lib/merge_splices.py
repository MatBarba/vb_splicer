#!env python3

import argparse
import logging

from Splice import SpliceCollection, SpliceDB


def merge_dbs(output_db, input_dbs):
    """Merge a list of SpliceDBs into one"""

    col = SpliceCollection()
    db = SpliceDB(output_db)
    db.init()
    logging.info("Created output db %s" % output_db)

    for input_db in input_dbs:
        logging.info("Merge db %s" % input_db)
        input_col = SpliceDB(input_db)
        in_col = input_col.get_collection()

        # Merge all splices
        col.add_splices(in_col.get_splices())

    db.add_collection(col)


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
