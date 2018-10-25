#!env python3

import argparse
import logging

from Splice import SpliceDB
from BCBio import GFF


def create_gff(input, output):
    input_db = SpliceDB(input)
    col = input_db.get_collection()

    if col.size == 0:
        logging.warn("Collection is empty, abort")
        return

    records = []
    for splice in sorted(col.get_splices(), key=lambda x: x.key):
        records.append(splice.get_gff_record())

    logging.info("%d records to write in GFF" % len(records))

    with open(output, 'w') as gff:
        GFF.write(records, gff)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Path to output the created gff file")
    parser.add_argument('input', help='Splice db to use')
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

    create_gff(args.input, args.output)
