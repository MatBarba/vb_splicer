#!env python3

import argparse
import logging

from Splice import Splice, SpliceCollection, SpliceDB
from SimpleGTF import SimpleGTF
from BCBio import GFF

import eHive
import os


class CreateGFF(eHive.BaseRunnable):

    def param_defaults(self):
        return {
            'category': 'known',
            'coverage': 1
        }

    def run(self):
        logging.basicConfig(level=logging.DEBUG)

        species = self.param_required('species')
        splice_db = self.param_required('splice_db')
        gtf_file = self.param_required('gtf_file')
        category = self.param_required('category')
        coverage = self.param_required('coverage')

        # Create the db file name
        filename = self.param_required('track_file')

        # Run it!
        CreateGFF.create_gff(splice_db, filename, category, coverage)

        self.dataflow({
            'species': species,
            'gff': filename,
            }, 2)

    def create_gff(input, filename, category, coverage=1):
        logging.info(
            "Import coverage filtered splices (coverage >= %d)" % coverage)
        input_db = SpliceDB(input)

        groups = {
            "all": [],
            "known": ["known"],
            "inbridge": ["inbridge"],
            "outbridge": ["outbridge"],
            "startends": ["left", "right"],
            "overlap": ["overlap"],
            "nocontact": ["nocontact"],
            "connected": ["inbridge", "outbridge", "left", "right"],
            "unconnected": ["overlap", "nocontact"],
            "duplicates": ["duplicates"],
        }
        antigroups = {
            "unknown": ["known", "duplicates"],
        }

        if category in groups:
            group = category
            tags = groups[category]

            logging.info("Writing group " + group)
            collection = input_db.get_collection(tags=tags, coverage=coverage)
            CreateGFF.print_gff(collection, filename)
        else:
            logging.warn("Unknown category %s" % category)

        if category in antigroups:
            group = category
            antitags = antigroups[category]

            logging.info("Writing group " + group)
            collection = input_db.get_collection(antitags=antitags, coverage=coverage)
            CreateGFF.print_gff(collection, filename)
        else:
            logging.warn("Unknown category %s" % category)

    def print_gff(collection, output):
        records = []
        for splice in collection.get_splices():
            records.append(splice.get_gff_record())

        logging.info("%d records to write in GFF %s" % (len(records), output))

        with open(output, 'w') as gff:
            GFF.write(records, gff)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Splice db to use')
    parser.add_argument('--known', dest='known',
                        help='output gff with only known splices')
    parser.add_argument(
        '--inbridge',
        dest='inbridge',
        help='output gff with splices that match different splices than known')
    parser.add_argument(
        '--outbridge',
        dest='outbridge',
        help='output gff with splices that match different genes')
    parser.add_argument(
        '--startends',
        dest='startends',
        help='output gff with splices that do not match any known introns')
    parser.add_argument(
        '--nocontact',
        dest='nocontact',
        help='output gff with all other splices')
    parser.add_argument(
        '--all',
        dest='all',
        help='output gff with all splices')
    parser.add_argument(
        '--unknown',
        dest='unknown',
        help='output gff with all splices that are not known')
    parser.add_argument(
        '--connected',
        dest='connected',
        help='output gff with all unknown splices that are connected to a known exon')
    parser.add_argument(
        '--unconnected',
        dest='unconnected',
        help='output gff with all unknown splices that are not connected to a known exon')
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

    outputs = {
        'all': args.all,
        'unknown': args.unknown,
        'known': args.known,
        'startends': args.startends,
        'inbridge': args.inbridge,
        'outbridge': args.outbridge,
        'nocontact': args.nocontact,
        'connected': args.connected,
        'unconnected': args.unconnected,
    }

    CreateGFF.create_gff(args.input, outputs, int(args.coverage))


if __name__ == "__main__":
    main()
