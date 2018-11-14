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
            'coverage': 1
        }

    def run(self):

        species = self.param_required('species')
        gff_dir = self.param_required('gff_dir')
        splice_db = self.param_required('splice_db')
        gtf_file = self.param_required('gtf_file')
        coverage = self.param('coverage')

        if not os.path.exists(gff_dir):
            os.makedirs(gff_dir)

        # Create the db file name
        gff_basename = os.path.join(gff_dir, species)
        outputs = {
            'known': gff_basename + '_known.gff',
            'startends': gff_basename + '_startends.gff',
            'inbridge': gff_basename + '_inbridge.gff',
            'outbridge': gff_basename + '_outbridge.gff',
            'nocontact': gff_basename + '_noncontact.gff',
        }

        # Run it!
        create_gff(splice_db, outputs, gtf_file, coverage)

        for gff in outputs.values():
            self.dataflow({
                'species': species,
                'gff': gff,
            }, 2)

    def create_gff(input, outputs, coverage=1):
        logging.info(
            "Import coverage filtered splices (coverage >= %d)" % coverage)
        input_db = SpliceDB(input)

        whole_collection = input_db.get_collection()
        if 'all' in outputs and outputs['all'] is not None:
            CreateGFF.print_gff(whole_collection, outputs['all'])

        groups = {
            "known": ["known"],
            "inbridge": ["inbridge"],
            "outbridge": ["outbridge"],
            "startends": ["left", "right"],
            "nocontact": ["nocontact"],
        }
        for group, tags in groups.items():
            if group in outputs and outputs[group] is not None:
                logging.info("Writing group " + group)
                collection = input_db.get_collection(tags=tags)
                CreateGFF.print_gff(collection, outputs[group])

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
        'known': args.known,
        'startends': args.startends,
        'inbridge': args.inbridge,
        'outbridge': args.outbridge,
        'nocontact': args.nocontact,
    }

    CreateGFF.create_gff(args.input, outputs, int(args.coverage))


if __name__ == "__main__":
    main()
