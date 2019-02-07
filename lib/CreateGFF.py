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
        logging.basicConfig(level=logging.DEBUG)

        species = self.param_required('species')
        gff_dir = self.param_required('gff_dir')
        splice_db = self.param_required('splice_db')
        gtf_file = self.param_required('gtf_file')

        if not os.path.exists(gff_dir):
            os.makedirs(gff_dir)

        # Create the db file name
        gff_basename = os.path.join(gff_dir, species)
        outputs = [
                { "category": 'known', "file": gff_basename + '_known.gff', "coverage": 1 },
                { "category": 'unknown', "file": gff_basename + '_unknown_all.gff', "coverage": 1 },
                { "category": 'unknown', "file": gff_basename + '_unknown_10.gff', "coverage": 10 },
                { "category": 'unknown', "file": gff_basename + '_unknown_100.gff', "coverage": 100 },
                { "category": 'unknown', "file": gff_basename + '_unknown_1000.gff', "coverage": 1000 },
                { "category": 'unknown', "file": gff_basename + '_unknown_10000.gff', "coverage": 10000 }
        ]

        # Run it!
        CreateGFF.create_gff(splice_db, outputs)

        for output in outputs:
            self.dataflow({
                'species': species,
                'gff': output["file"],
            }, 2)

    def create_gff(input, outputs, coverage=1):
        logging.info(
            "Import coverage filtered splices (coverage >= %d)" % coverage)
        input_db = SpliceDB(input)

        logging.info("Loade genes coverage")
        genes = input_db.get_genes_coverage()

        groups = {
            "all": [],
            "known": ["known"],
            "unknown": ["inbridge", "outbridge", "left", "right", "ingene", "outgene", "nocontact"],
            "inbridge": ["inbridge"],
            "outbridge": ["outbridge"],
            "startends": ["left", "right"],
            "ingene": ["ingene"],
            "outgene": ["outgene"],
            "nocontact": ["nocontact"],
            "connected": ["inbridge", "outbridge", "left", "right"],
            "unconnected": ["ingene", "outgene", "nocontact"],
        }

        for output in outputs:
            cat = output["category"]
            
            if cat in groups:
                group = cat
                tags = groups[cat]

                logging.info("Writing group " + group)
                filter_genes = {}
                filter_coverage = 1
                nointron_coverage = 1
                if group in ("startends", "ingene", "outgene", "connected"):
                    filter_genes = genes
                if group in ("ingene", "outgene", "unconnected"):
                    nointron_coverage = coverage
                if group in ("nocontact"):
                    filter_coverage = coverage
                if group in ("unknown"):
                    filter_coverage = output["coverage"]

                collection = input_db.get_collection(
                        tags=tags,
                        genes=filter_genes,
                        coverage=filter_coverage,
                        nointron_coverage=nointron_coverage
                        )

                if group in ("nocontact", "ingene", "outgene"):
                    collection.filter_by_overlap()
                
                CreateGFF.print_gff(collection, output["file"])

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
