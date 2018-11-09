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
            'links': gff_basename + '_links.gff',
            'others': gff_basename + '_others.gff',
        }

        # Run it!
        self.create_gff(splice_db, outputs, gtf_file, coverage)

        for gff in outputs.values():
            self.dataflow({
                'species': species,
                'gff': gff,
            }, 2)

    def create_gff(self, input, outputs, gtf_path, coverage=1):
        logging.info(
            "Import coverage filtered splices (coverage >= %d)" % coverage)
        input_db = SpliceDB(input)
        collection = input_db.get_collection()

        if collection.size == 0:
            logging.warn("No splices, abort")
            return

        # FILTER BY GENES
        if gtf_path is not None:
            filters = outputs.keys()
            cols = self.filter_splices(collection, gtf_path, filters, coverage)
            for output, col in cols.items():
                if output in outputs and outputs[output] is not None:
                    path = outputs[output]
                    self.print_gff(col, path)

        logging.info("Final number of splices: %d" % collection.size)

        if 'all' in outputs:
            self.print_gff(collection, outputs['all'])

    def filter_splices(self, collection, gtf_path, filters, coverage=1):

        introns = self.get_introns(gtf_path)

        outputs = {
            'known': SpliceCollection(),
            'startends': SpliceCollection(),
            'links': SpliceCollection(),
            'others': SpliceCollection(),
        }

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
                if 'known' in filters:
                    outputs['known'].add_splice(splice)
            elif left_ok and right_ok:
                stats['links'] += 1
                if 'links' in filters:
                    outputs['links'].add_splice(splice)
            elif left_ok:
                stats['left'] += 1
                if 'startends' in filters:
                    outputs['startends'].add_splice(splice)
            elif right_ok:
                stats['right'] += 1
                if 'startends' in filters:
                    outputs['startends'].add_splice(splice)
            else:
                stats['uncontacted'] += 1
                if 'others' in filters and splice.coverage >= coverage:
                    outputs['others'].add_splice(splice)

        logging.info("Splices filter stats:")
        for stat, num in stats.items():
            logging.info("%s : %d" % (stat, num))

        return outputs

    def get_introns(self, gtf_path):
        genes = self.get_genes(gtf_path)

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

    def get_genes(self, gtf_path):
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
    parser.add_argument('--known', dest='known',
                        help='output gff with only known splices')
    parser.add_argument(
        '--links',
        dest='links',
        help='output gff with splices that match different transcripts/genes')
    parser.add_argument(
        '--startends',
        dest='startends',
        help='output gff with splices that do not match any known introns')
    parser.add_argument(
        '--others',
        dest='others',
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
        'links': args.links,
        'others': args.others,
    }

    CreateGFF.create_gff(
        args.input, outputs, args.gtf, int(args.coverage))


if __name__ == "__main__":
    main()
