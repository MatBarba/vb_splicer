#!env python3

import argparse
import logging

import sqlite3
from Splice import SpliceCollection, SpliceDB

import os
import eHive


class MergeSplices(eHive.BaseRunnable):

    def run(self):

        logging.basicConfig(level=logging.INFO)

        splice_dir = self.param_required('splice_dir')
        species = self.param_required('species')
        splice_dbs = self.param_required('splice_dbs')
        force_merge = self.param('force_merge')

        # Create the db file name
        merged_db = os.path.join(splice_dir, species + '.sqlite')

        # Run it!
        if species in splice_dbs:
            if force_merge or not os.path.exists(merged_db):
                self.merge_dbs(merged_db, splice_dbs[species])
            else:
                logging.info("Splice db %s exists: reusing" % merged_db)
        else:
            self.warning("Can't find %s among the splice dbs species" % species)

        self.dataflow({
            'splice_db': merged_db
        }, 2)

    def merge_dbs(self, output_db, input_dbs):
        """Merge a list of SpliceDBs into one"""

        db = SpliceDB(output_db)
        db.init()
        logging.info("Created output db %s" % output_db)

        # First, get the list of chroms from all dbs
        chroms = get_chroms(input_dbs)
        logging.info("Total chroms: %d" % len(chroms.keys()))

        # If we have a lot of (small) chroms, iterating over every one of them
        # is not efficient so we can group all the small ones
        chrom_groups = group_chroms(chroms)
        logging.info("%d chrom groups to merge" % len(chrom_groups))

        size = 0
        n = 1

        dbs = {}
        for chrom_group in sorted(chrom_groups):
            logging.info(
                "%d Merging of %d groups (%s)"
                % (n, len(chrom_group), chrom_group[0]))
            col = SpliceCollection()
            for input_db in input_dbs:
                logging.debug("Merge db %s" % input_db)

                if input_db not in dbs:
                    dbs[input_db] = SpliceDB(input_db)
                in_db = dbs[input_db]
                in_col = in_db.get_collection(chroms=chrom_group)

                # Merge all splices
                col.add_splices(in_col.get_splices())
            logging.info("Load %d splices into db" % col.size)
            db.add_collection(col)
            size += col.size
            n += 1


def group_chroms(chroms):

    group_threshold = 50000

    small_groups = []
    groups = []
    small_groups_size = 0

    for chrom, num in chroms.items():
        if num > group_threshold:
            groups.append([chrom])
        else:
            small_groups.append(chrom)
            small_groups_size += num

            if small_groups_size >= group_threshold or len(small_groups) >= 500:
                groups.append(small_groups)
                small_groups = []
                small_groups_size = 0

    if len(small_groups) > 0:
        groups.append(small_groups)

    return groups


def get_chroms(dbs):

    chroms = {}
    for db in dbs:
        sql = "SELECT chrom, count(*) FROM splices group by chrom"
        conn = sqlite3.connect(db)
        c = conn.cursor()
        c.execute(sql)

        for row in c.fetchall():
            chrom, num = row[0], row[1]
            if chrom not in chroms:
                chroms[chrom] = num
            else:
                chroms[chrom] += num

    return chroms


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
    MergeSplices.merge_dbs(args.output, args.input_dbs)
