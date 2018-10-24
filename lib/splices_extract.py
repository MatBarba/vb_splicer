#!env python3

import argparse
import logging
import os.path

# To parse a BAM file
import HTSeq
import sqlite3


class Splice():
    """Splice junction found in a read"""

    sql_fields = ('chrom', 'start', 'end', 'strand',
                  'left_flank', 'right_flank', 'coverage',)

    def __init__(self, chrom, start, end, strand,
                 left_flank=1, right_flank=1, coverage=1):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.left_flank = left_flank
        self.right_flank = right_flank
        self.coverage = coverage

    def __str__(self):
        return "%s:%d-%d (%s) [%d-%d]\n" % (
            self.chrom,
            self.start,
            self.end,
            self.strand,
            self.left_flank,
            self.right_flank
        )

    def from_aln(aln):
        """Create a list of splices from a read alignment"""
        chrom = aln.iv.chrom
        strand = aln.iv.strand
        read_start = aln.iv.start

        # Collect gaps spanned by the read
        # We keep a record of the matched sequence in seq,
        # and the gaps (potential splice junctions) in gap
        # There should be +1 seqs than gaps
        seqs = [0]
        gaps = []

        cur_seq = 0
        for c in aln.cigar:
            if c.type == "N":
                gaps.append(c.size)
                seqs.append(0)
                cur_seq += 1
            else:
                seqs[cur_seq] += c.size

        splices = []
        seq_diff = 0
        for position, gap in enumerate(gaps):
            left_flank = seqs[position]
            right_flank = seqs[position+1]
            seq_diff += left_flank
            start = read_start + seq_diff
            end = start + gap
            splices.append(
                Splice(chrom, start, end, strand, left_flank, right_flank)
            )

        return splices

    def sql_values(self):
        return [getattr(self, f) for f in self.sql_fields]

    def sql_placeholder(self):
        return ",".join(['?' for f in self.sql_fields])


def extract_splices(bam_input, sqlite_output):
    """Extract splices from a bam file, and write them in an SQLite db file"""
    logging.info("Read bam file " + bam_input)

    # Read bam file
    bam_reader = HTSeq.BAM_Reader(bam_input)

    count = 0
    count_splices = 0

    cur_splices = []
    cur_chrom = ''
    for aln in bam_reader:
        count += 1

        new_splices = Splice.from_aln(aln)
        if len(new_splices) > 0:
            # New chromosome? Store it!
            splice_chrom = new_splices[0].chrom
            if cur_chrom != splice_chrom:
                logging.debug("From %s to %s: Writing to %s" % (
                    cur_chrom, splice_chrom, sqlite_output))
                store_splices(cur_splices, sqlite_output)
                cur_splices = []
                cur_chrom = splice_chrom

            cur_splices = cur_splices + new_splices

        if count_splices >= 10:
            break
    store_splices(cur_splices, sqlite_output)

    logging.info("Total read alignments: " + str(count))
    logging.info("Total splices: " + str(count_splices))


def store_splices(splices, output):

    # Init sqlite db if new
    init_sqlite_db(output)

    conn = sqlite3.connect(output)
    c = conn.cursor()

    for splice in splices:
        values = splice.sql_values()
        fields = ",".join(splice.sql_fields)
        placeholder = splice.sql_placeholder()
        sql = "INSERT INTO splices("+fields+") VALUES("+placeholder+")"
        logging.debug(sql)
        logging.debug(values)
        c.execute(sql, values)
    conn.commit()
    conn.close()


def init_sqlite_db(output):
    if os.path.isfile(output):
        os.remove(output)

    conn = sqlite3.connect(output)
    c = conn.cursor()

    c.execute('''CREATE TABLE splices (
                chrom text,
                start int,
                end int,
                strand text,
                left_flank int,
                right_flank int,
                coverage int
              )''')


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
