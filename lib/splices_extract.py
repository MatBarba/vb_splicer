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
        self.key = "|".join([self.chrom,
                             str(self.start),
                             str(self.end),
                             self.strand
                             ])

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

        # In some cases the read is not aligned, so skip it
        if aln.iv == None:
            return []

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

    def same_splice(self, other_splice):
        """Compare two splices: return True if they are at the same position"""

        if (self.chrom == other_splice.chrom and
                self.strand == other_splice.strand and
                self.start == other_splice.start and
                self.end == other_splice.end):
            return True
        else:
            return False

    def expand(self, new_splice):
        """Merge the splices if they are at the same position

        If they are not, return False
        """

        if self.same_splice(new_splice):
            self.merge(new_splice)
            return True
        else:
            return False

    def merge(self, new_splice):
        """Merge the coverage and flanks of two splices

        This should not be used directly, use expand instead to only merge is
        the splices are actually the same
        """

        # Merge coverages
        self.coverage += new_splice.coverage

        # Merge left flank
        if self.left_flank < new_splice.left_flank:
            self.left_flank = new_splice.left_flank

        # Merge right flank
        if self.right_flank < new_splice.right_flank:
            self.right_flank = new_splice.right_flank

    def to_gff(self):
        """Return a formatted string to be written in a gff file"""
        return ""


class SpliceCollection():
    """Collection of splices"""

    def __init__(self):
        self.splices = {}
        self.size = 0

    def add_splices(self, splices):
        for splice in splices:
            self.add_splice(splice)

    def add_splice(self, splice):
        if splice.key in self.splices:
            self.splices[splice.key].expand(splice)
        else:
            self.splices[splice.key] = splice
            self.size += 1

    def get_splices(self):
        return list(self.splices.values())


class SpliceDB():
    """Interface to an SQLite db where splices can be stored"""

    def __init__(self, path):
        self.path = path
        self.conn = None

    def drop_db(self):
        os.remove(self.path)

    def get_connection(self):
        if self.conn == None:
            self.conn = sqlite3.connect(self.path)
        return self.conn

    def init(self):
        if os.path.isfile(self.path):
            self.drop_db()

        conn = self.get_connection()
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

    def add_collection(self, collection):
        """Store all the Splices from a SpliceCollection in the SpliceDB"""

        conn = self.get_connection()
        c = conn.cursor()
        for splice in collection.get_splices():
            values = splice.sql_values()
            fields = ",".join(splice.sql_fields)
            placeholder = splice.sql_placeholder()
            sql = "INSERT INTO splices("+fields+") VALUES("+placeholder+")"
            logging.debug(sql)
            logging.debug(values)
            c.execute(sql, values)
        conn.commit()

    def get_collection(self, chrom=''):
        """Retrieve a SpliceCollection from the SpliceDB"""
        conn = sqlite3.connect(path)
        c = conn.cursor()

        sql = "SELECT " + ",".join(Splice.fields) + " FROM splices"
        data = []
        if chrom != '':
            sql += " WHERE chrom=?";
            data.append(chrom)
        c.execute(sql, data)

        col = SpliceCollection()
        for row in c.fetchall():
            s = {}
            for i, val in enumerate(row):
                field = Splice.fields[i]
                s[field] = val
            splice = Splice(
                s[chrom],
                s[start],
                s[end],
                s[strand],
                s[left_flank],
                s[right_flank],
                s[coverage]
            )
            col.add_splice(splice)

        return col


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
