#!env python3

import logging
import os.path

# To parse a BAM file
import sqlite3

# To create a GFF
# from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


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
        self.key = "-".join([self.chrom,
                             str(self.start),
                             str(self.end),
                             self.strand
                             ])
        self.leftkey = "-".join([self.chrom,
                                 str(self.start)
                                 ])
        self.rightkey = "-".join([self.chrom,
                                  str(self.end)
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

    def from_aln(aln, stranded=False):
        """Create a list of splices from a read alignment"""

        # In some cases the read is not aligned, so skip it
        if aln.iv is None:
            return []

        chrom = aln.iv.chrom
        if stranded:
            strand = aln.iv.strand
        else:
            strand = '.'
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

    def get_gff_record(self):
        """Return a formatted string to be written in a gff file"""

        gene_start = self.start - self.left_flank
        gene_end = self.end + self.right_flank
        strand = self.strand
        if (strand == '+'):
            strand = +1
        elif (strand == '-'):
            strand = -1
        else:
            strand = 0

        location = FeatureLocation(gene_start - 1, gene_end)
        gene_qualif = {
            "source": "rnaseq",
            "score": self.coverage,
            "Name": "x%d" % self.coverage,
            "ID": self.key
        }
        top_feat = SeqFeature(
            location,
            type="gene",
            strand=strand,
            qualifiers=gene_qualif
        )

        # Append the features
        sub_features = []
        sub_features.append(SeqFeature(
            FeatureLocation(gene_start - 1, self.start),
            'exon',
            strand=strand,
            qualifiers={"source": "rnaseq"}
        ))
#        sub_features.append(SeqFeature(
#            FeatureLocation(self.start, self.end),
#            'intron',
#            strand=strand,
#            qualifiers={"source": "rnaseq"}
#        ))
        sub_features.append(SeqFeature(
            FeatureLocation(self.end, gene_end),
            'exon',
            strand=strand,
            qualifiers={"source": "rnaseq"}
        ))
        top_feat.sub_features = sub_features

        # Create the actual record
        rec = SeqRecord(Seq(""), self.chrom)
        rec.features = [top_feat]

        return rec


class SpliceCollection():
    """Collection of splices"""

    def __init__(self, splices=[]):
        self.splices = []
        self.keys = {}
        self.leftkeys = {}
        self.rightkeys = {}
        self.size = 0
        self.add_splices(splices)

    def add_splices(self, splices):
        for splice in splices:
            self.add_splice(splice)

    def add_splice(self, splice):
        if self.is_known(splice):
            i = self.keys[splice.key]
            self.splices[i].expand(splice)

        else:
            self.splices.append(splice)
            i = self.size
            self.keys[splice.key] = i
            self.size += 1

            # Increase left key
            if splice.leftkey in self.leftkeys:
                self.leftkeys[splice.leftkey].append(i)
            else:
                self.leftkeys[splice.leftkey] = [i]

            # Increase right key
            if splice.rightkey in self.rightkeys:
                self.rightkeys[splice.rightkey].append(i)
            else:
                self.rightkeys[splice.rightkey] = [i]

    def get_splices(self):
        return self.splices

    def is_known(self, splice):
        return splice.key in self.keys

    def get_same_splice(self, splice):
        if self.is_known(splice):
            return self.get_splice_by_key(splice.key)

    def get_splice_by_key(self, key):
        if key in self.keys:
            i = self.keys[key]
            return self.splices[i]

    def get_splices_by_leftkey(self, leftkey):
        splices = []
        if leftkey in self.leftkeys:
            indexes = self.leftkeys[leftkey]
            for i in indexes:
                splices.append(splices[i])
        return splices

    def get_splices_by_rightkey(self, rightkey):
        splices = []
        if rightkey in self.rightkeys:
            indexes = self.rightkeys[rightkey]
            for i in indexes:
                splices.append(splices[i])
        return splices


class SpliceDB():
    """Interface to an SQLite db where splices can be stored"""

    def __init__(self, path):
        self.path = path
        self.conn = None

    def drop_db(self):
        os.remove(self.path)

    def get_connection(self):
        if self.conn is None:
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

    def get_collection(self, chrom='', coverage=1):
        """Retrieve a SpliceCollection from the SpliceDB"""

        splices = self.get_splices(chrom, coverage)

        col = SpliceCollection(splices)

        return col

    def get_splices(self, chrom='', coverage=1):
        """Retrieve a raw list of splices from the SpliceDB"""
        conn = self.get_connection()
        c = conn.cursor()

        sql = "SELECT " + ",".join(Splice.sql_fields) + " FROM splices"
        sql += " ORDER BY chrom, start, end, strand"
        conditions = []
        data = []
        if chrom != '':
            conditions.append("chrom=?")
            data.append(chrom)
        if coverage is not None and coverage > 1:
            conditions.append("coverage>=?")
            data.append(coverage)

        if len(conditions) > 0:
            sql += " WHERE " + " AND ".join(conditions)

        logging.debug(sql)
        c.execute(sql, data)

        splices = []
        for row in c.fetchall():
            s = {}
            for i, val in enumerate(row):
                field = Splice.sql_fields[i]
                s[field] = val
            splice = Splice(
                s['chrom'],
                s['start'],
                s['end'],
                s['strand'],
                s['left_flank'],
                s['right_flank'],
                s['coverage']
            )
            splices.append(splice)

        return splices
