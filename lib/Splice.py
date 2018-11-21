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
import statistics

class SpliceGene():
    """A simple gene representation for the spliceDB"""

    sql_fields = ('gene', 'chrom', 'start', 'end', 'strand',
                  'introns', 'covered_introns', 'splice_coverage')

    def __init__(self, gene, chrom, start, end, strand='.',
                 introns=0, covered_introns=0, splice_coverage=0):
        self.gene = gene
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.introns = introns
        self.covered_introns = covered_introns
        self.splice_coverage = splice_coverage
        self.key = gene
        self.min_coverage = 10
        self.coverages = []

    def sql_values(self):
        return [getattr(self, f) for f in self.sql_fields]

    def sql_placeholder(self):
        return ",".join(['?' for f in self.sql_fields])
    
    def add_known_intron(self, coverage):
        self.covered_introns += 1

        if self.covered_introns > self.introns:
            logging.warning("There are more covered introns (%d) than known introns (%d) in %s" % (self.covered_introns, self.introns, self.gene))

        cov = self.splice_coverage

        self.coverages.append(coverage)

        if cov == 0:
            self.splice_coverage = coverage
        else:
            self.splice_coverage = ((cov - 1) * self.covered_introns + coverage) / cov

    def compute_splice_coverage(self):
        mean = statistics.mean(self.coverage)
        
        # Is there any outrageous outliers?
        # I.e. with a very weak coverage and way below the mean and the others?
        # Check median

class SpliceGeneCollection():
    """Collection of splice genes"""

    def __init__(self, genes=[]):
        self.genes = []
        self.keys = {}
        self.size = 0
        self.add_genes(genes)

    def add_genes(self, genes):
        for gene in genes:
            self.add_gene(gene)

    def add_gene(self, gene):
        if self.is_known(gene):
            return

        else:
            # We're going to add 1 gene, use the previous size for the index
            i = self.size

            self.genes.append(gene)
            self.keys[gene.key] = i

            # Increment size only at the end
            self.size += 1

    def get_genes(self):
        return self.genes

    def get_gene(self, key):
        if key in self.keys:
            return self.genes[self.keys[key]]

    def is_known(self, gene):
        return gene.key in self.keys

    def from_GTF_genes(gtf_genes):
        splice_genes = []
        for ggene in gtf_genes:
            all_introns = {}
            for tr in ggene.transcripts:
                exons = tr.exons
                if len(exons) > 1:
                    for i, exon in enumerate(exons):
                        if i > 0:
                            intron = str(exons[i-1].start) + str(exons[i-1].end) + "-" + str(exons[i].start) + str(exons[i].end)
                            all_introns[intron] = 1

            splice_gene = SpliceGene(
                    ggene.id,
                    ggene.chrom,
                    ggene.start,
                    ggene.end,
                    ggene.strand,
                    introns=len(all_introns)
                    )
            splice_genes.append(splice_gene)
        sg_collection = SpliceGeneCollection()
        sg_collection.add_genes(splice_genes)
        return sg_collection

    def add_known_intron(self, gene_id, coverage):
        gene = self.get_gene(gene_id)
        gene.add_known_intron(coverage)


class Splice():
    """Splice junction found in a read"""

    sql_fields = ('chrom', 'start', 'end', 'strand',
                  'left_flank', 'right_flank', 'coverage', 'tag', 'gene')

    def __init__(self, chrom, start, end, strand='.',
                 left_flank=1, right_flank=1, coverage=1, tag=None, gene=None, id=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.left_flank = left_flank
        self.right_flank = right_flank
        self.coverage = coverage
        self.tag = tag
        self.id = id
        self.gene = gene

        # Generate keys
        self.key = "-".join([self.chrom,
                             str(self.start),
                             str(self.end)
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

        # Extract sequence length and gaps from CIGAR
        seqs, gaps = Splice.collect_gaps(aln)

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

    def collect_gaps(aln):
        """Collect gaps spanned by the read
        We keep a record of the matched sequence in seq,
        and the gaps (potential splice junctions) in gap
        There should be +1 seqs than gaps
        """
        seqs = [0]
        gaps = []

        cur_seq = 0
        for c in aln.cigar:
            if c.type == "N":
                gaps.append(c.size)
                seqs.append(0)
                cur_seq += 1
            elif c.type in ("M", "I", "S", "="):
                seqs[cur_seq] += c.size

        return seqs, gaps

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

    def set_tag(self, tag):
        self.tag = tag

    def set_gene(self, gene):
        self.gene = gene

    def is_in_gene(self, genes):
        for gene in genes:
            left = self.start - self.left_flank
            right = self.end + self.right_flank
            if (self.chrom == gene.chrom
                    and (left >= gene.start and left <= gene.end
                    or right >= gene.start and right <= gene.end)):
                return gene.id


class SpliceCollection():
    """Cllection of splices"""

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
            # We're going to add 1 splice, use the previous size for the index
            i = self.size

            self.splices.append(splice)
            self.keys[splice.key] = i

            # Init left/right keys if needed
            if splice.leftkey not in self.leftkeys:
                self.leftkeys[splice.leftkey] = []
            if splice.rightkey not in self.rightkeys:
                self.rightkeys[splice.rightkey] = []

            # Append splice index to left/right keys
            self.leftkeys[splice.leftkey].append(i)
            self.rightkeys[splice.rightkey].append(i)

            # Increment size only at the end
            self.size += 1

    def get_splices(self):
        return self.splices

    def is_known(self, splice):
        return splice.key in self.keys

    def left_is_known(self, splice):
        return splice.leftkey in self.leftkeys

    def right_is_known(self, splice):
        return splice.rightkey in self.rightkeys

    def get_same_splice(self, splice):
        if self.is_known(splice):
            return self.get_splice_by_key(splice.key)

    def get_splice_by_key(self, key):
        if key in self.keys:
            i = self.keys[key]
            return self.splices[i]

    def get_left_splices(self, splice):
        return self.get_splices_by_leftkey(splice.leftkey)

    def get_right_splices(self, splice):
        return self.get_splices_by_rightkey(splice.rightkey)

    def get_splices_by_leftkey(self, leftkey):
        splices = []
        if leftkey in self.leftkeys:
            indexes = self.leftkeys[leftkey]
            for i in indexes:
                splices.append(self.splices[i])
        return splices

    def get_splices_by_rightkey(self, rightkey):
        splices = []
        if rightkey in self.rightkeys:
            indexes = self.rightkeys[rightkey]
            for i in indexes:
                splices.append(self.splices[i])
        return splices
    
    def filter_by_gene_coverage(splices, genes, coverage=1):
        acceptable_ratio = 0.5

        if len(genes) == 0:
            return splices

        new_splices = []

        print("Filtering %d splices" % len(splices))

        for splice in splices:
            if splice.gene is not None:
                if splice.gene in genes:
                    vals = genes[splice.gene]

                    # Compare min gene value and this splice coverage
                    if splice.coverage > vals['min'] * acceptable_ratio:
                        new_splices.append(splice)
                else:
                    # This case covered_means that there is no intron in the current gene model
                    # We want to keep the splices that are in those genes in case there should be
                    # introns in their place
                    if splice.coverage >= coverage:
                        new_splices.append(splice)
        print("Filtering result: %d splices" % len(new_splices))

        return new_splices



class SpliceDB():
    """Interface to an SQLite db where splices can be stored"""

    def __init__(self, path):
        self.path = path
        self.conn = None

    def drop_db(self):
        os.remove(self.path)

    def get_connection(self):
        if self.conn is None:
            logging.debug("Connect to %s" % self.path)
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
                    coverage int,
                    tag text,
                    gene text
                )''')
        c.execute('''CREATE INDEX chrom_idx ON splices (chrom);''')
        c.execute('''CREATE INDEX tag_idx ON splices (tag);''')
        c.execute('''CREATE INDEX coverage_idx ON splices (coverage);''')
        c.execute('''CREATE INDEX splices_gene_idx ON splices (gene);''')

        c.execute('''CREATE TABLE genes (
                    gene text,
                    chrom text,
                    start int,
                    end int,
                    strand text,
                    introns int,
                    covered_introns int,
                    splice_coverage float
                )''')
        c.execute('''CREATE INDEX genes_gene_idx ON genes (gene);''')

    def add_genes(self, genes):
        """Store a list of genes in the SpliceDB"""

        conn = self.get_connection()
        c = conn.cursor()
        for gene in genes:
            values = gene.sql_values()
            fields = ",".join(gene.sql_fields)
            placeholder = gene.sql_placeholder()
            sql = "INSERT INTO genes("+fields+") VALUES("+placeholder+")"
            logging.debug(sql)
            logging.debug(values)
            c.execute(sql, values)
        conn.commit()

    def get_genes(self):
        """Retrieve a list of genes from the SpliceDB"""
        conn = self.get_connection()
        c = conn.cursor()

        get_fields = list(SpliceGene.sql_fields)
        get_fields.append('ROWID')

        data = []
        sql = "SELECT " + ",".join(get_fields) + " FROM genes"
        sql += " ORDER BY chrom, start, end, strand"

        logging.debug(sql)
        c.execute(sql, data)

        genes = []
        for row in c.fetchall():
            g = {}
            for i, val in enumerate(row):
                field = get_fields[i]
                s[field] = val
            gene = SpliceGene(
                s['gene'],
                s['chrom'],
                s['start'],
                s['end'],
                s['strand'],
                s['introns'],
                s['covered_introns'],
                s['splice_coverage'],
            )
            genes.append(gene)

        return genes

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

    def get_collection(self, chrom='', chroms=[], tags=[], coverage=1, nointron_coverage=1, genes={}):
        """Retrieve a SpliceCollection from the SpliceDB"""

        splices = self.get_splices(chrom, chroms, tags, coverage, nointron_coverage)

        splices = SpliceCollection.filter_by_gene_coverage(splices, genes, coverage)

        col = SpliceCollection(splices)

        return col

    def get_splices(self, chrom='', chroms=[], tags=[], coverage=1, nointron_coverage=1):
        """Retrieve a raw list of splices from the SpliceDB"""
        conn = self.get_connection()
        c = conn.cursor()

        get_fields = list(Splice.sql_fields)
        get_fields.append('ROWID')
        select_fields = map(lambda x: "s." + x, get_fields)

        sql = "SELECT " + ", ".join(select_fields) + " FROM splices s"
        conditions = []
        data = []
        left_join = ''
        if chrom != '':
            conditions.append("s.chrom=?")
            data.append(chrom)
        if len(chroms) > 0:
            chroms_cond = "s.chrom=?"
            chrom_condition = []
            for chrom in chroms:
                chrom_condition.append(chroms_cond)
                data.append(chrom)

            conditions.append("(" + " OR ".join(chrom_condition) + ")")
        if len(tags) > 0:
            tags_cond = "tag=?"
            tag_condition = []
            for tag in tags:
                tag_condition.append(tags_cond)
                data.append(tag)

            conditions.append("(" + " OR ".join(tag_condition) + ")")
        if coverage > 1:
            conditions.append("s.coverage>=?")
            data.append(coverage)
        if nointron_coverage > 1:
            left_join = "LEFT JOIN genes g USING(gene)"
            conditions.append("(g.covered_introns > 0 OR s.coverage > ?)")
            data.append(nointron_coverage)

        if left_join != '':
            sql += " " + left_join + " "

        if len(conditions) > 0:
            sql += " WHERE " + " AND ".join(conditions)

        sql += " ORDER BY s.chrom, s.start, s.end, s.strand"

        logging.debug(sql)
        c.execute(sql, data)

        splices = []
        for row in c.fetchall():
            s = {}
            for i, val in enumerate(row):
                field = get_fields[i]
                s[field] = val
            splice = Splice(
                s['chrom'],
                s['start'],
                s['end'],
                s['strand'],
                s['left_flank'],
                s['right_flank'],
                coverage=s['coverage'],
                tag=s['tag'],
                gene=s['gene'],
                id=s['ROWID']
            )
            splices.append(splice)

        return splices

    def detag(self):
        conn = self.get_connection()
        c = conn.cursor()

        c.execute('''UPDATE splices set tag=NULL, gene=NULL''')

    def tag_back(self, collection):
        """Apply the tags of splices back to database"""
        conn = self.get_connection()

        n = 1
        splices = collection.get_splices()
        for splice in splices:
#            logging.info("%d/%d %s" % (n, len(splices), str(splice).rstrip()))
            n += 1
            self.tag_splice(splice, conn.cursor())
        conn.commit()

    def tag_splice(self, splice, cursor):
        values = [splice.tag, splice.gene, splice.id]
        fields = ",".join(["tag", "gene", "ROWID"])
        conditions = "ROWID=?"
        sql = "UPDATE splices SET tag=?, gene=? WHERE " + conditions
        logging.debug(sql)
        logging.debug(values)
        cursor.execute(sql, values)

    def add_genes_collection(self, sg_collection):
        """Store the genes collection in the database"""
        conn = self.get_connection()
        c = conn.cursor()
        c.execute('''DELETE FROM genes''')

        for gene in sg_collection.get_genes():
            values = gene.sql_values()
            fields = ",".join(gene.sql_fields)
            placeholder = gene.sql_placeholder()
            sql = "INSERT INTO genes("+fields+") VALUES("+placeholder+")"
            logging.debug(sql)
            logging.debug(values)
            c.execute(sql, values)
        conn.commit()

    def get_genes_coverage(self):
        """Retrieve the splice coverage for all genes"""

        conn = self.get_connection()
        c = conn.cursor()

        sql = 'SELECT gene, count(*), min(coverage), sum(coverage) FROM splices WHERE tag="known" group by gene'
        c.execute(sql)

        genes = {}
        for row in c.fetchall():
            (gene, count, min_cov, sum_cov) = row
            genes[gene] = {
                    'min': min_cov,
                    'sum': sum_cov,
                    'count': count
                    }

        return genes

    def get_chroms(self):
        sql = "SELECT chrom, count(*) FROM splices group by chrom"
        conn = self.get_connection()
        c = conn.cursor()
        c.execute(sql)

        chroms = {}
        for row in c.fetchall():
            chrom, num = row[0], row[1]
            if chrom not in chroms:
                chroms[chrom] = num
            else:
                chroms[chrom] += num
        return chroms
