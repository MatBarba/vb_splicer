#!env python3

import argparse
import logging

# To parse a BAM file
import HTSeq


class Splice():
    """Splice junction found in by a read"""

    def __init__(self, chrom, start, end, strand, left_flank=1, right_flank=1):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.left_flank = left_flank
        self.right_flank = right_flank

    def __str__(self):
        return "%s:%d-%d (%s) [%d-%d]\n" % (
            self.chrom,
            self.start,
            self.end,
            self.strand,
            self.left_flank,
            self.right_flank
        )

    def from_read(read):
        """Create a list of splices from a read"""
        chrom = read.iv.chrom
        strand = read.iv.strand
        read_start = read.iv.start

        # Collect gaps spanned by the read
        # We keep a record of the matched sequence in seq,
        # and the gaps (potential splice junctions) in gap
        # There should be +1 seqs than gaps
        seqs = [0]
        gaps = []

        cur_seq = 0
        for c in read.cigar:
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


def extract_splices(bam_input, sqlite_output):
    """Extract splices from a bam file, and write them in an SQLite db file"""
    logging.info("Read bam file " + bam_input)

    # Read bam file
    bam_reader = HTSeq.BAM_Reader(bam_input)

    count = 0
    count_splices = 0

    cur_splices = []
    cur_chrom = ''
    for read in bam_reader:
        count += 1

        new_splices = Splice.from_read(read)
        if len(new_splices) > 0:
            # New chromosome? Store it!
            splice_chrom = new_splices[0].chrom
            if cur_chrom != splice_chrom:
                logging.info("Writing to " + sqlite_output)
                store_splices(cur_splices, sqlite_output)
                cur_splices = []

            cur_splices = cur_splices + new_splices

            for s in new_splices:
                count_splices += 1
                logging.info("Splice: " + str(s))

        if count_splices >= 10:
            break
    store_splices(cur_splices, sqlite_output)

    logging.info("Total reads: " + str(count))
    logging.info("Total splices: " + str(count_splices))


def store_splices(splices, output):
    with open(output, 'a') as out:
        for splice in splices:
            out.write(str(splice))

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
