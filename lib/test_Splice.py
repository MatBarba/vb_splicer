import unittest
from HTSeq import SAM_Alignment
from Splice import Splice, SpliceCollection


class TestSplice(unittest.TestCase):
    base_sam = [
        "READ1",
        "16",
        "SEQ_REGION",
        "14059",
        "1",
        "38M",
        "*",
        "0",
        "0",
        "TTTCTAGAGTAACATTTTGGCGCGAATACAAGAGTGAT",
        "IGGIIGIIHGIIIIIIIIIIIIHIIIIIIIIIIIIIII",
    ]

    def test_splice_exceptions(self):
        """Test splice creation failures"""

        with self.assertRaises(TypeError):
            Splice()
        with self.assertRaises(TypeError):
            Splice('SEQ_REGION')
        with self.assertRaises(TypeError):
            Splice('SEQ_REGION', 1)
        with self.assertRaises(TypeError):
            Splice('SEQ_REGION', 1, 10)

        with self.assertRaises(TypeError):
            Splice(
                chrom='SEQ_REGION',
            )
        with self.assertRaises(TypeError):
            Splice(
                chrom='SEQ_REGION',
                start=1,
            )
        with self.assertRaises(TypeError):
            Splice(
                chrom='SEQ_REGION',
                start=1,
                end=10,
            )

    def test_splice_creation(self):
        """Test splice creation"""

        Splice('SEQ_REGION', 1, 10, '-')
        Splice('SEQ_REGION', 1, 10, '+')
        Splice('SEQ_REGION', 1, 10, '-', 2)
        Splice('SEQ_REGION', 1, 10, '-', 2, 3)

        Splice(
            chrom='SEQ_REGION',
            start=1,
            end=10,
            strand='-',
        )
        Splice(
            chrom='SEQ_REGION',
            start=1,
            end=10,
            strand='-',
            left_flank=2,
        )
        Splice(
            chrom='SEQ_REGION',
            start=1,
            end=10,
            strand='-',
            left_flank=2,
            right_flank=3
        )

    def test_str(self):
        """Test splice string"""

        splice = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice_text = "SEQ_REGION:1-10 (-) [2-3]\n"
        self.assertEqual(str(splice), splice_text)

    def test_key(self):
        splice1 = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice2 = Splice('SEQ_REGION', 1, 10, '+', 10, 3)

        self.assertEqual(splice1.key, "SEQ_REGION|1|10|-")
        self.assertEqual(splice2.key, "SEQ_REGION|1|10|+")

    def test_from_align(self):
        """Test splice creation from an actual SAM line (no gap)"""

        aln = SAM_Alignment.from_SAM_line("\t".join(self.base_sam))
        splices = Splice.from_aln(aln)
        self.assertEqual(splices, [])

    def test_one_gap(self):
        """Test splice creation with one gap"""

        custom_sam = self.base_sam
        custom_sam[5] = "12M100N26M"
        aln = SAM_Alignment.from_SAM_line("\t".join(custom_sam))
        splices = Splice.from_aln(aln)

        # We have one splice
        self.assertEqual(len(splices), 1)

        # Check the splice object
        if len(splices) == 1:
            splice_text = "SEQ_REGION:14070-14170 (-) [12-26]\n"
            self.assertEqual(str(splices[0]), splice_text)

    def test_two_gaps(self):
        """Test splice creation with two gaps"""

        custom_sam = self.base_sam
        custom_sam[5] = "12M100N12M200N14M"
        aln = SAM_Alignment.from_SAM_line("\t".join(custom_sam))
        splices = Splice.from_aln(aln)

        # We have two splices
        self.assertEqual(len(splices), 2)

        # Check the splice objects
        if len(splices) == 2:
            splice_text_1 = "SEQ_REGION:14070-14170 (-) [12-12]\n"
            splice_text_2 = "SEQ_REGION:14082-14282 (-) [12-14]\n"
            self.assertEqual(str(splices[0]), splice_text_1)
            self.assertEqual(str(splices[1]), splice_text_2)

    def test_same_splice(self):
        splice1 = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice2 = Splice('SEQ_REGION', 1, 10, '+', 2, 3)
        splice3 = Splice('SEQ_REGION_2', 1, 10, '-', 2, 3)
        splice4 = Splice('SEQ_REGION', 2, 10, '-', 2, 3)
        splice5 = Splice('SEQ_REGION', 1, 5, '-', 2, 3)

        self.assertEqual(splice1.same_splice(splice1), True)
        self.assertEqual(splice1.same_splice(splice2), False)
        self.assertEqual(splice1.same_splice(splice3), False)
        self.assertEqual(splice1.same_splice(splice4), False)
        self.assertEqual(splice1.same_splice(splice5), False)

    def test_expand(self):
        splice1 = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice2 = Splice('SEQ_REGION', 1, 10, '+', 2, 3)

        self.assertEqual(splice1.expand(splice2), False)
        self.assertEqual(splice1.coverage, 1)
        self.assertEqual(splice1.expand(splice1), True)
        self.assertEqual(splice1.coverage, 2)

    def test_merge(self):
        splice1 = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice2 = Splice('SEQ_REGION', 1, 10, '-', 10, 3)
        splice3 = Splice('SEQ_REGION', 1, 10, '-', 5, 20)

        splice1.expand(splice2)
        self.assertEqual(splice1.coverage, 2)
        self.assertEqual(splice1.left_flank, 10)
        self.assertEqual(splice1.right_flank, 3)

        splice1.expand(splice3)
        self.assertEqual(splice1.coverage, 3)
        self.assertEqual(splice1.left_flank, 10)
        self.assertEqual(splice1.right_flank, 20)


class TestSpliceCollection(unittest.TestCase):

    def test_splice_collection_creation(self):
        """Test splice collection creation"""

        col = SpliceCollection()
        self.assertEqual(col.size, 0)
        self.assertEqual(col.get_splices(), [])

        splice1 = Splice('SEQ_REGION', 1, 10, '-')
        splice2 = Splice('SEQ_REGION', 1, 10, '-', 5, 5)
        splice3 = Splice('SEQ_REGION', 10, 20, '+')

        col.add_splice(splice1)
        self.assertEqual(col.size, 1)
        col.add_splice(splice2)
        self.assertEqual(col.size, 1)
        col.add_splice(splice3)
        self.assertEqual(col.size, 2)


if __name__ == '__main__':
        unittest.main()
