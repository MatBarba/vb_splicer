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
            Splice(
                chrom='SEQ_REGION',
            )
        with self.assertRaises(TypeError):
            Splice(
                chrom='SEQ_REGION',
                start=1,
            )

    def test_splice_creation(self):
        """Test splice creation"""

        Splice('SEQ_REGION', 1, 10)
        Splice('SEQ_REGION', 1, 10, '-')
        Splice('SEQ_REGION', 1, 10, '+')
        Splice('SEQ_REGION', 1, 10, '-', 2)
        Splice('SEQ_REGION', 1, 10, '-', 2, 3)

        Splice(
            chrom='SEQ_REGION',
            start=1,
            end=10,
        )
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
        splice3 = Splice('SEQ_REGION', 1, 10, '.', 10, 3)

        self.assertEqual(splice1.key, "SEQ_REGION-1-10")
        self.assertEqual(splice2.key, "SEQ_REGION-1-10")
        self.assertEqual(splice3.key, "SEQ_REGION-1-10")

    def test_leftkey(self):
        splice1 = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice2 = Splice('SEQ_REGION', 1, 10, '+', 10, 3)
        splice3 = Splice('SEQ_REGION', 2, 10, '.', 10, 3)

        self.assertEqual(splice1.leftkey, "SEQ_REGION-1")
        self.assertEqual(splice2.leftkey, "SEQ_REGION-1")
        self.assertEqual(splice3.leftkey, "SEQ_REGION-2")

    def test_rightkey(self):
        splice1 = Splice('SEQ_REGION', 1, 10, '-', 2, 3)
        splice2 = Splice('SEQ_REGION', 1, 15, '+', 10, 3)

        self.assertEqual(splice1.rightkey, "SEQ_REGION-10")
        self.assertEqual(splice2.rightkey, "SEQ_REGION-15")

    def test_from_align(self):
        """Test splice creation from an actual SAM line (no gap)"""

        custom_sam = self.base_sam
        custom_sam[5] = "38M"
        aln = SAM_Alignment.from_SAM_line("\t".join(custom_sam))
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
            splice_text = "SEQ_REGION:14070-14170 (.) [12-26]\n"
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
            splice_text_1 = "SEQ_REGION:14070-14170 (.) [12-12]\n"
            splice_text_2 = "SEQ_REGION:14082-14282 (.) [12-14]\n"
            self.assertEqual(str(splices[0]), splice_text_1)
            self.assertEqual(str(splices[1]), splice_text_2)

    def test_cigar(self):
        """Test extraction of gaps from cigar"""

        cigars = [
            ["12M", [12], []],
            ["12M100N12M", [12, 12], [100]],
            ["12M100N12M200N14M", [12, 12, 14], [100, 200]],
            ["12M1I12M", [25], []],
            ["12M1D12M", [24], []]
        ]

        for cigar in cigars:
            custom_sam = self.base_sam
            custom_sam[5] = cigar[0]
            aln = SAM_Alignment.from_SAM_line("\t".join(custom_sam))
            seqs, gaps = Splice.collect_gaps(aln)
            self.assertEqual(seqs, cigar[1])
            self.assertEqual(gaps, cigar[2])

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
    test_splices = [
        Splice('SEQ_REGION', 1, 10, '-'),
        Splice('SEQ_REGION', 1, 10, '-', 5, 5),
        Splice('SEQ_REGION', 10, 20, '+')
    ]

    def test_splice_collection_creation(self):
        """Test splice collection creation"""

        col = SpliceCollection()
        self.assertEqual(col.size, 0)
        self.assertEqual(col.get_splices(), [])

    def test_splice_collection_add_splice(self):
        """Test splice collection add splice"""

        col = SpliceCollection()

        col.add_splice(self.test_splices[0])
        self.assertEqual(col.size, 1)
        col.add_splice(self.test_splices[1])
        self.assertEqual(col.size, 1)
        col.add_splice(self.test_splices[2])
        self.assertEqual(col.size, 2)

    def test_splice_collection_add_splices(self):
        """Test splice collection add a list of splices"""

        col = SpliceCollection()
        col.add_splices(self.test_splices)
        self.assertEqual(col.size, 2)

    def test_splice_collection_get_same_splice(self):
        """Test splice collection add a list of splices"""

        col = SpliceCollection()
        col.add_splices(self.test_splices)

        same1 = col.get_same_splice(self.test_splices[1])
        self.assertEqual(same1.key, "SEQ_REGION-1-10")
        same2 = col.get_same_splice(Splice("A", 1, 10))
        self.assertEqual(same2, None)

    def test_splice_collection_is_known(self):
        """Test splice collection is known"""

        col = SpliceCollection()
        col.add_splices(self.test_splices)

        known1 = col.is_known(self.test_splices[0])
        self.assertEqual(known1, True)
        known2 = col.is_known(Splice("A", 1, 10))
        self.assertEqual(known2, False)

    def test_splice_collection_left_is_known(self):
        """Test splice collection left is known"""

        col = SpliceCollection()
        col.add_splices(self.test_splices)

        known1 = col.left_is_known(self.test_splices[0])
        self.assertEqual(known1, True)
        known2 = col.left_is_known(Splice("SEQ_REGION", 1, 100))
        self.assertEqual(known2, True)
        known3 = col.left_is_known(Splice("SEQ_REGION", 2, 10))
        self.assertEqual(known3, False)

    def test_splice_collection_right_is_known(self):
        """Test splice collection right is known"""

        col = SpliceCollection()
        col.add_splices(self.test_splices)

        known1 = col.right_is_known(self.test_splices[0])
        self.assertEqual(known1, True)
        known2 = col.right_is_known(Splice("SEQ_REGION", 2, 10))
        self.assertEqual(known2, True)
        known3 = col.right_is_known(Splice("SEQ_REGION", 1, 100))
        self.assertEqual(known3, False)

if __name__ == '__main__':
        unittest.main()
