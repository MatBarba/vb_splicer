import unittest
from SimpleGTF import SimpleGTF


class TestSimpleGTF(unittest.TestCase):

    def test_parse_GTF(self):
        lines = (
            'KB671544	VectorBase	gene	43848	44327	.	+	.	gene_id "AEPI014344";',
            'KB671544	VectorBase	transcript	43848	44327	.	+	.	gene_id "AEPI014344"; transcript_id "AEPI014344-RA";',
            'KB671544	VectorBase	exon	43848	44327	.	+	.	gene_id "AEPI014344"; transcript_id "AEPI014344-RA"; exon_number "1";',
            'KB671544	VectorBase	CDS	43848	44324	.	+	0	gene_id "AEPI014344"; transcript_id "AEPI014344-RA"; exon_number "1";',
            'KB671544	VectorBase	start_codon	43848	43850	.	+	0	gene_id "AEPI014344"; transcript_id "AEPI014344-RA"; exon_number "1";',
            'KB671544	VectorBase	stop_codon	44325	44327	.	+	0	gene_id "AEPI014344"; transcript_id "AEPI014344-RA"; exon_number "1";'
        )

        genes = SimpleGTF.parse_GTF(lines)
        self.assertEqual(len(genes), 1)

        gene = genes[0]
        self.assertEqual(len(gene.transcripts), 1)

        transcript = gene.transcripts[0]
        self.assertEqual(len(transcript.exons), 1)


if __name__ == '__main__':
        unittest.main()
