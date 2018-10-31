import os
import eHive


class SpeciesFilter(eHive.BaseRunnable):
    """Only keep species with bam files"""

    def run(self):
        species = self.param_required('species')
        bam_dir = self.param_required('bam_dir')

        if species in os.listdir(bam_dir):
            self.dataflow({}, 2)
