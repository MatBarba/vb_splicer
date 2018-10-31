import os
import eHive


class BamFactory(eHive.BaseRunnable):
    """Create jobs to generate a list of bam files"""

    def run(self):
        # Get the list of gffs created
        species = self.param_required('species')
        bam_dir = self.param_required('bam_dir')

        # Get the list of all bam files
        species_bam_dir = os.path.join(bam_dir, species)

        if os.path.exists(species_bam_dir):
            bams = os.listdir(species_bam_dir)

            for bam_file in bams:
                bam_path = os.path.join(species_bam_dir, bam_file)
                self.dataflow({
                    'bam_file': bam_path
                }, 2)
