import os
import eHive


class SpliceDBFactory(eHive.BaseRunnable):
    """Create jobs to generate a list of bam + splice db files"""

    def run(self):
        # Get the list of gffs created
        species = self.param_required('species')
        bam_dir = self.param_required('bam_dir')
        splice_dir = self.param_required('splice_dir')

        # Prepare the splice dirs
        splice_species_dir = os.path.join(splice_dir, species)
        if not os.path.exists(splice_species_dir):
            os.makedirs(splice_species_dir)

        # Get the list of all bam files
        bam_species_dir = os.path.join(bam_dir, species)

        if os.path.exists(bam_species_dir):
            bams = filter(
                lambda x: x.endswith('.bam'), os.listdir(bam_species_dir)
            )

            for bam_file in bams:
                bam_path = os.path.join(bam_species_dir, bam_file)

                # Prepare the splice db file
                splice_filename = os.path.basename(bam_file)
                splice_filename = splice_filename.replace('.bam', '.sqlite')
                splice_db = os.path.join(splice_species_dir, splice_filename)

                self.dataflow({
                    'bam_file': bam_path,
                    'splice_db': splice_db
                }, 2)
