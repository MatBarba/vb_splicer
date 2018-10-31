import os
import eHive


class GtfFactory(eHive.BaseRunnable):
    """Create jobs to dump a GTF for each species"""

    def run(self):
        # Get the list of gffs created
        species = self.param_required('species')
        gtf_dir = self.param_required('gtf_dir')
        out_file_stem = self.param_required('out_file_stem')
        force_gtf = self.param_required('force_gtf')

        gtf_file = os.path.join(gtf_dir, species + '.' + out_file_stem)

        # Do we need to recreate a gtf file?
        if not force_gtf and os.path.exists(gtf_file):
            self.warning("GTF file already exits: skip dumper")
        else:
            # Send a dumper job
            self.dataflow({}, 2)

        # Continue
        self.dataflow({'gtf_file': gtf_file}, 3)
