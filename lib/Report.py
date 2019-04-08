import eHive


class Report(eHive.BaseRunnable):
    """Write a report for a given species"""

    def run(self):
        gffs = {}

        # Get the list of gffs created
        if self.param_exists('files'):
            files = self.param_required('files')
        species = self.param_required('species')

        report = self.create_report(species, files)
        self.param('report_text', report)

    def create_report(self, species, files):
        n_files = 0
        text_lines = []

        if species in files:
            sp_files = files[species]
            n_files = len(sp_files)
            for file in sp_files:
                text_lines.append("\t%s" % file)

        text_lines.insert(0, '%d GFF files for %s' % (n_files, species))
        return"\n".join(text_lines)

    def write_output(self):
        self.dataflow({
            'species': self.param('species'),
            'report': self.param('report_text'),
        }, 2)
