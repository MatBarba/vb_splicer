import eHive


class Report(eHive.BaseRunnable):
    """Write a report for a given species"""

    def run(self):
        gffs = {}

        # Get the list of gffs created
        if self.param_exists('gffs'):
            gffs = self.param_required('gffs')
        species = self.param_required('species')

        report = self.create_report(species, gffs)
        self.param('report_text', report)

    def create_report(self, species, gffs):
        n_files = 0
        text_lines = []

        if species in gffs:
            gff_files = gffs[species]
            n_files = len(gff_files)
            for gff in gff_files:
                text_lines.append("\t+ %s" % gff)

        text_lines.insert(0, '%d GFF files for %s' % (n_files, species))
        return"\n".join(text_lines)

    def write_output(self):
        self.dataflow({
            'species': self.param('species'),
            'report': self.param('report_text'),
        }, 2)
