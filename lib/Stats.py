#!env python3

import argparse
import logging

from Splice import SpliceDB

import eHive
from glob import glob
import os.path

class Stats(eHive.BaseRunnable):

    queries = [
            ['genes', "SELECT count(*) FROM genes"],
            ['genes_without_intron', "SELECT count(*) FROM genes WHERE introns=0"],
            ['genes_completely_covered', "SELECT count(*) FROM genes WHERE introns=covered_introns AND introns > 0"],
    #            ['genes_completely_covered_10', "SELECT count(*) FROM genes WHERE introns=covered_introns AND introns > 0 AND splice_coverage >= 10"],
            ['genes_partially_covered', "SELECT count(*) FROM genes WHERE covered_introns > 0 AND covered_introns < introns AND introns > 0"],
            ['genes_not_covered', "SELECT count(*) FROM genes WHERE covered_introns = 0 AND introns > 0"],
            ['introns', "SELECT sum(introns) FROM genes"],
            ['known_introns', "SELECT SUM(covered_introns) FROM genes"],
            ['known_introns_%genes', "SELECT SUM(covered_introns)*1.0/sum(introns) FROM genes"],
            ['known_introns_genes_coverage', "SELECT sum(splice_coverage)*1.0/sum(covered_introns) FROM genes"],
    #            ['splices', "SELECT count(*) FROM splices"],
    #            ['genes_with_weak_splices', "SELECT count(*) FROM genes g LEFT JOIN splices s USING(gene) WHERE s.coverage < g.splice_coverage/100 AND s.coverage < 10"],
            ]
        
    # Do the same using transcripts instead of genes
    duplicated_queries = []
    for q in queries:
        if 'gene' in q[0] or 'gene' in q[1]:
            newq = [q[0].replace('gene', 'transcript'), q[1].replace('gene', 'transcript')]
            duplicated_queries.append(newq)
    queries += duplicated_queries

    def run(self):
        logging.basicConfig(level=logging.INFO)

        splice_dbs = self.param_required('splice_dbs')
        stats_path = self.param_required('stats_path')

        # Run it
        Stats.make_stats(splice_db, stats_path)

    def make_stats(inputs, output):
        stats = []
        paths = sorted(glob(inputs))
        for input in paths:
            stats.append(Stats.get_stats(input))
        Stats.print_stats(stats, output)

    def get_stats(input):
        input_db = SpliceDB(input)
        
        # Make a few queries...
        path = os.path.basename(input).replace('.sqlite', '')
        stats = {
                'db': path
                }
        for q in Stats.queries:
            name = q[0]
            query = q[1]
            #query = query.replace('gene', 'transcript')
            stats[name] = input_db.query_count(query)

        return stats
    
    def print_stats(stats, output):
        
        names_order = [q[0] for q in Stats.queries]
        header = "#" + "\t".join(names_order)
        
        stat_lines = [header]
        for stat in stats:
            line = []
            for name in names_order:
                line.append(str(stat[name]))
            stat_line = "\t".join(line)
            stat_lines.append(stat_line)

        if output != '-':
            with open(output, 'w') as out:
                for stat_line in stat_lines:
                    out.write(stat_line + "\n")
        else:
            for stat_line in stat_lines:
                print(stat_line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Splice db to use')
    parser.add_argument('output', help='Output to this tab file')
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

    Stats.make_stats(args.input, args.output)


if __name__ == "__main__":
    main()
