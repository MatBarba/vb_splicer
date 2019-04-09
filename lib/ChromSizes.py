#!env python3

import argparse
import logging

import eHive
import os
import json

import requests, sys

class ChromSizes(eHive.BaseRunnable):
    def param_default(self):
        return {
                rest_server: "https://www.vectorbase.org/rest",
        }

    def run(self):
        logging.basicConfig(level=logging.INFO)

        species = self.param('species')
        rest_server = self.param('rest_server')
        tmp_dir = self.param('tmp_dir')

        size_file = os.path.join(tmp_dir, species + ".sizes")
        sizes, version = self.get_sizes(rest_server, species)
        self.print_sizes(sizes, size_file)

        self.dataflow({
            'sizes': size_file,
            'version': version
        }, 2)

    def get_sizes(self, rest_server, species):
        ext = "/info/assembly/%s?" % (species)
        logging.info(rest_server + ext)
        r = requests.get(rest_server + ext, headers={ "Content-Type" : "application/json"})
        
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        
        json_data = json.loads(r.text)
        sizes = []
        for seq in json_data['top_level_region']:
            sizes.append([seq["name"], str(seq["length"])])
        logging.info("Got %d seq_regions", len(sizes))

        version = json_data['assembly_name']
        
        return(sizes, version)

    def print_sizes(self, sizes, size_file):
        logging.info("Print to %s", size_file)
        with open(size_file, 'w') as size_output:
            for seq_size in sizes:
                size_output.write("\t".join(seq_size) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Splice db to use')
    parser.add_argument(
        '--gtf', dest='gtf', help='GFT file with genes features')
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

    Tagger.tag_splices(args.input, args.gtf)


if __name__ == "__main__":
    main()
