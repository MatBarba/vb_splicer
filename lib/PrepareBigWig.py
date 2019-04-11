#!env python3

import argparse
import logging

import eHive
import os
import json

import requests, sys

class PrepareBigWig(eHive.BaseRunnable):
    def param_default(self):
        return {
                rest_server: "https://www.vectorbase.org/rest",
        }

    def run(self):
        logging.basicConfig(level=logging.INFO)

        species = self.param('species')
        rest_server = self.param('rest_server')
        tmp_dir = self.param('tmp_dir')
        bigwig_dir = self.param('bigwig_dir')
        summary_dir = self.param('summary_dir')
        ftp_summary_dir = self.param('ftp_summary_dir')
        json_dir = self.param('json_dir')
        force_bigwig = self.param('force_bigwig')

        size_file = os.path.join(tmp_dir, species + ".sizes")

        # Get the size file for the bed -> bigwig conversion
        sizes, version = self.get_sizes(rest_server, species)

        # If only one bigwig file, only copy
        sp_bigwig_dir = os.path.join(bigwig_dir, species)
        bigwigs = [name for name in os.listdir(sp_bigwig_dir) if ".bw" in name]

        # Check if the final bigwig already exists
        do_bigwig = 1
        final_bw = os.path.join(summary_dir, species, species + ".bw")
        if os.path.isfile(final_bw) and not force_bigwig:
            do_bigwig = 0
        
        one_bigwig = 0
        if do_bigwig:
            if len(bigwigs) == 1:
                one_bigwig = 1
            else:
                # Print to file if it doesn't exist already
                if not os.path.isfile(size_file):
                    self.print_sizes(sizes, size_file)
        
        # Flow
        if do_bigwig:
            if one_bigwig:
                self.dataflow({
                    'version': version,
                    'file': final_bw,
                }, 2)
            else:
                self.dataflow({
                    'sizes': size_file,
                    'version': version,
                    'file': final_bw,
                }, 3)
        
        # Whatever happens, create the json file for this file
        ftp_file = os.path.join(ftp_summary_dir, species, species + ".bw")
        json_file = os.path.join(json_dir, species + ".bw.json")
        species_name = species.replace("_", " ")
        species_name = species_name[0].upper() + species_name[1:]
        description = "Merged bigwig track for %s" % species_name
        label = version + "_bw"
        
        self.dataflow({
            'file': final_bw,
            'version': version,
            'ftp_file': ftp_file,
            'description': description,
            'label': label,
            'track_file': ftp_file,
            'json_file': json_file
        }, 4)

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
