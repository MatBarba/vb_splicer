#!env python3

import argparse
import logging

import eHive
from glob import glob
import os.path
import json

class Apollo_track(eHive.BaseRunnable):

    def run(self):
        logging.basicConfig(level=logging.INFO)

        infile = self.param_required('track_file')
        outfile = self.param_required('json_file')
        label = self.param_required('label')
        version = self.param_required('version')
        description = self.param_required('description')

        # Run it
        self.make_json(infile, outfile, label, description, version)

    def make_json(self, infile, outfile, label, description, version):
        print("From %s create %s" % (infile, outfile))

        caption = label
        category = "Summary track"
        source_url = infile
        source_type = 'gff'
        if ".bw" in source_url:
            source_type = 'bigwig'

        data = {
                'label': label,
                'metadata' : {
                    'caption': caption,
                    'category': category,
                    'description': description,
                    'display': 'off',
                    'source_type': source_type,
                    'source_url': source_url,
                    'version': version
                    }
                }

        with open(outfile, 'w') as out:
            out.write(json.dumps(data))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='File to create as a track for Apollo')
    parser.add_argument('output', help='Output json file')
    parser.add_argument('label', help='Label of the track')
    parser.add_argument('description', help='Description of the track')
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

    Apollo_track.make_json(args.input, args.output, args.label, args.description)


if __name__ == "__main__":
    main()
