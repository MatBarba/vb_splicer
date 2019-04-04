#!env python3

import argparse
import logging

import eHive
from glob import glob
import os.path

class Apollo_tracks(eHive.BaseRunnable):

    def run(self):
        logging.basicConfig(level=logging.INFO)

        infile = self.param_required('in_file')
        outfile = self.param_required('out_json')

        # Run it
        Stats.make_json(infile, outfile)

    def make_json(infile, outfile):
        print("From %s create %s" % (infile, outfile))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='File to create as a track for Apollo')
    parser.add_argument('output', help='Output json file')
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

    Apollo_tracks.make_json(args.input, args.output)


if __name__ == "__main__":
    main()
