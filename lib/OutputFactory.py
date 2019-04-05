#!env python3

import argparse
import logging
import eHive
import os
import requests, sys
import json

from Splice import SpliceDB

class OutputFactory(eHive.BaseRunnable):
    def run(self):
        logging.basicConfig(level=logging.INFO)

        species = self.param_required('species')
        gff_dir = self.param_required('gff_dir')
        json_dir = self.param_required('json_dir')
        splice_db = self.param_required('splice_db')
        ftp_gff_dir = self.param_required('ftp_gff_dir')

        min_splices = 10000
        max_splices = 500000

        if not os.path.exists(gff_dir):
            os.makedirs(gff_dir)
        if not os.path.exists(json_dir):
            os.makedirs(json_dir)

        
        outputs = [
                { "category": 'known', "coverages": [1] },
#                { "category": 'unknown', "coverage": 1 },
#                { "category": 'unknown', "coverage": 2 },
#                { "category": 'unknown', "coverage": 5 },
#                { "category": 'unknown', "coverage": 10 },
#                { "category": 'unknown', "coverage": 100 },
#                { "category": 'unknown', "coverage": 1000 },
#                { "category": 'unknown', "coverage": 10000 },
#                { "category": 'unknown', "coverage": 100000 },
                { "category": 'unknown', "coverages": [100000, 10000, 1000, 100, 10, 5, 2, 1] },
                { "category": 'duplicates', "coverages": [1] },
        ]

        rest_server = self.param_required('rest_server')
        version = self.get_version(rest_server, species)
        db = SpliceDB(splice_db)
        
        for output in outputs:
            for cov in output['coverages']:

                # Check that there is enough data to create this
                n = 0
                if output['category'] == 'unknown':
                    n = db.count_tag(antitag='unknown', coverage=cov)
                else:
                    n = db.count_tag(tag=output['category'], coverage=cov)
                logging.info("%s splices to extract for %s : %s (%s)" % (str(n), species, output['category'], cov))

                # No splices, no track
                if n == 0:
                    logging.info("No splices to extract, skip this file")
                    continue

                # Not enough splices to make an informative track (filtered unknown only)
                if output['category'] == 'unknown' and cov > 1 and n < min_splices:
                    logging.info("Not enough splices to extract, skip this file")
                    continue

                if output['category'] == 'duplicates' and n == 0:
                    logging.info("Not enough splices to extract, skip this file")
                    continue

                # Too much splices to extract
                if output['category'] == 'unknown' and n > max_splices:
                    logging.info("Too much splices to extract, end this category")
                    break

                file_parts = [species, version, output['category']]
                if cov > 1:
                    file_parts.append(str(cov))
                file_name = "_".join(file_parts)
                track_name = file_name + ".gff"
                json_name = file_name + ".json"
                track_file = os.path.join(gff_dir, track_name)
                ftp_track_file = os.path.join(ftp_gff_dir, track_name)
                json_file = os.path.join(json_dir, json_name)

                label = self.make_label(version, output['category'], cov)
                description = self.make_description(species, output['category'], cov)

                # GFF
                self.dataflow({
                    'species': species,
                    'category': output["category"],
                    'coverage': cov,
                    'track_file': track_file
                }, 2)

                # Json
                self.dataflow({
                    'label': label,
                    'description': description,
                    'track_file': ftp_track_file,
                    'json_file': json_file,
                    'version': version,
                }, 3)
    
    def make_label(self, version, cat, cov):
        parts = [version, cat]
        if cov > 1:
            parts.append(str(cov))
        return "_".join(parts)

    def make_description(self, species, cat, cov):
        coverage = ""
        if cov > 1:
            coverage = ' (above %sx coverage)' % cov
        str_species = species[0].upper() + species[1:]
        str_species = str_species.replace('_', ' ')
        desc = "Summary track for %s, showing %s splice-junctions" % (str_species, cat)
        desc += coverage
        return desc

    def get_version(self, rest_server, species):
        ext = "/info/species"
        r = requests.get(rest_server + ext, headers={ "Content-Type" : "application/json"})
        
        # Too bad if we fail, but not catastrophic
        if not r.ok:
            r.raise_for_status()
        #    sys.exit()

        data = json.loads(r.text)
        for sp_data in data['species']:
            if sp_data['name'] == species:
                return sp_data['assembly']

