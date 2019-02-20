#!env python3

import argparse
import logging
import eHive

class GFFFactory(eHive.BaseRunnable):
    def run(self):
        species = self.param_required('species')
        
        outputs = [
                { "category": 'known', "coverage": 1 },
                { "category": 'unknown', "coverage": 1 },
                { "category": 'unknown', "coverage": 10 },
                { "category": 'unknown', "coverage": 100 },
                { "category": 'unknown', "coverage": 1000 },
                { "category": 'duplicates', "coverage": 1 },
        ]

        for output in outputs:
            self.dataflow({
                'species': species,
                'category': output["category"],
                'coverage': output["coverage"],
            }, 2)

