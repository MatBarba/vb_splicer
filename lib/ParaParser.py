#!env python3

import argparse
import logging
import os

import eHive
import xml.etree.ElementTree as ET

class ParaParser(eHive.BaseRunnable):
    def run(self):
        do_paralogs = self.param('do_paralogs')
        if not do_paralogs:
            return

        logging.basicConfig(level=logging.INFO)
        xml_dir = self.param('homologs_dir')
        output_dir = self.param('paralogs_dir')

        ParaParser.parse(xml_dir, output_dir)
    
    def parse(xml_dir, output_dir):
        # Get the list of xml files
        xml_files = ParaParser.get_xmls(xml_dir)
        cur_path = ''
        paralogs = {}
        for xml_path, xml_name in xml_files:
            if cur_path != xml_path:
                cur_path = xml_path
                logging.info("Parsing xmls in " + cur_path)
            xml_file = os.path.join(xml_path, xml_name)
            logging.debug("Parsing file " + xml_file)

            tree = ET.parse(xml_file)
            root = tree.getroot()

            genes = {}
            gene_species = {}
        
            # Get species
            for child in root:
                if "species" in child.tag:
                    species = child.attrib["name"]
                    for db in child:
                        for gene in db[0]:
                            gid = gene.attrib["id"]
                            gene_name = gene.attrib["geneId"]
                            genes[gid] = gene_name
                            gene_species[gid] = species
                elif "groups" in child.tag:
                    for group in child:
                        if "paralogGroup" in group.tag:
                            para_genes = []
                            for prop in group:
                                if "geneRef" in prop.tag:
                                    gid = prop.attrib["id"]
                                    perc_id = prop[0].attrib["value"]
                                    para_genes.append({"gid": gid, "perc_id": float(perc_id)})

                            # Are all paralogs above the identity threshold?
                            min_id = 50
                            if len(para_genes) == 2:
                                if para_genes[0]["perc_id"] > min_id and para_genes[0]["perc_id"] > min_id:
                                    gid1 = para_genes[0]["gid"]
                                    gid2 = para_genes[1]["gid"]
                                    gene1 = genes[gid1]
                                    gene2 = genes[gid2]
                                    species = gene_species[gid1]

                                    if not species in paralogs:
                                        paralogs[species] = {}

                                    if not gene1 in paralogs[species]:
                                        paralogs[species][gene1] = []

                                    if not gene2 in paralogs[species]:
                                        paralogs[species][gene2] = []

                                    paralogs[species][gene1].append(gene2)
                                    paralogs[species][gene2].append(gene1)

        # Write the result
        ParaParser.print_paralogs(output_dir, paralogs)

    def print_paralogs(output_dir, paralogs):
        os.makedirs(output_dir, exist_ok=True)

        for species in paralogs:
            sp_file = output_dir + "/" + species + ".txt"

            logging.info("Write %d paralogs in %s" % (len(paralogs[species]), sp_file))

            with open(sp_file, 'w') as out:
                for gene1 in paralogs[species]:
                    for gene2 in paralogs[species][gene1]:
                        out.write("%s\t%s\n" % (gene1, gene2))

    def get_xmls(path):
        xmls = []
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename.endswith(".xml"):
                    xmls.append((dirpath, filename))
        return xmls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to homologs xml root dir")
    parser.add_argument("output", help="Path to paralogs dir where files will be created")
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

    ParaParser.parse(args.input, args.output)


if __name__ == "__main__":
    main()
