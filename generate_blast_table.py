#!/usr/bin/env python3

# *-----------------------------------------------------------------------------------------------
# | PROGRAM NAME: generate_blast_table.py
# | DATE: 26/08/22 
# | CREATED BY: Lila Maciel Rodriguez Perez
# | PROJECT FILE: https://github.com/limrp/generate_blast_table
# *-----------------------------------------------------------------------------------------------

# *-------------------------------------  Libraries ----------------------------------------------
import re
import argparse
import logging
import time
from Bio import SearchIO
from pprint import pprint
from pprint import pformat
import pandas as pd
from alive_progress import alive_bar

# pip install tqdm
# pip install openpyxl
# pip install progress
# pip install alive-progress

def parse_args():
    parser = argparse.ArgumentParser(description='Parse BLAST XML results and generate summary table.')
    parser.add_argument('-i', '--input', required=True, help='Input BLAST XML file.')
    parser.add_argument('-o', '--output_excel', required=True, help='Output Excel file.')
    parser.add_argument('-o2', '--output_csv', help='Output CSV file (optional).')
    parser.add_argument('-p', '--protein_families', help='File containing custom list of protein families.')
    parser.add_argument('-e', '--evalue', type=float, default=1e-10, help='E-value threshold.')
    parser.add_argument('-a', '--aln_span', type=int, default=50, help='Alignment span threshold.')
    parser.add_argument('-l', '--log', default='log.txt', help='Log file name (optional). Default is "log.txt".')
    return parser.parse_args()

def read_protein_families(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]

def find_families(hit_description, protein_families):
    families_found = []
    for family in protein_families:
        # Use \b to enforce that we only get whole word matches
        pattern = r'\b' + re.escape(family.lower()) + r'\b'
        # Checking if a pattern (family) was found
        if re.search(pattern, hit_description.lower()):
            families_found.append(family)
    return families_found

def pass_filter(hsp, e_value_threshold, aln_span_threshold):
    return hsp.evalue <= e_value_threshold and hsp.aln_span >= aln_span_threshold

def main():
    args = parse_args()
    # Get the log file name from parsed arguments
    log_file_name = args.log
    
    # Initialize logging (To log only to a file)
    logging.basicConfig(level=logging.INFO, filename=log_file_name, filemode='w')

    # To ALSO print the log messages to the console
    # logging.basicConfig(level=logging.INFO, filename=log_file_name, filemode='w')
    # console = logging.StreamHandler()  # Create console handler
    # console.setLevel(logging.INFO)  # Set level of console handler
    # logging.getLogger('').addHandler(console)  # Add console handler to the root logger

    # Constants
    E_VALUE_THRESHOLD = args.evalue
    ALN_SPAN_THRESHOLD = args.aln_span

    if args.protein_families:
        protein_families = read_protein_families(args.protein_families)
    else:
        # Variables: List of protein families
        protein_families = [
            "Adhesin",
            "Cadherin",
            "Catenin",
            "Collagen",
            "Fibronectin",
            "Alpha catulin",
            "Adhesion",
            "Claudin",
            "GAIN domain",
            "Immunoglobulin",
            "Integrin",
            "Intercellular adhesion molecule",
            "Invasin",
            "Laminin",
            "Lectin",
            "Plekstrin Homology Domain",
            "Selectin",
            "Triple Helix Repeat",
            "Vimentin",
            "Protocadherin",
            "Dystrophin",
            "Cell-cell adhesion",
            "Bacteria adhesin",
            "Aggrecan",
            "Adhesion GPCr",
        ]
    
    # Initialize dictionaries counters
    total_hits_counter_dict = {family: 0 for family in protein_families}
    filtered_hits_counter_dict = {family: 0 for family in protein_families}
    # Initialize dictionary to store unique hit ids using sets
    unique_hits_dict = {family: set() for family in protein_families}

    # Reading and Parsing BLAST XML
    print()
    with alive_bar(total = 1, title = 'P1: Reading XML file', unknown='classic') as bar:
        qresult = SearchIO.parse(args.input, "blast-xml")
        # Update the progress bar after finishing the process
        bar()
    # Logging
    logging.info("Finished reading the XML file")

    # Logging the start of the iteration
    print()
    logging.info("Starting the iteration through queries and hits")

    # Iterate through queries and hits
    print()
    with alive_bar(title = 'P2: Processing queries and hits', unknown='squares', theme = 'smooth') as bar:
        for query in qresult:
            for hit in query.hits:
                # Logging hit
                logging.info(f"$ Processing {hit.id} {hit.description}")
                # logging.info(f"$ Processing ID: {hit.id}, DES: {hit.description} ...")

                # Check if hit belongs to any family
                families = find_families(hit.description, protein_families)
                logging.info(f"List of families: {families}")

                for family in families:
                    if family in total_hits_counter_dict:
                        total_hits_counter_dict[family] += 1
                        logging.info(f"{family} match: {hit.id} is from the {family} family.")
                        logging.info(f"{family} TOTAL COUNT: {total_hits_counter_dict[family]}")
                    else:
                        logging.info(f"{family} is not in total_hits_counter_dict")
                    
                    # Check if hit passes the filter
                    if any(pass_filter(hsp, E_VALUE_THRESHOLD, ALN_SPAN_THRESHOLD) for hsp in hit.hsps):
                        filtered_hits_counter_dict[family] += 1
                        unique_hits_dict[family].add(hit.id)
                    else:
                        pass
            # Update the progress bar (after finishing with a query)
            bar()


    # Generate DataFrame
    print()
    with alive_bar(total = 1, title='P3: Generating DataFrame') as bar:
        data_dict = {
            "Protein Family": protein_families,
            "Total BLAST Hits": [total_hits_counter_dict[family] for family in protein_families], # Check if this column remains the same!
            "Filtered BLAST Hits": [filtered_hits_counter_dict[family] for family in protein_families],
            "Unique BLAST Hits": [len(unique_hits_dict[family]) for family in protein_families],
        }
        # Convert dictionary to DataFrame
        df = pd.DataFrame(data_dict)
        # Update the progress bar after finishing the process
        bar()
    # Logging after finishing the process
    # logging.info("DataFrame generated.")
    
    # Display DataFrame
    print()
    print()
    print(df.to_string(index=False))
    # To also store in log file
    logging.info("")
    logging.info(df.to_string(index=False))
    
    # Save to Excel
    print()
    with alive_bar(total = 1, title='P4: Saving to Excel') as bar:
        df.to_excel(args.output_excel, index=False)
        # Update the progress bar after finishing the process
        bar()
    # Logging after finishing the process
    logging.info(f"Excel file saved to {args.output_excel}\n")
    print()
    
    # Optionally save to CSV
    if args.output_csv:
        with alive_bar(total = 1, title='P5: Saving to CSV') as bar:
            df.to_csv(args.output_csv, index=False)
            # Update the progress bar after finishing the process
            bar()
        # Logging after finishing the process
        logging.info(f"CSV file saved to {args.output_csv}")
        print()

    for family in protein_families:
        logging.info(f"{family} final total count: {total_hits_counter_dict[family]}")

    # logging.info(f"Total hit counts by family")
    # logging.info(f"Adhesion TOTAL COUNT: {adhesion_total_hits_counter}")
    # logging.info(pprint(total_hits_counter_dict)) # printed to the console, not to the file
    # logging.info((total_hits_counter_dict)) # didn't print to the console but printed ugly version to the file
    dict_pretty_str = pformat(total_hits_counter_dict)
    logging.info("")
    logging.info(f"Dictionary with TOTAL counts of hits by family:")
    logging.info(dict_pretty_str)
    
if __name__ == '__main__':
    main()
