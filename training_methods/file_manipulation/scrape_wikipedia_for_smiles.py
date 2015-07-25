__author__ = 'ddeconti'


import random
import re
import sys
import time
import wikipedia

from bs4 import BeautifulSoup

'''
Scraping Wikipedia for SMILES for drugbank csv from Aqeel.
Deprecated...
'''

def scrape_for_smiles(name):
    '''
    Scrapes wikipedia for smiles string of given drug name
    :param name: name of drug
    :return: smiles string
    '''
    max_wait = 7
    time.sleep(max_wait*random.random())
    search_results = wikipedia.search(name)
    wiki_page = wikipedia.page(search_results[0])
    wiki_html = wiki_page.html()
    soup = BeautifulSoup(wiki_html, "html.parser")
    smiles_tag = soup.find("div",
                           {"style":
                                "word-wrap:break-word; text-indent:-1.5em"})
    if not smiles_tag:
        return None
    smiles_tag = str(smiles_tag)
    first = re.search("-1.5em", smiles_tag).start()
    last = re.search("/div", smiles_tag).start()
    smiles = smiles_tag[first+8:last-1]
    return smiles


def get_smiles(drug_dict):
    '''
    Puts smiles value from wikipedia to "smiles" key from drug dict
    :param drug_dict: dictionary of drugs with keys for smiles and counts
    :return: updated drug_dict with smiles key and value
    '''
    for drug_name in drug_dict:
        smiles = scrape_for_smiles(drug_name)
        if not smiles:
            continue
        drug_dict[drug_name]["smiles"] = smiles
    return drug_dict


def print_dict(drug_dict):
    '''
    Prints drug dict_dict keys and values to tab delimited stdout
    :param drug_dict:
    :return: None
    '''
    out_str = ["drug_name", "se_count", "gene_count", "smiles"]
    sys.stdout.write('\t'.join(out_str) + "\n")
    for drug_name in drug_dict:
        se_count = drug_dict[drug_name]["se_count"]
        gene_count = drug_dict[drug_name]["gene_count"]
        if "smiles" not in drug_dict[drug_name]:
            continue
        smiles = drug_dict[drug_name]["smiles"]
        out_list = [drug_name, se_count, gene_count, smiles]
        out_list = map(lambda x: str(x), out_list)
        out_str = '\t'.join(out_list)
        sys.stdout.write(out_str + '\n')


def parse_csv(filename):
    '''
    Parse Aqeel csv summary of database with drug names, # of side effects
    and gene expression interactions
    :param filename: name of csv file
    :return: drug as k,v with name and k,v counts
    '''
    drug_dict = {}
    try:
        handle = open(filename, 'rU')
    except IOError as e:
        sys.stderr.write("IOError: " + str(e) +
                         "\nError in parse_csv()\n")
        sys.exit()
    handle.readline()
    for line in handle:
        line = line.strip('\n').split(',')
        drug_name = line[1]
        try:
            gene_count = int(line[2])
            se_count = int(line[3])
        except ValueError as e:
            sys.stderr.write("ValueError: " + str(e) +
                             "\nError in parse_csv()\n")
            sys.exit()
        drug_dict[drug_name] = {"se_count": se_count, "gene_count": gene_count}
    handle.close()
    return drug_dict


def main(sa):
    csv_filename = sa[0]
    drug_dict = parse_csv(csv_filename)
    drug_dict = get_smiles(drug_dict)
    print_dict(drug_dict)

if __name__ == "__main__":
    main(sys.argv[1:])
