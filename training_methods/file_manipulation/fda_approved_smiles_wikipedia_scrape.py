__author__ = 'ddeconti'
# Goes through fda approved drug list from Joe
# Scrapes wikipedia for the drug smiles
# Gets unique smiles string from set
# prints to screen for analysis outside script

import random
import re
import sys
import time
import wikipedia

from bs4 import BeautifulSoup

def scrape_for_smiles(name, max_wait=2):
    '''
    Scrapes wikipedia given drug name for smiles string
    :param name: name of drug as string
    :param max_wait: optional default wait to not overload search requests
    :return: smiles string
    '''
    time.sleep(max_wait*random.random())
    try:
        search_results = wikipedia.search(name)
    except:
        return None
    if len(search_results) < 1:
        return None
    try:
        wiki_page = wikipedia.page(search_results[0])
    except:
        return None
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


def parse_csv(filename):
    '''
    Parses given comma-delimited csv file to get nda type drug names to scrape
    for smiles
    :param filename: name of csv
    :return: set of smiles strings
    '''
    try:
        handle = open(filename, 'rU')
    except IOError as e:
        sys.stderr.write("IOError: " + str(e) +
                         "\nError in parse_csv()\n")
        sys.exit()
    smiles_set = set([])
    handle.readline()
    for line in handle:
        arow = line.strip('\n').split(',')
        app_type = arow[-2]
        if app_type != "nda":
            continue
        drug_name = arow[0]
        alt_name = arow[1]
        if re.search("nda", drug_name):
            end = re.search("nda", drug_name).start()-2
            drug_name = drug_name[:end]
        smiles = scrape_for_smiles(drug_name)
        if not smiles:
            smiles = scrape_for_smiles(alt_name)
        if smiles != None and smiles not in smiles_set:
            smiles_set.add(smiles)
            sys.stdout.write(drug_name + "\t" + smiles + '\n')
            sys.stdout.flush()
    handle.close()


def print_smiles(smiles_set):
    '''
    Prints set of smiles string to stdout
    :param smiles_set: set of unique fda approved smiles strings
    :return:
    '''
    for smiles in smiles_set:
        sys.stdout.write(smiles + '\n')
        sys.stdout.flush()


def main(sa):
    '''
    Parses arg

    :param sa:
    :return:
    '''
    csv_filename = sa[0]
    parse_csv(csv_filename)
    #print_smiles(smiles_set)


if __name__ == "__main__":
    main(sys.argv[1:])
