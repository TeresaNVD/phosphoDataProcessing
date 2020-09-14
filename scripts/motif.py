import csv
import sys
import re
import os

#Script functions

def hprd_to_regex(motif_kinase):
    
    """
    Converts a given HPRD motif (natural language) to regular expression.

    Important for finding potential kinase-substrate relationships in the vicinity of phosphosites identified in the experiment.
    
    Parameters
    ----------
    arg1 : list
        list of length 2 containing motif and kinase.
        E.g., RXRXX[pS/pT][F/L]	Akt kinase substrate motif.

    Returns
    -------
    list
        returns a list of length 2 containing the motif as regex and the kinase.

    """
    
    motif = motif_kinase[0]
    kinase = motif_kinase[1]

    #replaces separators and indicators of phosphorylation / *
    string_to_replace = ['*', '/']
    for string in string_to_replace:
        motif = motif.replace(string, '')
    
    #replaces any amino acid indicator X with .
    motif = motif.replace('X', '.')
    
    #replaces the indicator of phosphorylation pS/pT/pY for S/T/Y
    motif = motif.replace('pS', 'S')
    motif = motif.replace('pT', 'T')
    motif = motif.replace('pY', 'Y')

    #replaces 'substrate motif' to shorten the kinase information
    kinase = kinase.replace(' substrate motif', '')

    kinase_motif = [kinase, motif]
    return(kinase_motif)

def kinase_motif_scan(motif_regex, sequences_list):
    
    """
    Scans a list of sequences for a given motif.

    The scan is not centered and is meant to find all the possible kinase-substrate 
    relationships in the vicinity of a given phosphosite.
    
    Parameters
    ----------
    arg1 : list
        list of length 2 containing motif as regex and kinase.
    arg2 : list
        list of sequences centered and of a given length.

    Returns
    -------
    list
        returns a list of length 2 of kinase-substrate relationships in the vicinity of a given phosphosite.

    """
    
    motif_regex = re.compile(motif_regex)
    
    motif_match = list()
    sequences_match = list()

    for item in sequences_list:
        #performs scan
        if motif_regex.findall(item):
            motif_match.append(motif_regex.findall(item))
            sequences_match.append(item)
    
    if not motif_match:
        motif_sequences_match = [motif_regex.pattern, 'NO match']
    else:
        motif_sequences_match = [motif_match, sequences_match]
    return(motif_sequences_match)

def sequence_logo(project_path, input_sequence_path, output_logo_path):
    
    #graphical representation of an amino acid or nucleic acid multiple sequence alignment.
    #http://weblogo.threeplusone.com/manual.html#API

    #stdin stdout
    os.system('weblogo --size large --color-scheme chemistry --format SVG < ' + project_path + input_sequence_path + ' > ' + project_path + output_logo_path)

#Script body

def main():
    
    #project_path = '/home/tesa/MyPhD/bioinformatics_pathology/final_project'
    project_path = sys.argv[1]

    #file containing the sequences to scan
    seq_file = sys.argv[2]

    with open(project_path + '/results/{0}.tsv' .format(seq_file), 'r') as tsvfile:

        #uploads sequences found in the experiment, centered and at a given sequence window
        reader = csv.reader(tsvfile, delimiter='\t')

        sequences = list()
        for row in reader:
            #row is a list of length 1
            sequence = ''.join(row)
            #list of sequences to scan later on
            sequences.append(sequence)

    with open(project_path + '/data/HPRD_motifs_Jul2018.tsv', 'r') as tsvfile:

        #read the HPRD kinase/phosphatase motif list one by one
        #http://hprd.org/serine_motifs
        reader = csv.reader(tsvfile, delimiter='\t')
        
        with open(project_path + '/results/HPRD_motifs_Jul2018_as_regex.tsv', 'w') as out_tsvfile, \
        open(project_path + '/results/kinase_motif_scan_{0}.tsv' .format(seq_file), 'w') as out_scan_tsvfile:

            #creates writer to save the HPRD motif converted to regex
            motif_regex_writer = csv.writer(out_tsvfile, delimiter='\t')
            header = ['Kinase', 'HPRD motif as regex']
            motif_regex_writer.writerow(header)

            #creates writer to save the motif scan results
            motif_scan_writer = csv.writer(out_scan_tsvfile, delimiter='\t')
            header = ['Kinase', 'HPRD motif', 'HPRD motif as regex', 'Motif hit', 'Sequence']
            motif_scan_writer.writerow(header)

            for row in reader:
                #processes motif list in natural language to regex in python
                kinase_motif_regex = hprd_to_regex(row)
                #saves motif conversion to file
                motif_regex_writer.writerow(kinase_motif_regex)

                #scans the list of sequences with each HPRD motif as regex one at a time
                kinase = kinase_motif_regex[0]
                motif = kinase_motif_regex[1]

                motif_sequence_hits = kinase_motif_scan(motif, sequences)
                if motif_sequence_hits[1] == 'NO match':
                    scan_result = [kinase, row[0], kinase_motif_regex[1], 
                                    motif_sequence_hits[0], motif_sequence_hits[1]]
                    motif_scan_writer.writerow(scan_result)
                else:
                    for i, value in enumerate(motif_sequence_hits[0]):
                        scan_result = [kinase, row[0], kinase_motif_regex[1], 
                                        ';'.join(motif_sequence_hits[0][i]), motif_sequence_hits[1][i]]
                        motif_scan_writer.writerow(scan_result)
    
    #generates sequence logo for all the sequences identified
    #reads sequences centered and of a given length generated in exp_summary.py
    input_sequence_path = '/results/{0}.tsv' .format(seq_file)
    output_logo_path = '/results/{0}.svg' .format(seq_file)
    sequence_logo(project_path, input_sequence_path, output_logo_path)

    #generates sequence logo for all the sequences identified that match the CK2 consensus motif
    #uses grep command for generating a file containing the sequences that match the CK2 consensus motif
    input_kinase_motif_scan = '/results/kinase_motif_scan_{0}.tsv' .format(seq_file)
    input_sequence_path = '/results/ck2_hit_{0}.tsv' .format(seq_file)
    output_logo_path = '/results/ck2_hit_{0}.svg' .format(seq_file)

    os.system('grep \'Casein Kinase II\' ' + project_path + input_kinase_motif_scan + ' | cut -f5 | grep -v "NO match" > ' + project_path + input_sequence_path)
    sequence_logo(project_path, input_sequence_path, output_logo_path)

if __name__ == '__main__':
    main()    