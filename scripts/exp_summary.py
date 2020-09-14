import csv
import sys

#Script functions

def file_header(dict_reader, project_path):
    """
    Returns and saves to file the header of the input file.
    
    Parameters
    ----------
    arg1 : dict
        reader with mapped information into a dict
    arg2 : str
        local path to the project folder

    Returns
    -------
    list
        list of dictionary keys
    """

    #get column names
    header = dict_reader.fieldnames
    with open(project_path + '/results/header_info.tsv', 'w') as tsvfile:
        for col_number, col_name in enumerate(header):
            #save column names to file
            col_number_name = (col_number + 1, col_name)
            output = "{0[0]}\t{0[1]}\n".format(col_number_name)
            tsvfile.write(output)
    #return column names
    return(header)

def summary_col(dict, col_tag, mode):

    """
    Returns a list of unique elements found in one or more columns.

    It provides a quick summary of the data in those columns.

    Parameters
    ----------
    arg1 : dict
        a given row of the file being analyzed
    arg2 : str
        part of the name or name of a column(s) you want to summarize.
        E.g.,use 'raw' for summarizing data from the columns: Best localization raw file, 
        Best score raw file, Best PEP raw file.
    arg3: str
        match mode: 'partial' for summarizing the information on one or more columns containing tag,
        or 'complete' for summarizing the information on a specific column.
    
    Returns
    -------
    list
        list of unique values in the selected column(s).

    """

    mode_list = ['partial', 'complete']
    if mode in mode_list:
        if mode == 'complete':
            value = dict.get(col_tag, None)
            value = value.split(";")
            col_info = value
        elif mode == 'partial':
            col_info = list()
            for key, value in dict.items():
                #case insensitive
                if col_tag in key.lower():
                    col_info.append(value)
        #return a list of unique elements
        if not col_info:
            print('Column tag not found')
            col_info = None
    else:
        print('Mode incorrect; must be \'partial\' or \'complete\'')
        col_info = None
    
    return(col_info)
    
def flatten_list(list_of_lists):
    
    """
    flatten a list of lists and returns a list of unique elements.

    """

    list_unique = list()
    for x in list_of_lists:
        for y in x:
            list_unique.append(y)
    list_unique = list(set(list_unique))
    return(list_unique)

def sequence_window(sequences_list, arm):
    
    """
    Subsets a list of sequences of the same length centered at the modified residue.

    Parameters
    ----------
    arg1 : list
        list of sequences of same length centered at the modified residue.
    arg2 : int
        how many characters at either end of the modified residue to include.
        
    Returns
    -------
    list
        list of subsetted sequences.
        E.g., for xxxxxxxSxxxxxxx input of length 7, if arm = 6 it gives xxxxxxSxxxxxx of length 6.
        arm = 0, returns modified residue;
        arm = 6, returns PhosphoSitePlus database format (13 amino acids);
        arm = 7, returns MaxQuant short sequence format (15 amino acids).

    """

    #check that input sequences are all same length
    sequences_len = list(set(map(len, sequences_list)))
    original_len = int(sequences_len[0])
    if len(sequences_len) == 1 and original_len%2 != 0:
        #check that the arm parameter will generate a subsetted sequence smaller than input
        sequences_recentered_list = list()
        if (arm * 2) + 1 < original_len:
            #subset the sequences
            start = ((original_len - 1)/2)-arm
            #account for right-open intervals
            end = (((original_len-1)/2)+arm)+1
            for sequence in sequences_list:
                sequence_recentered = sequence[start:end]
                sequences_recentered_list.append(sequence_recentered)
        else:
            print('Nothing to do, length for output sequences is greater than : ' + str(original_len))
        return(sequences_recentered_list)
    else:
        print('Input sequences are not centered, length is : ' + str(original_len))

#Script body

def main():

    #project_path = '/home/tesa/MyPhD/bioinformatics_pathology/final_project'
    project_path = sys.argv[1]

    #assumes that the input file is located under the data subfolder in the final_project folder 
    with open(project_path + '/data/Phospho(STY)Sites.txt', 'r') as tsvfile:

        #reads MaxQuant output for phosphosites: Phospho(STY)sites.txt
        reader = csv.DictReader(tsvfile, dialect = 'excel-tab')
        
        #saves header to a new file and returns header as column names
        col_names = file_header(reader, project_path)

        #check if all columns to summarize exist in file this columns are always found in the input file
        summary_col_names = ['Potential contaminant', 'Gene names', 
                            'Proteins', 'Modification window', 'Sequence window']
        if set(summary_col_names).issubset(col_names) == False:
            print('Important columns to summarize not found check validity of input file')
            sys.exit()
    
    #open experiment file again as the previous with statement closed the file handler
    #create a new data set filtered from contaminants and reverse hits (input for the other scripts)
    with open(project_path + '/data/Phospho(STY)Sites.txt', 'r') as tsvfile, \
    open(project_path + '/data/Phospho(STY)Sites_filtered.txt', 'w') as filter_tsvfile:

        reader = csv.DictReader(tsvfile, dialect = 'excel-tab')

        #writes header/column names to file
        filter_writer = csv.writer(filter_tsvfile, delimiter='\t')
        filter_writer.writerow(col_names)

        #summarizes data on the following columns
        raw_files_list = list()
        contaminants_list = list()
        reverse_list = list()
        genes_list = list()
        proteins_list = list()
        sequences_list = list()
        modifications_list = list()
        
        #iterate over data
        for row in reader:
            
            #raw files used as input for MaxQuant and thus used in generating the data
            raw_files_list.append(summary_col(row, 'raw', 'partial'))
            #a proper and complete input file should always have at least one column containing raw file information
            if not any(raw_files_list):
                sys.exit()

            #contamination is expected to occur thus we filter the contaminants first
            if row.get('Potential contaminant') == "+":
                #contaminants identified in the experiment
                contaminants = summary_col(row, 'Protein names', 'complete')
                if any(contaminants):
                    contaminants_list.append(contaminants)

            #false positives are expected so FDR is estimated by target-decoy(reverse)
            #reverse hits are filtered
            elif row.get('Reverse') == "+":
                #reverse hits identified in the experiment
                reverse = summary_col(row, 'Protein', 'complete')
                if any(reverse):
                    reverse_list.append(reverse)

            else:
                #writes the contaminants and reverse filtered dataset
                filter_line = list()
                for col in col_names:
                    filter_line.append(row.get(col))
                filter_writer.writerow(filter_line)

                #gene names identified in the experiment
                gene_names = summary_col(row, 'Gene names', 'complete')
                if any(gene_names):
                    genes_list.append(gene_names)
                
                #proteins identified in the experiment
                proteins = summary_col(row, 'Proteins', 'complete')
                if any(proteins):
                    proteins_list.append(proteins)

                #modification types identified
                modifications = summary_col(row, 'Modification window', 'complete')
                if any(modifications):
                    modifications_list.append(modifications)

                #sequences identified in the experiment
                sequences = summary_col(row, 'Sequence window', 'complete')
                if any(sequences):
                    sequences_list.append(sequences)

   
    #saves to file the summary for each column
    #all lists contain unique elements
    data_summary = [raw_files_list, contaminants_list, reverse_list, 
                    genes_list, proteins_list, modifications_list, 
                    sequences_list]
    data_name = ['raw_files', 'contaminants', 'reverse_hits', 
                 'gene_names', 'proteins', 'modifications', 
                 'sequences']
    for i, data_type in enumerate(data_summary):
        #generates files names with unique names provided by the data_name list
        f = open(project_path + '/results/{0}.tsv' .format(data_name[i]), 'w')
        data_type_flat = flatten_list(data_type)
        for item in data_type_flat:
            f.write(str(item)+'\n')
        f.close()

    #generates list of sequences subsetted in two frequently used lengths (13 and 15 aminoacids)
    sequences_list_flat = flatten_list(sequences_list)
    common_centered_len = [6, 7]
    common_centered_dict = dict()
    for len in common_centered_len:
        common_centered_dict[len] = sequence_window(sequences_list_flat, len)
    for key in common_centered_dict.keys():
        #writes subsetted sequences to file with unique file names
        file = open(project_path + '/results/sequence{0}.tsv' .format(key), 'w')
        for item in common_centered_dict[key]:
            file.write(str(item)+'\n')
        file.close()

if __name__ == '__main__':
    main()