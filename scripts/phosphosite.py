import csv
import sys

#Script functions

def data_to_sites(proteins_positions_aa, mode):
    
    """
    Converts a list of three elements: protein, position, and amino acid into a phosphosite element.
    
    Parameters
    ----------
    arg1 : list
        list of length 3 containing protein, position, and amino acid.
        E.g., Q9Y4H2, 779, T
    arg2 : str
        str indicating the format of the phosphosite element. 
        Four possible modes exist:
        psp: PhosphoSitePlus database format e.g., Q9Y4H2   T779-p
        iptmnet: iPTMnet database format e.g., Q9Y4H2   T   779
        networkin: iPTMnet database format e.g., Q9Y4H2   779   T
        classical: format used to report phophosites in the literature and elsewhere, 
        also used by PhosphoSitePlus e.g., Q9Y4H2   T779

    Returns
    -------
    list
        returns a list of phosphosites elements as in MaxQuant each column containing the input 
        for this function can have more than one item separated by ;

    """

    #separate possible proteins sharing the same sequence and thus grouped under the same phosphosite
    list_split = list()
    for item in proteins_positions_aa:
        list_split.append(item.split(';'))

    #site formatting
    count = 0
    sites_list = list()

    if mode.lower() == 'psp':
        for i, item in enumerate(list_split[0]):
            mode_site = list_split[2][0] + list_split[1][i] + '-p'
            site = [list_split[0][i], mode_site]
            sites_list.append(site)
    elif mode.lower() == 'iptmnet':
        for i, item in enumerate(list_split[0]):
            site = [list_split[0][i], list_split[2][0], list_split[1][i]]
            sites_list.append(site)
    elif mode.lower() == 'networkin':
        for i, item in enumerate(list_split[0]):
            site = [list_split[0][i], list_split[1][i], list_split[2][0]]
            sites_list.append(site)
    else:
        #gives the flexibility of defaulting any other ptm mode to classic
        #PTM site list defaulting to classic mode, e.g., Q9Y4H2 T779
        for i, item in enumerate(list_split[0]):
            mode_site = list_split[2][0] + list_split[1][i]
            site = [list_split[0][i], mode_site]
            sites_list.append(site)

    return(sites_list)

def kinase_substrate_psp(site, kinase_subs_writer, path):
    
    """
    Finds a given phosphosites in classical format in the Kinase_Substrate_Dataset.txt dataset downloaded from PhosphoSitePlus
    
    Important for finding known kinase-substrate relationships.

    Parameters
    ----------
    arg1 : list
        list of length 2 containing the protein identifier and the modified residue.
        E.g., Q9Y4H2, T779
    arg2 : writer
        file handle for saving the kinase-substrate relationships found.

    Returns
    -------
    list
        saves the kinase-substrate relationships found for a given phosphosite to file

    """

    #open file for each iteration as with closes the file and we need to compare each site to the whole list
    with open(path + '/data/Kinase_Substrate_Dataset.txt', 'r') as phospho_tsvfile:
        kinase_subs_reader = csv.DictReader(phospho_tsvfile, dialect = 'excel-tab')
        kinase_subs_header = kinase_subs_reader.fieldnames
        kinase_reader = csv.DictReader(phospho_tsvfile, dialect = 'excel-tab')
        for kinase_subs_row in kinase_subs_reader:
            #comparisson for finding known kinase-substrate relationships fro the given phosphosite
            if site[0] == kinase_subs_row.get('SUB_ACC_ID') and site[1] == kinase_subs_row.get('SUB_MOD_RSD'):
                kinase_subs_data = list()
                for col in kinase_subs_header:
                    kinase_subs_data.append(kinase_subs_row.get(col))
                kinase_subs_writer.writerow(kinase_subs_data)           

#Script body

def main():
    
    #project_path = '/home/tesa/MyPhD/bioinformatics_pathology/final_project'
    project_path = sys.argv[1]

    #defines mode of phosphosite elements
    ptm_mode = sys.argv[2]

    print(ptm_mode)

    with open(project_path + '/data/Phospho(STY)Sites_filtered.txt', 'r') as tsvfile, \
    open(project_path + '/results/sites_{0}.tsv' .format(ptm_mode), 'w') as site_tsvfile, \
    open(project_path + '/data/Kinase_Substrate_Dataset.txt', 'r') as phospho_tsvfile, \
    open(project_path + '/results/kinase_substrate.tsv', 'w') as kinase_tsvfile:

        #upload MaxQuant output for phosphosite data (Phospho(STY)sites table)
        reader = csv.DictReader(tsvfile, dialect = 'excel-tab')

        #upload PhosphoSitePlus kinase-substrate relationship data
        kinase_reader = csv.DictReader(phospho_tsvfile, dialect = 'excel-tab')
        kinase_header = kinase_reader.fieldnames

        #phosphosite writer
        site_writer = csv.writer(site_tsvfile, delimiter='\t')

        #kinase_substrate writer
        kinase_writer = csv.writer(kinase_tsvfile, delimiter='\t')
        #save kinase-substrate relationship header
        kinase_writer.writerow(kinase_header)

        #iterate over Phospho(STY)sites file
        for row in reader:
            
            #the format of the 'Proteins' column in MaxQuant is one or a set of UniProt identifiers separated by ';'
            proteins = row.get('Proteins', None)

            #the format of the 'Positions within proteins' column in MaxQuant is a number or a set of numbers separated by ';'
            #depending on whether the 'Proteins' column contains one or more identifiers 
            positions = row.get('Positions within proteins', None)
            
            #amino acid (S, T, or Y) in the modified site, same for every identifier and positions in 
            #the 'Proteins' 'Positions within proteins' columns
            aa = row.get('Amino acid',None)

            #checks if there is information for phosphosite in the current row of the file
            phosphosite_info = [proteins, positions, aa]
            if not any(phosphosite_info):
                break

            sites = data_to_sites(phosphosite_info, ptm_mode)
            for site in sites:
                #save site information for later use with other tools within the script
                site_writer.writerow(site)
            
            #looks for kinase-substrate relationships in data downloaded from PhosphoSitePlus using the classical formatting
            #writes the relationships found to a file
            kinase_sites = data_to_sites(phosphosite_info, 'classic')
            for kinase_site in kinase_sites:
                kinase_substrate_psp(kinase_site, kinase_writer, project_path)
                
if __name__ == '__main__':
    main()