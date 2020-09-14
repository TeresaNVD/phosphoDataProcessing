import Bio
from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import ExPASy, SwissProt
import urllib,urllib2
import time
import csv
from datetime import date, datetime

#Script functions

def fetch_uniprot_fasta(uniprot_acc):
    
    """
    Fetches the corresponding fasta sequence of a given UniProt accession 
    or identifiers (SwissProt and Trembl).
    
    Parameters
    ----------
    arg1 : str
        UniProt accession or identifier from both SwissProt and Trembl.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        holds sequence and information about it.

    """

    if not uniprot_acc:
        print('Please input a UniProt accession or identifier')
    else:
        url = "http://www.uniprot.org/uniprot/" + uniprot_acc + ".fasta"
        while True: 
            try:
                request = urllib2.Request(url)
                #email address here to help us debug in case of problems.
                contact = "tnunezde@uwo.ca"
                request.add_header('User-Agent', 'Python %s' % contact)
                response = urllib2.urlopen(request)
                #fetch fasta sequence for a given UniProt identifiers (isoform sensitive)
                fasta_iterator = SeqIO.parse(response, "fasta")
                for seq in fasta_iterator:  
                    #sequence = str(seq.seq)
                    return(seq)
            except urllib2.HTTPError, e:
                if e.code == 404:
                    print('Please input a valid UniProt accession or identifier')
                    sys.exit()
                #catch 503 Service Unavailable. The server is currently unable to handle the request 
                #due to a temporary overloading or maintenance of the server.
                if e.code == 503:
                    wait_time = int(e.hdrs.get('Retry-after', 0))
                    print('Sleeping %i seconds...' % (wait_time,))
                    #waits for a while and coninues the request
                    time.sleep(wait_time)
                    continue 
            break

def dating(date_to_convert):
    
    """
    Converts a date to ISO format.
    
    Parameters
    ----------
    arg1 : str
        date to conver.

    Returns
    -------
    str
        date in ISO format: yyyymmdd.
    """
    #Converts either collection or record deposition dates to ISO format

    #returns None if date to convert is empty
    if date_to_convert is None:
        return None

    #Converts date to ISO by matching a list of possible date formats to the date
    #Found date types: ['2015'], ['01-Jun-2015'], ['Feb-2017'], ['2016-11-22'], ['2016-08'], None
    for fmt in ('%Y','%d-%b-%Y','%b-%Y', '%Y-%m-%d', '%Y-%m'):
        try:
            dt = datetime.strptime(date_to_convert, fmt)
            return(dt.date().isoformat())
        except ValueError:
            pass

def fetch_swp_expasy(uniprot_acc): 
    
    """
    Fetch information on SwissProt accession (manually reviewed UniProt entry).
    
    http://biopython.org/DIST/docs/api/Bio.SwissProt.Record-class.html

    Parameters
    ----------
    arg1 : str
        SwissProt accession or identifier.

    Returns
    -------
    list
        list of length 2 with the name of the attributes found and their values.
    """

    #generates record object with information regarding SwissProt identifier
    handle = ExPASy.get_sprot_raw(uniprot_acc)
    record = SwissProt.read(handle)

    #checks all the attributes possibles for the record object generated and their type
    #attributes are of type: str, tuple, or list
    #attribute list found here: http://biopython.org/DIST/docs/api/Bio.SwissProt.Record-class.html
    attrib_names = ['accessions', 'data created', 'date created (ISO)','organism','gene names', 
                    'description', 'comments', 'keywords']
    swp_info_list = [record.accessions, record.created[0], dating(record.created[0]), record.organism, 
                    record.gene_name, record.description, record.comments, record.keywords]
    return(attrib_names, swp_info_list)

#Script body

project_path = sys.argv[1]
id = sys.argv[2]

def main():
    
    with open(project_path + '/results/sequences.fasta', 'a') as sequence_fasta, \
    open(project_path + '/results/swp_info.tsv', 'a') as swp_tsvfile:

        #fetch UniProt sequence in fasta format
        fseq = fetch_uniprot_fasta(id)
        sequence_fasta.write(fseq.format("fasta"))

        #writes to file a summary of the protein's function
        #which includes identifiers, gene name, functional information
        info_names, swp_info = fetch_swp_expasy(id)
        for i, name in enumerate(info_names):
            if type(swp_info[i]) is list:
                if name == 'comments':
                    row = '\n'.join(swp_info[i]) + '\n'
                else:
                    row = info_names[i] + ': ' + ';'.join(swp_info[i]) + '\n'
            else:
                row = info_names[i] + ': ' + swp_info[i] + '\n'
            swp_tsvfile.write(row)

if __name__ == '__main__':
    main()
