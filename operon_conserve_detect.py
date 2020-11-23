'''
@author: ichaudr

Determines how conserved a specific reference operon is within a specific taxonomic clade. 

'''

from Bio import pairwise2, Seq, Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Blast import NCBIWWW, NCBIXML
from features import AnnotatedHit
from genome_fragment import GenomeFragment
from species import Species
from tqdm import tqdm
import datetime
import time
import json
import traceback
import csv
import random
import xml.dom.minidom
import re
import sys
import csv
import threading
import uuid
import os
import shutil


#****INPUT JSON FILE PATH****#
INPUT_JSON = "input.json" 
#****************************#

##################################################
# Parameters are loaded from the input json file #
##################################################

#Entrez request parameters
REQUEST_LIMIT = 5
SLEEP_TIME = .5
EMAIL = ""
E_API = ""

#Blast parameters
tax_include = []
tax_exclude = []
database = 'ref_prok_rep_genomes'
e_val = 10-10
coverage_min = .8
max_hits = 500
max_blast_attemts = 3
blast_search_mult_factor = 2
annotate = True
extensive_search = True
reverse_blast = True

#Hit feature detection parameters
margin_limit = 15
max_feature_detect_attempts = 10
feature_search_mult_factor = 4

#Operon assembly parameters
use_reference_threshold = True
ref_threshold_margin = .5
feature_limit = 3
intergenic_limit = 1500

#Other paramters
thread_limit = 10
species_percent_id_limit = 0.45
color_code = {}

#Input records
input_records = []
reference_genome_accession = ''
reference_assembly_accession = ''
reference_genome_name = ''

#Output parameters
cache_dir = './cache/'
output_dir = './output/{run_id}/'
run_id = str(datetime.date.today()) + '_' + str(uuid.uuid4())

##################################################
##################################################

#Data related to the reference operon
ref_genome_frag = None
ref_features = []

#Reverse BLAST db directories
reverse_blast_root = './reverse_blast/{ref_assembly_accession}'
reference_total_protein = './reverse_blast/{ref_assembly_accession}/{ref_assembly_accession}.fasta'
reference_blast_db = './reverse_blast/{ref_assembly_accession}/{ref_assembly_accession}_blastdb'


#The list of GenomeFragment objects that the BLAST hits get sorted into
genome_frags = []

#The list of species
species = []

def search_blast(input_records, db='ref_prok_rep_genomes', max_attempts=3, search_mult_factor=2, max_hits=50, e_cutoff=10E-10, tax_incl=[], tax_excl=[], annotate = True, min_cover=None, extensive_search=True):
    '''
    Performs blast search for a set of records. 

    Parameters
    ----------
    input_records : list[string]
        The list of accession numbers to conduct the BLAST search.
    db: string, optional
        The database to use for the BLAST search. The refseq_representative_genomes is set as the default.
    max_attempts: int, optional
        The BLAST search will be repeated this many times until all hits withtin the e-value cutoff are returned.
    search_mult_factor: int, optional
        The max_hits will be multiplied by this factor following every incomplete BLAST search. The default is 2. 
    max_hits : int, optional
        Tne starting number of max number of hits to return. The default is 50. This will change if the E-cutoff is not met by the last hit returned
    e_cutoff : float, optional
        The threshold for the E-value of the hits. The default is 10E-10.
    tax_incl : list[int], optional
        The taxa to include in the BLAST search. The default is None.
    tax_excl : list[int], optional
        The taxa to exclude in the BLAST search. The default is None.
    annotate:
        Specifies whether or not to return a AnnotatedHit object or the tuple formatted data, see below for return. 
    min_cover: float, optional
        The minimum coverage of the hits
    extensive_search: bool
        Specifies if the BLAST search should be exhaustive or if it should only iterate once. 


    Returns
    -------
    hits : list[(input_record, hit_record ,alignment_object)] (if annotate=False)
        A list containing the associated input record, the accession for the hit, and the alignment object for the hit.
    
    annotated_hits: list[AnnotateHit.object]
        A list of AnnotatedHit.objects that hold metadata for each of the BLAST hits. 

    '''

    #Holding the max attempts in a temp variable so that it can reset for each input record passed in
    temp_max_hits = max_hits
    
    #If extensive search is off, the max_attempts is set to 1
    if not extensive_search:
        max_attempts = 1


    print("|~> BLAST search: " + str(input_records) + "...")
    
    #check the legnth of the input_records. 
    if len(input_records) < 1:
        raise Exception("Need at least one protein record to conduct BLAST search.")
    
    #The final list of BLAST hits. May be AnnotatedHit objects or tuples depending on the annotate parameter. See above. 
    return_hits = []
    
    
    #Gets the accession numbers for all the hits in the BLAST search and
    # appends hits[] with every unique record. 
    for input_record in input_records:

        #Reset the max_hits for each input record
        max_hits = temp_max_hits
        
        #Fetches the protein record based off the accession number 
        print("\t|~> Getting protein record for " + str(input_record) + "...")
        
        
        for i in range(REQUEST_LIMIT):
            
            try:
            
                handle = Entrez.efetch("protein", id=input_record, rettype="fasta",
                                       retmode="text")
                
                time.sleep(SLEEP_TIME)
                break
                
            except:
                
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                time.sleep(SLEEP_TIME)
                
                if i == (REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
        
        
        #Fetches the protein sequence to be used in the BLAST search
        print("\t|~> Getting protein sequence for " + str(input_record) + "...")
        input_seq = (SeqIO.read(handle, "fasta")).format("fasta")
        

        #Keeps track of the current number of attempts made to complete the BLAST search
        current_attempts = 0

        #If false, the BLAST search has not returned the last hit within the e-value cutoff
        # and the max_attempts has not been reach.
        blast_complete = False

        #Continue repeating BLAST searches until one of the following is met:
        #   - max_attempts is reached
        #   - or all hits within the e-value cutoff are returned. 
        while not blast_complete:

            #Check if the max_attempts has been reached
            if current_attempts == max_attempts:
                print("\t\t|~> Max number of BLAST search attempts has been reached.")
                blast_complete = True
                continue
            
            print("\t|~> Performing BLAST search " + str(current_attempts + 1))

            
            #Performs the appropriate BLAST search based
            if len(tax_incl) > 0 or len(tax_excl) > 0:

                #Holds the parameter for the taxonomic limitation set
                taxon = ""


                #Appends the taxon parameter with all included taxa IDs
                if len(tax_incl) == 1:

                    taxon = "txid" + str(tax_incl[0]) + "[orgn]"

                elif len(tax_incl) > 1:

                    #Goes through each of the taxa limits appends to the overall entrez_query parameter
                    for i in range(len(tax_incl) - 1):
                        taxon = taxon + "txid" + str(tax_incl[i]) + "[orgn]" + " AND "

                    taxon = taxon + "txid" + str(tax_incl[-1]) + "[orgn]"

                
                #Appends the taxon parameter with all the excluded taxa IDs
                if len(tax_excl) == 1:

                    taxon = taxon + " NOT " + "txid" + str(tax_excl[0]) + "[orgn]"

                elif len(tax_excl) > 1:

                    taxon = taxon + " NOT "

                    #Goes through each of the taxa limits appends to the overall entrez_query parameter
                    for i in range(len(tax_excl) - 1):
                        taxon = taxon + "txid" + str(tax_excl[i]) + "[orgn]" + " NOT "

                    taxon = taxon + "txid" + str(tax_excl[-1]) + "[orgn]"
                
                

                #Send BLAST request with all the parameters
                for i in range(REQUEST_LIMIT):
                
                    try:
                    
                        result_handle = NCBIWWW.qblast("tblastn", db ,input_seq, 
                                                entrez_query=taxon, expect=e_cutoff,
                                                hitlist_size=max_hits)
                        
                        #Parses the resulting hits as a list
                        print("\t\t|~> Getting BLAST result records")
                        blast_records = list(NCBIXML.parse(result_handle))
                        
                        time.sleep(SLEEP_TIME)
                        break
                        
                    except:
                        
                        print("\t\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                        time.sleep(SLEEP_TIME)
                        
                        if i == (REQUEST_LIMIT - 1):
                            print("\t\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

            #Perform BLAST search if there is no taxon limit set
            else:            
                
                for i in range(REQUEST_LIMIT):
                
                    try:
                    
                        result_handle = NCBIWWW.qblast("tblastn", db ,input_seq, 
                                                expect=e_cutoff,
                                                hitlist_size=max_hits)
                        
                        #Parses the resulting hits as a list
                        print("\t\t|~> Getting BLAST result records")
                        blast_records = list(NCBIXML.parse(result_handle))
                        
                        time.sleep(SLEEP_TIME)
                        break
                        
                    except:
                        
                        print("\t\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                        time.sleep(SLEEP_TIME)
                        
                        if i == (REQUEST_LIMIT - 1):
                            print("\t\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
                

            #If the number of returned hits is zero, continue the loop
            if len(blast_records[0].alignments) == 0:
                print("\t|~> BLAST search returned no hits. Reattempting...")
                current_attempts = current_attempts + 1
                continue

            #If the number of returned hits is less than the max_hits requested, then all hits within the e-value threhold should have been returned and the BLAST search is complete
            if len(blast_records[0].alignments) < max_hits:
                blast_complete = True
                print("\t|~> BLAST search was successful")
                continue 
            
            #Check the E-value of the last hit returned. If it is within the threhold, the BLAST search is complete. Otherwise, the search repeated with a larger max_hit. 
            if blast_records[0].alignments[-1].hsps[0].expect > e_cutoff:
                blast_complete = True
                print("\t|~> BLAST search was successful")
                continue
            else:
                print("\t|~> BLAST search was not complete. Increasing max_hits...")
                max_hits = max_hits * search_mult_factor
                print("\t|~> Max hits: " + str(max_hits) + ". Reattempting...")
                current_attempts = current_attempts + 1
        


        #Adds each unique accession number to hits[]
        print("\t|~> Extracting hits from BLAST results...")
        for record in blast_records[0].alignments:
            
            current_hit_def = re.sub('[^A-Za-z0-9]+', '_', record.hit_def)
            curr_hit_rec = record.hit_id.split('|')[-2]
            print("\t\t|~> Analyzing hit " + str(curr_hit_rec))
            
            #Iterate through the hits
            for hit in record.hsps:

                #Initiates a AnnotatedHit object if set by the parameters. 
                if annotate:
                    a_hit = AnnotatedHit(query_accession=input_record, hit_accession=curr_hit_rec, genome_fragment_name=current_hit_def, align_start=hit.sbjct_start, alignment_seq=hit.sbjct, 
                        align_end=hit.sbjct_end, strand=hit.frame[1], percent_identity=(hit.identities/hit.align_length), req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME)
                        
                #Checks if hit meets the minimum coverage if provided 
                if min_cover:

                    #Calculate the coverage for the current hit                  
                    cov = (hit.query_end - hit.query_start + 1) / (len(input_seq))
                    
                    if(cov >= min_cover):

                        #Appends the AnnotatedHit object if requested
                        if annotate:
                            return_hits.append(a_hit)

                        else:
                            return_hits.append((input_record, curr_hit_rec ,record))
                            
                    #Prints error if the minimum coverage is not met    
                    else:
                        print("\t\t|~> Hit did not meet coverage requirement: " + str(curr_hit_rec))
                else:
                    
                    #Appends the AnnotatedHit object if requested
                    if annotate:

                        if(len(return_hits) == 0):
                            print("\t\t|~> Adding first hit: " + str(curr_hit_rec))
                            return_hits.append(a_hit)
                        elif(not(a_hit in return_hits)):
                            print("\t\t|~> Adding hit: " + str(curr_hit_rec))
                            return_hits.append(a_hit)

                    else:
                        #Check if the hit is already in the return list
                        if len(return_hits) == 0:
                            print("\t\t|~> Adding first hit: " + str(curr_hit_rec))
                            return_hits.append((input_record, curr_hit_rec ,record))
                        elif (not (curr_hit_rec in list(zip(*return_hits))[1])):
                            print("\t\t|~> Adding hit: " + str(curr_hit_rec))
                            return_hits.append((input_record, curr_hit_rec ,record))
                

    print("\t|~> Returning " + str(len(return_hits)) + " unique hits")
    return return_hits

def genome_fragment_exists(fragment_name, genome_accession):
    '''
    Checks if an GenomeFragment object with name fragment_name and gene_accession has already been constructed.

    Parameters
    ----------
    fragment_name: string
        The name of the GenomeFragment of interest.
    gene_accession: string
        The accession code the genome record of the GenomeFragment
    
    Returns
    -------
    does_exist: (bool, GenomeFragment)
        Tuple where the first element is True if exists, false otherwise. The second element is the matching GenomeFragment object, or None if it doesnt exist. 

    '''

    does_exist = (False, None)

    global genome_frags

    for fragment in genome_frags:
        if fragment.name == fragment_name and fragment.genome_accession == genome_accession:
            does_exist = (True, fragment)
    

    return does_exist

def species_exits(assembly_accession):
    '''
    Checks if species requested exists in the species list

    Parameters
    ----------
    species_name: string
        The name of the requested species
    
    Returns
    -------
    does_exit: (bool, Species)
        Tuple where the first element is True if exists, false otherwise. The second element is the matching Species object, or None if it doesnt exist. 
    '''

    does_exist = (False, None)

    global species

    for sp in species:
        if sp.assembly_accession == assembly_accession:
            does_exist = (True, sp)   

    return does_exist

def load_input_file(filename):
    '''
    Loads all the paramters from the input JSON.

    Parameters
    ----------
    filename: string
        The input file located inside the /input directory that is in the same directory as this .py script. 
    '''

    #Open the JSON file
    file_reader = json.load(open(filename))

    #Get all global variables and assign them
    
    #Entrez request parameters
    global REQUEST_LIMIT 
    REQUEST_LIMIT = file_reader['entrez'][0]['request_limit']

    global SLEEP_TIME
    SLEEP_TIME = file_reader['entrez'][0]['sleep_time']

    global EMAIL
    EMAIL = file_reader['entrez'][0]['email']

    global E_API
    E_API = file_reader['entrez'][0]['api_key']

    #Blast parameters
    global tax_include
    tax_include = file_reader['blast'][0]['tax_include']

    global tax_exclude
    tax_exclude = file_reader['blast'][0]['tax_exclude']

    global database
    database = file_reader['blast'][0]['database']

    global e_val
    e_val = file_reader['blast'][0]['e_val']

    global coverage_min
    coverage_min = file_reader['blast'][0]['coverage_min']

    global max_hits
    max_hits = file_reader['blast'][0]['max_hits']

    global max_blast_attempts
    max_blast_attempts = file_reader['blast'][0]['max_attempts']

    global blast_search_mult_factor
    blast_search_mult_factor = file_reader['blast'][0]['search_mult_factor']

    global annotate
    annotate = file_reader['blast'][0]['annotate']

    global extensive_search
    extensive_search = file_reader['blast'][0]['extensive_search']

    global reverse_blast
    reverse_blast = file_reader['blast'][0]['reverse_blast']

    #Hit feature detection parameters
    global margin_limit
    margin_limit = file_reader['hit_feature_detection'][0]['margin_limit']

    global max_feature_detect_attempts
    max_feature_detect_attempts = file_reader['hit_feature_detection'][0]['max_attempts']

    global feature_search_mult_factor
    feature_search_mult_factor = file_reader['hit_feature_detection'][0]['mult_factor']

    #Operon assembly parameters
    global feature_limit
    feature_limit = file_reader['operon_assembly'][0]['feature_limit']

    global intergenic_limit
    intergenic_limit = file_reader['operon_assembly'][0]['intergenic_limit']

    global use_reference_threshold
    use_reference_threshold = file_reader['operon_assembly'][0]['use_ref_limit']

    global ref_threshold_margin
    ref_threshold_margin = file_reader['operon_assembly'][0]['ref_limit_margin']

    #Other parameters
    global thread_limit
    thread_limit = file_reader['thread_limit']

    global species_percent_id_limit
    species_percent_id_limit = file_reader['species_percent_id_limit']

    global color_code
    color_code = file_reader['color_code'][0]

    #Input records
    global input_records
    input_records = file_reader['input_records']

    global reference_genome_accession
    reference_genome_accession = file_reader['reference_genome_accession']

    global reference_genome_name
    reference_genome_name = file_reader['reference_genome_name']

    global reference_assembly_accession
    reference_assembly_accession = file_reader['reference_genome_assembly']

    #Output paramaters
    global cache_dir
    cache_dir = file_reader['cache_dir']

def write_all_out(filename='output.csv'):
    '''
    Creates the output directory and outputs everything to that folder. The final output directory will have the following:
        - A CSV file holding the species name, structure similarity score, percent id for each reference gene, taxid, and assembly accession
        - An exact copy of the input JSON
        - A folder named SVG that will hold the diagrams for each species.  

    Parameters
    ----------
    filename: string
        The output file that will be written in the output/ directory
    
    Returns
    -------
    None
    '''

    #Set the full path to the output file
    filename = output_dir + filename

    #Make the output stream to the file
    ofile_stream = csv.writer(open(filename, mode='w'))

    #Setup and write the header row to the file 
    header_row = ['Species Name','Structural Similarity', 'Average Percent Amino Acid Identity']

    for record in input_records:
        header_row.append(str(record))

    header_row.append('Taxonomic ID')
    header_row.append('Genome Assembly Accession')
    ofile_stream.writerow(header_row)

    #Write info for all the species in the file
    for sp in species:
        sp_row = [sp.species_name, (str(sp.sim_score*100) + "%"), (str(100 * sum(sp.query_percent_ids.values())/len(sp.query_percent_ids)) + '%')]

        #Pulls the query percent IDs
        for query in input_records:
            sp_row.append(str(sp.query_percent_ids[query] * 100) + '%')
        
        #Pulls the taxid
        sp_row.append(sp.taxid)

        #Pulls the genome assembly accession
        sp_row.append(sp.assembly_accession)
        
        ''''
        #Pulls all the genome accessions associated with the species
        for accession in sp.genome_fragments_accessions:
            sp_row.append(accession)'''
        
        #Write row to file
        ofile_stream.writerow(sp_row)

def pre_process_fragments():
    '''
    Will fetch the assembly accessions, species name, and taxids for all the fragments.
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    '''
    taxids = {}
    species_names = {}
    assembly_accessions = {}

    accessions = []
    for frag in genome_frags:
        accessions.append(frag.genome_accession)


    curr_search_acc = accessions
    cont = True
    total_processed = 0
    total_to_process = len(accessions)
    missed_fragments = 0
    print('Total number of fragments: ' + str(total_to_process))

    while cont:
        print('Downloading records...')
        for i in range(REQUEST_LIMIT):

                try:

                    handle = Entrez.efetch(db="nuccore", id=curr_search_acc, seq_start=0, seq_stop=1, rettype='gb', retmode='XML')
                    records = list(Entrez.read(handle, 'xml'))
                    time.sleep( SLEEP_TIME)

                    break

                except:

                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

                    if i == (REQUEST_LIMIT - 1):
                            print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")


        for record in records:
            curr_acc = record['GBSeq_accession-version']

            assembly_found = False
            for info in record['GBSeq_xrefs']:
                if info['GBXref_dbname'] == 'Assembly':
                    total_processed = total_processed + 1
                    assembly_found = True
                    assembly_accessions[curr_acc] = info['GBXref_id']

            if not assembly_found:
                missed_fragments = missed_fragments + 1
                print('Assembly not found for: ' + curr_acc)

            for feature in record['GBSeq_feature-table']:
                if feature['GBFeature_key'] == 'source':
                    for qual in feature['GBFeature_quals']:
                        if qual['GBQualifier_name'] == 'db_xref':
                            taxids[curr_acc] = (qual['GBQualifier_value'].split(':')[1])
                        if qual['GBQualifier_name'] == 'organism':
                            species_names[curr_acc] = (qual['GBQualifier_value'])
            
            last_acc = curr_acc

            if last_acc == accessions[len(accessions)-1]:
                cont = False
            elif last_acc == accessions[len(accessions)-2]:
                curr_search_acc = accessions[len(accessions) - 1]
            else:
                curr_search_acc = accessions[accessions.index(last_acc): len(accessions) - 1]
            
            if (total_processed + missed_fragments) >= total_to_process:
                cont = False

    #print(assembly_accessions)

    #If no assembly accession is available, the fragment will be removed
    frags_to_remove = []
    for fragment in genome_frags:
        if fragment.genome_accession in assembly_accessions.keys():
            fragment.assembly_accession = assembly_accessions[fragment.genome_accession]
        else:
            frags_to_remove.append(fragment)

        if fragment.genome_accession in species_names.keys():
            fragment.species_name = species_names[fragment.genome_accession]

        if fragment.genome_accession in taxids.keys():
            fragment.taxid = taxids[fragment.genome_accession]
    
    for frag in frags_to_remove:
        genome_frags.remove(frag)

def process_frag(fragment, lock):
    '''
    Will complete the GenomeFragment object that is passed in by:
        1. Fetching all the features for the fragment
        2. Assign the features for the hit
        3. Assemble the operon
    
    Note: This function was implemented so that the fragments could be processed on multiple threads. 
    
    Parameters
    ----------
    fragment: GenomeFragment object
        The fragment to be processed.
    lock: ThreadLock
        Needed to print out to screen in sync
    '''
    fragment.fetch_record()
    fragment.fetch_features()
    fragment.fetch_hit_features(margin_limit=margin_limit, max_attempts=max_feature_detect_attempts, mult_factor=feature_search_mult_factor)
    fragment.purge_hits()
    fragment.assemble_operons(feature_limit=feature_limit, intergenic_limit=intergenic_limit)
    lock.acquire()
    tqdm.write("Completed:\n" + str(fragment) + "-"*50)
    lock.release()

    #Clear up memory by deleting the full features list
    del fragment.full_record
    fragment.full_record = None

def process_reference():
    '''
    Pulls information for the reference operon. Creates a GenomeFragment object and isolates the features that are apart of the operon. 
    The reference nucelotide accession and protein accessions should have already been parsed from the input JSON file. 

    Parameters
    -----------
    None

    Returns
    -------
    None
    '''
    
    global ref_genome_frag
    global ref_features

    #Make a GenomeFragment object for the reference nucelotide record
    ref_genome_frag = GenomeFragment(name=reference_genome_name, genome_fragment_accession=reference_genome_accession, req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME, cache_directory=cache_dir)
    ref_genome_frag.fetch_features()

    #Go through all of the features for this nucleotide accession and pull the GenomeFeature objects for the reference operon
    reference_features = []
    for record in input_records:
        for feat in ref_genome_frag.all_features:
            if feat.protein_accession == record:
                reference_features.append(feat)
    
    ref_features = reference_features

def get_reference_intergenic_distance():
    '''
    Determines the maximum intergenic distance in the reference operon. 

    Parameters
    ----------
    None - already parsed from the input JSON

    Returns
    -------
    intrgnc_dist: int
    The maximum distance between two genes in the reference operon.
    '''

    global ref_genome_frag
    global ref_features
    global ref_threshold_margin

    #Check if the data for the reference operon has been fetched. If not, fetch the information
    if ref_genome_frag == None and len(ref_features) == 0:
        process_reference()
    
    #Determine the maximum intergenic distance
    max_intergenic_distance = 0

    for i in range(len(ref_features)-1):
        #Calculate the intergenic distance between the current feature, i, and the next feature, i+1
        temp_ig_dist = ref_features[i].get_intergenic_distance(ref_features[i+1])
        if  temp_ig_dist > max_intergenic_distance:
            max_intergenic_distance = temp_ig_dist
    
    #Set the intergenic limit and considers the margin
    global intergenic_limit
    intergenic_limit = (1 + ref_threshold_margin) * max_intergenic_distance

def calculate_percent_ids(sp):
    '''
    Determines the percent identity of each of the hits in a Species object by aligning the entire feature sequence with the sequence of the reference hit. This replaces the BLAST percent identity.

    Parameters
    ----------
    None

    Returns
    -------
    None
    '''

    global ref_features

    #A list holding the symbols of the traditional amino acids
    trad_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    #Iterate through all of the features in every fragment in each species
    
    for frag in sp.genome_fragments:
        for feat in frag.hits:
            
            #Pull the sequence for this feature
            feat_seq = feat.aa_sequence

            #Hold the reference sequence
            ref_seq = 'None'

            #Pull the sequence for the reference gene this hit came from
            for ref_feat in ref_features:
                if feat.query_accession == ref_feat.protein_accession:
                    ref_seq = ref_feat.aa_sequence
            
            #Purge the reference sequence and the hit feature sequence of any non traditional amino acids
            purged_ref_seq = ''
            purged_feat__seq = ''

            if ref_seq == None:
                ref_seq = ''
            
            if feat_seq == None:
                feat_seq = ''

            for c in ref_seq:
                if c in trad_aa:
                    purged_ref_seq = purged_ref_seq + c

            for c in feat_seq:
                if c in trad_aa:
                    purged_feat__seq = purged_feat__seq + c
            
            ref_seq = purged_ref_seq
            feat_seq = purged_feat__seq

            
            #print('Performing alignment between: ' + str(feat_seq) + ' and ' + str(ref_seq))
            
            #Perform the pairwise alignment
            alignments = pairwise2.align.globaldx(ref_seq, feat_seq, blosum62)

            #Calculate the percent identity for this alignment
            percent_matches = []
    
            for align in alignments:
                
                #The format of an alignment:
                    # [0] : first sequence aligned
                    # [1] : second sequenced aligned
                    # [2] : alignment score
                    # [3] : starting position of the alignment
                    # [4] : ending position of the alignment 
                
                ref_seq_aligned = align[0]
                feat_seq_aligned = align[1]
                matches = 0
                
                #For a global alignment, the start position is always 0
                align_size = align[4]
                size_adj = 0
                
                for x in range(align_size):
                    
                    #size_adj is subtracted from the length of the alignment to remove
                    # gapped positions. 
                    if ref_seq_aligned == "-" or feat_seq_aligned == "-":
                        size_adj += 1
                        continue
                    
                    if ref_seq_aligned[x] == feat_seq_aligned[x]:
                        matches += 1
                        
                #The size of the alignment is adjusted for the gapped positions
                current_percent_matches = (matches/(align_size-size_adj))
                
                percent_matches.append(current_percent_matches)
                #print("The percent match: ", current_percent_matches)
            
            #The average percent match is calculated for the ensemble of alignments 
            # that was returned. 
            total_match_percentage = 0
            
            for match_percentage in percent_matches:
                total_match_percentage += match_percentage
            if len(percent_matches) != 0:
                average_percent_similar = total_match_percentage / len(percent_matches)
            else:
                average_percent_similar = 0
            
            #print(average_percent_similar)
            
            feat.percent_identity = average_percent_similar

def make_reference_blastdb():
    '''
    Makes the local BLAST database out of the reference to be used for the reverse blast filter. 

    Parameters
    ---------
    None

    Returns
    --------
    None
    '''

    global reference_assembly_accession
    global reverse_blast_root
    global reference_total_protein
    global reference_blast_db

    #Make the neccessary directories: the root directory in the 
    if not os.path.exists(reverse_blast_root.format(ref_assembly_accession=reference_assembly_accession)):
        os.mkdir(reverse_blast_root.format(ref_assembly_accession=reference_assembly_accession))

        #Gets all genomes associated with the reference genome assembly by searching the nucleotide DB with the genome assembly accession
        for i in range(REQUEST_LIMIT):

            try:
                handle = Entrez.esearch(db='nuccore', term=reference_assembly_accession, retmode='xml', retmax='5000')
                search_records = Entrez.read(handle)
                time.sleep(SLEEP_TIME)
            except:
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

                if i == ( REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str( REQUEST_LIMIT) + " attempts")
        
        #Download all of the protein sequences from each of the nucelotide records associated with this assembly and write them to a single fasta file to be used to compile the database. 
        #Stores all the sequence objected for all the proteins in the species
        sequences_to_write = []

        #Fetches the records
        for i in range(REQUEST_LIMIT):

            try:
                handle = Entrez.efetch(db='nuccore', id=search_records['IdList'], start='begin', stop='end', rettype='fasta_cds_aa')
                records = list(SeqIO.parse(handle, 'fasta'))
                time.sleep(SLEEP_TIME)
                break
            except:
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

                if i == ( REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str( REQUEST_LIMIT) + " attempts")

        sequences_to_write.extend(records)
        
        if len(sequences_to_write) > 0:
            #Write all the sequences to the output file
            print('\tWriting all CDS features to file for ' + reference_assembly_accession)
            SeqIO.write(sequences_to_write, reference_total_protein.format(ref_assembly_accession=reference_assembly_accession), 'fasta')
        else:
            print('\tNo sequences found for ' + reference_assembly_accession)
        

        #Make the blast database to be used for the reverse blast searches
        cmd = 'makeblastdb -in {input_file} -out {output_file} -dbtype prot'

        #Set up formatted file directories
        input_file = reference_total_protein.format(ref_assembly_accession=reference_assembly_accession)
        output_file = reference_blast_db.format(ref_assembly_accession=reference_assembly_accession)

        os.system(cmd.format(input_file=input_file, output_file=output_file))

def check_reverse_blast(query_accession, annoted_hit):
    '''
    Performs a local BLAST search of the hit against the query to see if the same feature is returned. 

    Parameters
    ----------
    query_accession: string
        The accession for the query gene that the hit resulted from.
    annotated_hit: AnnotatedHit object
        The hit that is being tested. 
    '''

    global reference_blast_db
    global reference_assembly_accession
    global reference_genome_accession

    #Template command
    cmd = 'blastp -query {q_path} -out {out_path} -db {db_path} -outfmt 5'

    #The temp files for this BLAST search
    temp_out = './reverse_blast/temp_out'
    temp_in = './reverse_blast/temp_in.fasta'

    #Clean the alignment sequence of the '-' from the gaps. This avoids a warning when running the local BLAST
    purged_alignment_seq = ''

    for c in annoted_hit.alignment_seq:
        if not c == '-':
            purged_alignment_seq = purged_alignment_seq + c

    #Write the sequence of the hit sequence to the temp_in file
    SeqIO.write(SeqRecord(Seq.Seq(purged_alignment_seq), id='temp'), temp_in, 'fasta')

    #The path to the database
    db_path = reference_blast_db.format(ref_assembly_accession=reference_assembly_accession)

    #Conduct the BLAST search
    os.system(cmd.format(q_path=temp_in, out_path=temp_out, db_path=db_path))

    #Get the BLAST results from the temp_out file and determine if there was a hit back to the query
    hit_returned_query = False
    with open(temp_out) as file:

        #Parse the output file from the BLAST search
        result_handle = list(NCBIXML.parse(file))

        if len(result_handle[0].alignments) > 0:
            record  = result_handle[0].alignments[0]
            #tqdm.write(query_accession + '\t' + reference_genome_accession + '\t' + record.hit_def)
            if (query_accession in record.hit_def) and (reference_genome_accession in record.hit_def):
                hit_returned_query = True
                
        return hit_returned_query

def calc_operon_cons():
    '''
    The main function in operon_conserve_detect.py. 
    Pipeline:
        1. Get parameters from the input file
        2. Set up the output directory and copy the input file to it
        3. Conduct BLAST search for each of the individual input records and store them all in a list of hits
        4. Group the hits into their respective GenomeFragment objects
        5. For each GenomeFragment, detect the feature for each hit it contains and assemble the operons for each GenomeFragement
        6. Group the GenomeFragments into Species
        7. Filter the Species based of the species_percent_id_limit
            |-> This will remove any species where no BLAST results are greater than the species_percent_id_limit
        8. Calculate the score for each Species
        9. Output all the operons into a CSV file.
    '''

    global genome_frags
    global species
    global output_dir
    global run_id
    global use_reference_threshold
    global reverse_blast
    
    #Take note of the start time
    start_time = datetime.datetime.now()

    ##Load parameters
    input_file = ''
    if(len(sys.argv) >= 2):
        input_file = 'input/' + sys.argv[1]
        load_input_file(filename=input_file)
    else:
        input_file = 'input/' + INPUT_JSON
        load_input_file(input_file)
    
    #Set Entrez parameters
    Entrez.email = EMAIL
    Entrez.api_key = E_API

    ##Setup the output directory
    output_dir = output_dir.format(run_id=run_id)
    os.mkdir(output_dir)
    os.mkdir(output_dir + 'svg/')
    
    #Copy the input file to the ouput directory
    shutil.copy(input_file, output_dir)


    #Process reference 
    tqdm.write('Processing reference...')
    make_reference_blastdb()
    process_reference()
    if use_reference_threshold:
        tqdm.write('Getting reference intergenic distance...')
        get_reference_intergenic_distance()
        tqdm.write(str(intergenic_limit) + '\n')

    ##Conduct BLAST search

    #Holds the hits from all of the BLAST searches
    final_hits = []

    for record in input_records:
        record_hits = search_blast(
            input_records=[record],
            db=database,
            max_attempts=max_blast_attemts,
            search_mult_factor=blast_search_mult_factor, 
            min_cover=coverage_min, 
            max_hits=max_hits, 
            tax_incl=tax_include, 
            tax_excl=tax_exclude, 
            e_cutoff=e_val, 
            annotate=annotate, 
            extensive_search=extensive_search
        )

        #Filters based off the reverse BLAST check if the user requested it
        if(reverse_blast):
            purged_record_hits = []

            for hit in tqdm(record_hits, desc='Reverse Blast'):
                if check_reverse_blast(record, hit):
                    purged_record_hits.append(hit)

            final_hits.extend(purged_record_hits)
        else:
            final_hits.extend(record_hits)

        print('Running total of hits returned: ' + str(len(final_hits)))
        time.sleep(3)
    
    ##Group all the hits into GenomeFragments
    for hit in tqdm(final_hits, desc='Grouping Hits to Fragments'):

        #Check if hits are annotated
        if not isinstance(hit, AnnotatedHit):
            print("The annotate option needs to be set to true for the GenomeFragment grouping to work.")
            return
        
        #Check if there is a genome fragment for the hit
        genome_fragment_check = genome_fragment_exists(fragment_name=hit.genome_fragment_name, genome_accession=hit.genome_accession)

        #Adds the hit to appropriate GenomeFragment or makes a new GenomeFragment object if one does not exist
        if genome_fragment_check[0]:
            genome_fragment_check[1].add_hit(hit)
        else:
            new_genome_fragment = GenomeFragment(name=hit.genome_fragment_name, genome_fragment_accession=hit.genome_accession, req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME, cache_directory=cache_dir)
            new_genome_fragment.add_hit(hit)
            genome_frags.append(new_genome_fragment)


    #Get all the assembly accessions for the hits
    tqdm.write('Pre-processing fragments...')
    pre_process_fragments()


    #Grouping all GenomeFragments to species
    for fragment in tqdm(genome_frags, desc='Grouping Fragments to Species'):

        #Add the fragment to the appropriate Species object or make a new Species object
        species_check = species_exits(fragment.assembly_accession)

        if species_check[0]:
            species_check[1].add_genome_fragment(fragment)
        else:
            new_species = Species(assembly_accession=fragment.assembly_accession, genome_fragments=[fragment])
            species.append(new_species)
   
    '''
    ##Filter the species based off the species_percent_id_limit
    to_remove = []
    for sp in tqdm(species, desc='Filtering Species'):

        #Calculate the percent id for each of the querery genes
        tqdm.write("\tCalculating percent IDs for " + str(sp.assembly_accession) + " ...")
        sp.get_query_percent_ids(input_records)

        #Determine if any of the quereies are above the limit
        above_limit = False 

        for val in sp.query_percent_ids.values():
            if val >= species_percent_id_limit:
                above_limit = True
        
        #Remove species object if it is not above the limit
        if not above_limit:
            to_remove.append(sp)
    
    for sp in to_remove:
        species.remove(sp)'''


    ##Detect the features for the hits in each GenomeFragment, assemble the operons, and sort into a Species

    #How many threads to run at once
    thread_limit = 10

    #Holds threads for each iteration
    processing_threads = []

    #Proccess one species at a time
    to_remove = []
    for sp in tqdm(species, desc='Species Processing'):

        #Iterate through each fragment in the species and process it
        for fragment in sp.genome_fragments:
            lock = threading.Lock()
            temp_thread = threading.Thread(target=process_frag, args=(fragment, lock))
            temp_thread.start()
            processing_threads.append(temp_thread)

            if len(processing_threads) == thread_limit:
                #Waits for all threads to finish before continuing 
                for t in processing_threads:
                    t.join()
                del processing_threads
                processing_threads = []
            
        #Waits for all threads to finish before continuing 
        for t in processing_threads:
            t.join()

        #Calculate the percent id for each of the hits in every fragment in every species 
        tqdm.write('Caclulating percent ids for ' + str(sp.assembly_accession) + '...')
        calculate_percent_ids(sp)
        sp.get_query_percent_ids(input_records)

        #Determine if any of the quereies are above the limit
        above_limit = False 

        for val in sp.query_percent_ids.values():
            if val >= species_percent_id_limit:
                above_limit = True
        
        if above_limit:

            #Extract the species name
            tqdm.write("Extracting species name for " + str(sp.assembly_accession) + " ...")
            sp.extract_names()
            tqdm.write('\t' + str(sp.species_name))
            
            #Extract the taxIds
            tqdm.write("Extracting txid for " + str(sp.species_name) + " ...")
            sp.extract_taxids()
            tqdm.write('\t' + str(sp.taxid))

            #Extract the genome accessions
            tqdm.write("Extracting genome accessions for " + str(sp.species_name) + " ...")
            sp.extract_genome_accessions()
            tqdm.write('\t' + str(sp.genome_fragments_accessions))

            #Calculate the score for the species
            tqdm.write("Calculating structural similarity score for " + str(sp.species_name) + " ...")
            sp.measure_sim(input_records)
            tqdm.write('\t' + str(sp.sim_score))

            #Drawing operons
            tqdm.write("Drawing " + str(sp.species_name) + " ...")
            sp.draw_figure(color_code=color_code, output_dir=output_dir + 'svg/')

            #Clear up the species object so it is not using up memory
            tqdm.write("Cleaning up " + str(sp.species_name) + " ...")
            sp.clean()
        else:
            #Clear up the species object so it is not using up memory
            tqdm.write("Did not meet percent ID threshold. Cleaning up " + str(sp.assembly_accession) + " ...")
            sp.clean()

            to_remove.append(sp)



    ##Output the results
    '''
    output_filename_base = ''
    if sys.argv[2]:
        output_filename_base = sys.argv[2]
    else:
        output_filename_base = 'output_' +str(datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S"))'''

    
    #Remove all species that did not meet the percent ID requirement
    for sp in to_remove:
        species.remove(sp)

    print("Writing output now...")
    write_all_out()

    #Get time elapsed
    end_time = datetime.datetime.now()
    elapsed_time = end_time - start_time
    print("*"*10 + " Finished in " + str(elapsed_time) + "*"*10)



calc_operon_cons()