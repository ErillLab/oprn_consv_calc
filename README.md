# oprn_consv_calc

### Overview

The operon_conserve_detect.py script takes a reference operon as a list of genome accessions and will determine its structural conservation (i.e. gene order) within a specific taxonomic clade.

### Pipeline
![Workflow](/extra/operon_detect_pipeline.svg)

### Object Definitions

* **GenomeFeature (feature.py)**
	* Holds information about a feature
	* Member variables:
		* Genome accession 
		* Coding start and end
		* Strand (+ or -)
		* Protein accession
		* Locus tag
		* Amino acid sequence
	* Also requires the following parameters to be passed in when initializing the object:
		* Request limit for the Entrez searches
		* Sleep time between consecutive Entrez searches
	* Main Functions:
		* get_intergenic_distance(self, other): Determines the distance to another GenomeFeature
		* \_\_str\_\_() and \_\_eq\_\_() are overridden string and equivalence functions


* **AnnotatedHit (feature.py)**
	* A subclass of GenomeFeature. Contains additional information pertaining to the BLAST hit results
	* Member variables:
		* The query accession for the BLAST search that resulted in this hit
		* The alignment start and end positions
		* The percent amino acid identity
		* The alignment sequence
	* Main Functions:
		* fetch_feature(self, record, margin_limit, max_attempts, mult_factor): Determines the feature that corresponds to this hit. It works by finding a feature in the full genome accession record that contains the alignment start and end positions in bounds with a set margin. 
		* \_\_str\_\_() and \_\_eq\_\_() are overriden string and equivalence functions

* **Operon (operon.py)**
	* Holds a set of GenomeFeature (and/or by extension AnnotatedHit) objects that have been assigned to the same operon.
	* Member variables:
		* features[]: A list of GenomeFeature (and/or by extension AnnotatedHit) objects
		* genome_accession: The nucleotide accession for the GenomeFragment this operon belongs to
		* genome_features[]: A list of all the features in the GenomeFragment this operon belongs to.
		* strand: Whether the operon is on the plus or minus strand
	* Main Functions:
		* add_feature(self, feature): Checks if the feature belongs to the GenomeFragment this operon belongs to, and adds it to the list of features. It then sorts the current list of features from 5' to 3' 
		* \_\_str\_\_() is a overridden string function

* **GenomeFragment (genome_fragment.py)**
	* Represents a single nucleotide accession (i.e. plasmid, chromosome, contig, ...etc.)
	* Member variables:
		* hits[]: The AnnotatedHits associated with this GenomeFragment
		* all_features[]: All the GenomeFeatures that are parsed out of the full record for this GenomeFragment
		* operons[]: Operon objects associated with this GenomeFragment
		* name: The name of this GenomeFragment
		* genome_accession: The genome accession for this GenomeFragment
		* assembly_accession: The accession for the assembly this nucleotide record belongs to
		* taxid: the taxonmic id
		* full_record: The entire nucleotide record. Downloaded locally in the cache directory.
	* Main Functions:
		* fetch_features(self): Obtains all coding features in the genome and saves them as a list of GenomeFeature obejcts
		* fetch_record(self): Obtains the full genome record for this nucleotide accession. Reads in from a file or downloads it.
		* fetch_hit_features(self, margin_limit, max_attempts, mult_factor): Fetches the features associated with each of the hits associated with this GenomeFragment
		* purge_hits(self): Removes any duplicate hits
		* add_hit(self, a_hit): Adds an AnnotatedHit object to the current list of hits.
		* assemble_operons(self, feature_limit, intergenic_limit): Takes the list of hits and organizes them into putative operons (Operon objects)


* **Species (species.py)**
	* Holds all the GenomeFragment objects that belong to the same species
	* Member variables:
		* assembly_accession: The assembly accession for this Species
		* species_name: The name of the assembly
		* genome_fragments[]: A list of GenomeFragment objects associated with this Species
		* genome_fragment_accessions[]: A list of the accessions for all the GenomeFragments in genome_fragments[]
		* sim_score: The structural similarity score for this Species relative to the reference
		* query_percent_ids[]: The amino acid percent identity for each of the reference genes.
		* taxid: The taxonomic id for this assembly
	* Main Functions
		* add_genome_fragment(self, genome_fragment): Adds a GenomeFragment object to the species
		* measure_sim(self, reference_operon): Calculates the structural similarity of the Species relative to the reference operon. 
		* extract_taxids(self), extract_genome_accessions(self), extract_names(self): These functions set the respective member variables so that the information is stored when the Species is cleaned from memory. 
		* draw_figure(self, color_code, output_dir): Draw a genome diagram of each of the Operons to SVG files. 
		* clean(self): Deletes all the GenomeFragments to clear up memory. 


### Sample Input

See the template input JSON file in /input/. 

Parameter | Description 
---|---
*entrez* Parameters |
`request_limit` | The max number of attempts that should be made to complete an Entrez request
`sleep_time` | The amount of time that should be waited between consective attempts at an Entrez request (ms)
`email` | Email to use for the Entrez api
`api_key` | Entrez API key
---|---
*blast* Parameters|
`tax_include` and `tax_exclude` | These will restrict the BLAST to include (`tax_include`) and exclude (`tax_exclude`) specific taxa. 
`database` | The database to conduct the BLAST search against.
`blast_type` | Can be set to `remote` or `local`. If `local` , use `local_db_path` to set the path. 
`e-val` | The maximum e-value threshold for the BLAST searches. 
`coverage_min` | The coverage threshold to restrict the BLAST results. 
`max_hits` | The max hits to recieve from a BLAST search. 
`max_attempts` | The maximum number of attempts that should be made to pull all the BLAST hits within the e-value cutoff. This is only applicable is the `extensive_search` option is set to `true`.
`search_mult_factor` | The `max_hits` is increased by this factor through each iteration if the `extensive_search` option is set to `true`. 
`annotate` | Can be `true` or `false`. Determines if the BLAST search will return `AnnotatedHit` objects. Current functionality requires that it MUST BE `true`. 
`extensive_search` | If `true`, multiple BLAST searches will be conducted (if neccessary) to pull all hits within the `e-value` limit. If `false`, a single BLAST search will be conducted regardless of whether all hits witin the `e-value` threshold are returned. 
`reverse_blast` | If `true`, all return hits will be tested for true homology by BLAST searching the reference genome with the hit. 
---|---
*hit\_feature\_detection* |
`margin_limit` | The acceptable margin between the alignment positions and the positions of the annotated feature. 
`max_attempts`| The maximum number of attempts that will be made to identify a feature for a hit.
`mult_factor` | The factor by which the `margin_limit` will increase for every attempt at identifying the feature for a hit. 
---|---
*operon\_assembly* |
`feature_limit` | Maximum number of features allowed between two hits on an operon.
`intergenic_limit` | Maximum distance in bp allowed between two genes in an operon.
`use_ref_limit` and `ref_limit_margin` | If `true`, the `intergenic_limit` will be determined by the max intergenic distance in the refernce operon multiplied by the margin. 
---|---
`thread_limit` | A multithreaded approach is implemented when processing each of the Species objects. This is the maximum number of threads.
`species_percent_id_limit` | For a species to be outputted, it must have at least one hit with an amino acid percent identity above this limit.
`input_records` | The protein accessions for the genes in an operon. The order the genes are representative of the order of the genes in the reference genome. 
`reference_genome_accession`| The nucleotide accession of the nucleotide record from the reference speices. This nucleotide record contains the reference operon. 
`reference_genome_assembly`| The genome assembly accession from which the reference operon comes from. Needed to conduct the reverse BLAST locally. 
`reference_genome_name` | The name of the reference assembly.
`cache_dir` | The cache directory where all nucleotide records will be downloaded locally. 
`color_code` | A dictionary that defines the colors for each of the genes in the reference operon. These colors are used to draw the opeon diagram for the putative operons. `intergenic` refers to features that are inserted into the operon, but were not a hit to any of the reference genes. 


### Output

Results are outputted to a CSV and each operon is drawn to a SVG: 

### Usage
Once the `conda` environment is set up, the script can be ran via command line: <br><br>
`python operon_conserve_detect.py [input file name].json` <br> <br>
where the input file name is the name of the input file, not including the directory. The input file is expected to be in `/input/`. 


















