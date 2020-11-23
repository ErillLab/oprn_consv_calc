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

### Sample Input

### Output

### Usage

