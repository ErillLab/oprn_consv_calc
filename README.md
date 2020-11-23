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
	* Functions:
		* get_intergenic_distance(self, other): Determines the distance to another GenomeFeature
		* \_\_str\_\_() and \_\_eq\_\_() are overridden string and equivalence functions


* **AnnotatedHit (feature.py)**
	* A subclass of GenomeFeature. Contains additional information pertaining to the BLAST hit results
	* Member variables:
		* The query accession for the BLAST search that resulted in this hit
		* The alignment start and end positions
		* The percent amino acid identity
		* The alignment sequence
	* Functions:
		* fetch_feature(self, record, margin_limit, max_attempts, mult_factor): Determines the feature that corresponds to this hit. It works by finding a feature in the full genome accession record that contains the alignment start and end positions in bounds with a set margin. 
		* \_\_str\_\_() and \_\_eq\_\_() are overriden string and equivalence functions

* **Operon**
	* Holds a set of GenomeFeature (and/or by extension AnnotatedHit) objects that have been assigned to the same operon.
	* Member variables:
		* features[]: A list of GenomeFeature (and/or by extension AnnotatedHit) objects
		* genome_accession: The nucleotide accession for the GenomeFragment this operon belongs to
		* genome_features[]: A list of all the features in the GenomeFragment this operon belongs to.
		* strand: Whether the operon is on the plus or minus strand
	* Functions:
		* add_feature(self, feature): Checks if the feature belongs to the GenomeFragment this operon belongs to, and adds it to the list of features. It then sorts the current list of features from 5' to 3' 
		* \_\_str\_\_() is a overridden string function

* **GenomeFragment**
	* Represents a single nucleotide accession (i.e. plasmid, chromosome, contig, ...etc.)
	* Member variables:
		* 

* **Species**

### Sample Input

### Output

### Usage

