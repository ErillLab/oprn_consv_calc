# oprn_consv_calc

### Overview

The operon_conserve_detect.py script takes a reference operon as a list of genome accessions and will determine its structural conservation (i.e. gene order) within a specific taxonomic clade.

### Pipeline
![Workflow](/extra/operon_detect_pipeline.svg)

### Object Definitions

* **GenomeFeature**
	* Holds information about a feature including:
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
		* __str__() and __eq__() are overriden string and equivalence functions


* **AnnotatedHit**

* **Operon**

* **GenomeFragment**

* **Species**

### Sample Input

### Output

### Usage

