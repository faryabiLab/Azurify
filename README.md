
```
  ___                     _   __        
 / _ \                   (_) / _|       
/ /_\ \ ____ _   _  _ __  _ | |_  _   _ 
|  _  ||_  /| | | || '__|| ||  _|| | | |
| | | | / / | |_| || |   | || |  | |_| |
\_| |_//___| \__,_||_|   |_||_|   \__, |
                                   __/ |
                                  |___/ 
```

Azurify aims to classify the pathogencity of small genomic variants by leveraging machine learning on a feature set of resources that can be used in the clinical classification of somatic variants for the purpose of cancer precision mediicine. 

Azurify aggregates data from CiVIC, ClinVar, gnomAD, COSMIC, KEGG, PubMed, Uniprot and over 15,000 clinical classifications to create a model that can determinee the pathogencity of small genomic variants (SNVs & Indels < 50bp).
The output classes being pathogenic, Likely pathogenic, uncertain significance (VUS), likely benign, and benign. 

## Installation

Installation of Azurify and its dependencies are made easy and can be found within the setup.py file. Any dependencies associated with model and figure generation are outside of Azurify and will need to be installed manually.
```
git clone https://github.com/faryabiLab/Azurify.git
cd Azurify
python setup.py
```
## Usage

```
python azurify.py -i /path/to/input.tsv -o /path/to/output.tsv
```
Azurify expects the following columns as input: CHROM, POS, REF, ALT, FAF, GENE, PCHANGE, EFFECT, EXON_Rank. All values can be derived when annotating a VCF with [snpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/). VCF and hg38 support are scheduled for the next release. Currently Azurify works using hg19 locations.

Example Input:

|CHROM|POS|REF|ALT|FAF|GENE|PCHANGE|EFFECT|EXON_Rank|
|:----|:----|:----|:----|:----|:----|:----|:----|:----|
|chr7|66993288|C|T|50.85714286|SBDS|p.Val130Met|missense_variant|3|
|chr4|1961074|G|A|41.69453735|WHSC1|p.Glu1099Lys|missense_variant|18|
|chr12|57102878|T|C|32.51928021|STAT6|p.Asp419Gly|missense_variant|12|


You may include additional columns to the rquired tab-delimtied input and they will be appended to your final results.

## Runtime

Azurify annotates 100,000 variants in approximately 30 minutes of runtime. You should expect to need more than 16GB of memory and longer runtimes as record volume increases.

## The Azurify Project

In addition to Azurify, we have provided all of the code used to develop as well as evaluate the Azurify model in an effort to ease the accesiblity of variant classification across the field. 

As the Azurify project expands we hope to include more resources as features, automate model generation to keep pace with the field and generate disease specific models. 
We would also love to include training data beyond the borders of our instituion and seek to provide a framework where variant classification models can be shared to mitigate variability, so collaborators and contributors are encouraged to reach out. 

## License
Azurify classifies the pathogencity of small genomic variants using 
predictive machine learning through clinically relevant features.
Copyright (C) 2024 Ashkan Bigdeli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU General Public License](https://www.gnu.org/licenses/) for more details.

## Author
Please direct inqueries, questions and bug reports to Ashkan Bigdeli ashkan.bigdeli@pennmedicine.upenn.edu
