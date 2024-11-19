
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

Azurify aims to classify the pathogencity of small genomic variants by leveraging the power of machine learning to operate on a feature-set derived from data sources recommended by professional societies and clinically annotated somatic variant datasets.

Azurify aggregates data from CiVIC, ClinVar, gnomAD, COSMIC, KEGG, PubMed, Uniprot and over 15,000 clinical classifications to create a model that can determine the pathogencity of small genomic variants (SNVs & Indels < 50bp).
Azurify outputs variant classes and probabilities for pathogenic, likely pathogenic, uncertain significance (VUS), likely benign, and benign variants.

## Installation

Azurify is written in python3 and can easily be installed via pip and git. Azurify is also dependent on a valid [snpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/) installation with hg19 configured. Any dependencies associated with the Azurify publication may need to be installed manually and are clearly marked at the top of all corresponding notebooks. 
```
git clone https://github.com/faryabiLab/Azurify.git
cd Azurify
pip install -r requires.txt
```
## Usage

```
python3 azurify.py -i hg38_example.txt -o /path/results/ -s /path/snpEff/snpEff.jar -g hg38
```
Azurify expects the following columns as input: CHROM, START, STOP, REF, ALT, VAF in tab-delimited format.  The current version of Azurify is based on hg19, but accepts hg38 positions via the -g hg38 parameter. 

Example Input:

|CHROM|START|STOP|REF|ALT|VAF|
|:----|:----|:----|:----|:----|:----|
|chr7|66993288|66993289|C|T|50.85|
|chr4|1961074|1961075|G|A|7.69|
|chr12|57102878|57102879|T|C|42.51|


You may include additional columns to the required tab-delimited input and they will be appended to your final results.

## Creating your own model

To create your own model using the Azurify feature set, simply run Azurify and follow the Juypter Notebook located within the repo under publication_code/build_model.iypnb. Please make sure your internal classifications are under the defined column "CATEGORIZATION". 

## Runtime

Azurify annotates 100,000 variants in approximately 30 minutes of runtime. You should expect to need more than 16GB of memory and longer runtimes as record volume increases. Considering chunking your input into smaller files if you are on a low memory machine. 

## The Azurify Project

In addition to Azurify, we have provided all of the code used to develop as well as evaluate the Azurify model should you want to generate your own data using our methodologies. 

As the Azurify project expands we hope to include more resources as features, automate model generation to keep pace with the field, and generate disease specific models.

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
