
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

Azurify aggregates data from CiVIC, ClinVar, gnomAD, COSMIC, KEGG, Pubmed and Uniprot to create a feature set that is combined with over 15,000 clinical classifications to create a model that can classify small variants (SNVs & Indels < 50bp).
The output classes being pathogenic, Likely pathogenic, uncertain significance (VUS), likely benign, and benign. 


    Installation
    Usage
    Documentation
    Contributing
    License

## Installation

Installation of Azurify and its dependencies are made easy and can be found within the setup.py file. Any dependencies associated with model and figure generation are outside of Azurify and will need to be installed manually.
```
git clone https://github.com/faryabiLab/Azurify.git
cd Azurify
python setup.py
```
bash

pip install -r requirements.txt

Usage

Explain how to use the project. Include examples and sample code snippets if necessary.

python

python your_script.py

Documentation

Link to or include documentation for the project. This can be in the form of inline comments in the code, a separate documentation file, or a link to an external documentation site.
Contributing

## A note to Users

In addition to Azurify, we have provided all of the code used to develop as well as evaluate the Azurify model in an effort to ease the accesiblity of variant classification across the field. 

As the Azurify project expands we hope to include more resources as features, automate model generation to keep pace with the field and generate disease specific models. 
We would also love to include training data beyond the borders of our instituion and seek to provide a framework where variant classification models can be shared to mitigate variability, so collaborators and contributors are encouraged to reach out. 

## License
Azurify classifies the pathogencity of small genomic variants using 
predictive machine learning through clinically relevant features.
Copyright (C) 2024  Ashkan Bigdeli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU General Public License](https://www.gnu.org/licenses/) for more details.

## Author
Please direct questions and bug reports to Ashkan Bigdeli ashkan.bigdeli@pennmedicine.upenn.edu
