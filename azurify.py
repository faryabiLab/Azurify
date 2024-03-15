#I need domain, cosmic, clinvar, civic, gnomad, litvar
#the user should input a file that has a sample string if needed, chrom, start, ref, alt, effect, exon#, 


#hardcoded files to remove from git
domain_file_path = "C:\\dev\\phd\\dt\\resource_data\\uniprot_hg19_domain.txt"
input_file_path = 'C:\\dev\\phd\\cornell_data\\cornell_withexons_parsed.tsv'


#imports
import pandas as pd
import argparse
import catboost
import os
                                                                                                                                                                                                                                                           
def add_domain(df):
    domains_df = pd.read_csv(domain_file_path, sep='\t', header=None, names=['chrom', 'start', 'stop', 'name'], dtype={'chrom': str, 'start': int, 'stop': int, 'name': str})
    result_df = pd.DataFrame(columns=['chrom', 'pos', 'DomainName'])

                                                                                                                                                                                                                                                              
def merge_keys(df):
    for filename in os.listdir('features/'):
        k = pd.read_csv(('C:\\dev\\phd\\dt\\resource_data\\resdev\\keys\\' + filename), sep='\t',low_memory=False)
        k.set_index("KEY", inplace = True)
        df = ex.join(k)
    return(df)

def print_ascii():
    print(r"""
  ___                     _   __        
 / _ \                   (_) / _|       
/ /_\ \ ____ _   _  _ __  _ | |_  _   _ 
|  _  ||_  /| | | || '__|| ||  _|| | | |
| | | | / / | |_| || |   | || |  | |_| |
\_| |_//___| \__,_||_|   |_||_|   \__, |
                                   __/ |
                                  |___/ 
          """)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Azurify classifies the pathogencity of small variants based on clinical training labels.")
    parser.add_argument("-s", "--sample_column_name", required=False, metavar="COLUMN_NAME", type=str, help="Specify the sample column name.")
    parser.add_argument("-l", "--litvar", action="store_true", required=False, help="Should we ping LitVar to pull publications. Takes 1 second per ping")
    parser.add_argument("-i", "--input_file", metavar="FILE_PATH", required=False, type=str, help="Specify the input file.")
    parser.add_argument("-o", "--output_directory", metavar="DIR_PATH", required=False,type=str, help="Specify the output directory.")
    parser.add_argument("-f", "--output_filename", metavar="OUTPUT_FILENAME", required=False,type=str, help="Specify the output filename.")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Check for required arguments
    if not args.input_file or not args.output_directory or not args.output_filename:
        parser.error("Input file, output directory, and output filename are required.")

    return args

def main():
    # Example usage
    print_ascii()
    args = parse_arguments()
    print("hello")
    # Access the values of the arguments
    #print("Sample Column Name:", args.sample_column_name)
    #print("Include LitVar:", args.litvar)
    #print("Input File:", args.input_file)
    #print("Output Directory:", args.output_directory)
    #print("Output Filename:", args.output_filename)

if __name__ == '__main__':
    main()




import time
#pchange = []
#pmidc = []
out="C:\\dev\\phd\\dt\\litvar_pchanges.tab"
with open(out, 'w+') as o:
    for i in pchanges:
        url="https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1/entity/search/" + i
        time.sleep(.5)
        response = requests.get(url)
        code = response.status_code
        if code == 200:
            jr = response.json()
            pmid_count = jr[0]['pmids_count']
        else:
            pmid_count = 0
        
        o.write(i + "," + str(pmid_count) + '\n')
    #pchange.append(i)
    #pmidc.append(pmid_count)
#df2 = pd.DataFrame({'PCCHANGE': pchange, 'PMID_COUNT':pmidc})
print("thats a win")



#parse the domain file to meet the below specs

domain_file_path = "C:\\dev\\phd\\dt\\resource_data\\uniprot_hg19_domain.txt"
input_file_path = 'C:\\dev\\phd\\cornell_data\\cornell_withexons_parsed.tsv'

# Read domain data into a DataFrame
#domains_df = pd.read_csv(domain_file_path, sep='\t', header=None, names=['chrom', 'start', 'stop', 'name'], dtype={'chrom': str, 'start': int, 'stop': int, 'name': str})

# Create a DataFrame to store the results
#result_df = pd.DataFrame(columns=['chrom', 'pos', 'DomainName'])

with open(input_file_path, 'r') as input_file:
    header = next(input_file).rstrip('\n')

    for line in input_file:
        fields = line.rstrip('\n').split("\t")
        chrom, pos = fields[0], int(fields[1])

        # Merge domain information based on the chromosome
        merged_df = pd.merge_asof(result_df, domains_df[domains_df['chrom'] == chrom].sort_values('start'),
                                  left_on='pos', right_on='start', direction='forward')

        # Check if the position is within the domain range
        mask = (pos >= merged_df['start']) & (pos <= merged_df['stop'])

        # Update the 'DomainName' column with the domain name or 'NA'
        merged_df.loc[mask, 'DomainName'] = merged_df.loc[mask, 'name']

        # Append the current line's data to the result DataFrame
        result_df = pd.concat([result_df, pd.DataFrame({'chrom': [chrom], 'pos': [pos], 'DomainName': ['NA']})])

# Print or further process the result DataFrame                                                                                                                                                                                                                                                                                                                                                                                                                
print(result_df)