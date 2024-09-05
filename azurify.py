import os
import argparse
import subprocess
from collections import OrderedDict
from tqdm import tqdm
import pandas as pd
from catboost import CatBoostClassifier
from liftover import get_lifter


#global get working directory
cwd = os.getcwd()

#global progress bar
progress_bar = tqdm(total=100, desc="Progress", unit="iter")


def split_df(df, in_cols):
    included_df = df[in_cols]
    excluded_columns = [col for col in df.columns if col not in in_cols]
    excluded_df = df[excluded_columns]
    return included_df, excluded_df

def convert38(input_file, output_file):
    converter = get_lifter('hg38', 'hg19', one_based=True)
    with open(input_file, 'r') as f:
        lines = f.readlines()

    processed_lines = [lines[0]] 
    for line in lines[1:]:
        values = line.strip().split('\t')
        chrom, pos = values[:2]
        try:
            x = converter.query(chrom, int(pos))[0]
        except IndexError:
            continue
        values[0], values[1] = x[0], str(x[1])
        processed_lines.append('\t'.join(values))

    with open(output_file, 'w') as f:
        for line in processed_lines:
            f.write(line + '\n')

    return open(output_file)

def run_snpeff(snpeff_jar_path, input_file, output_file): 
    snpeff_command = f"java -jar {snpeff_jar_path} -canon hg19 {input_file} > {output_file}"
    
    subprocess.run(snpeff_command, shell=True, check=True)
    
    return output_file

def format_annotations(input_file, output_file):
    with open(input_file, 'r') as anno_in, open(output_file, 'w') as anno_out:
        annod = OrderedDict.fromkeys(['CHROM', 'POS', 'REF', 'ALT', 'FAF', 'GENE', 'PCHANGE', 'EFFECT', 'EXON_Rank'])
        anno_out.write('\t'.join(annod.keys()) + '\n')

        for line in anno_in:
            if line.startswith('##') or line.startswith('#'):
                continue
            values = line.strip().split('\t')
            annod['CHROM'] = 'chr' + str(values[0])
            annod['POS'] = values[1]
            annod['REF'] = values[3]
            annod['ALT'] = values[4]
            annod['FAF'] = values[5]

            sample_ann = values[7]
            anns = sample_ann.split(';')[1].split('|')
            annod['EFFECT'] = anns[1]
            annod['GENE'] = anns[3]
            annod['EXON_Rank'] = anns[8].split('/')[0]
            annod['PCHANGE'] = anns[10]

            anno_out.write('\t'.join(annod.values()) + '\n')
    
    return output_file

def add_domain(df):
    domain_file_path = os.path.join(cwd, 'utils/uniprot_hg19_domain_parsed.txt')
    domains_df = pd.read_csv(domain_file_path, sep='\t', header=None, names=['CHROM', 'START', 'STOP', 'Domain'], dtype={'CHROM': str, 'START': int, 'STOP': int, 'Domain': str})

    # # Merge based on chromosome and then intersect loci and update dataframe with the appropriate domain
    merged_df = pd.merge(df, domains_df, on='CHROM', how='left')
    mask = (merged_df['POS'] >= merged_df['START']) & (merged_df['POS'] <= merged_df['STOP'])
    df.loc[mask, 'Domain'] = merged_df.loc[mask, 'Domain']

    return(df)

def add_litvar(df):
    litvar_file_path = os.path.join(cwd, 'utils/litvar_pchanges.txt')
    litvar_df = pd.read_csv(litvar_file_path, sep='\t', header=None, names=['PCHANGE', 'PMID_COUNT'], dtype={'PCHANGE': str, 'PMID_COUNT': int})

    merged_df = pd.merge(df, litvar_df, on="PCHANGE", how='left')
    return(merged_df)

def add_kegg(df):
    kegg_file_path = os.path.join(cwd, 'utils/kegg_pathways.txt')
    kegg_df = pd.read_csv(kegg_file_path, sep='\t', header=None, names=['GENE'], dtype={'PCHANGE': str})

    kegg_df['IN_KEGG'] = 1
    plus_kegg = pd.merge(df,kegg_df[['GENE', 'IN_KEGG']],on='GENE', how='left')
    plus_kegg["KEGG"] = plus_kegg['IN_KEGG'].isnull()*1
    plus_kegg = plus_kegg.drop(['IN_KEGG'], axis=1)
    return(plus_kegg)

def merge_keys(df):
    utilp = os.path.join(cwd, 'utils')

    df['KEY'] = df['CHROM'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']

    for filename in os.listdir(utilp):
        if filename.endswith('_key.tsv.gz'):
            k = pd.read_csv((os.path.join(utilp, filename)), sep='\t',low_memory=False)
            df = df.merge(k, on='KEY', how='left')
            progress_bar.update(10)
    return(df)

def run_azurify(df, drug_targets):
    azurify = CatBoostClassifier()
    if drug_targets:
        azurify.load_model((os.path.join(cwd,'models/azurify_hg19.0.99.json')), format='json')
    else:
        azurify.load_model((os.path.join(cwd,'models/azurify_no_drug_hg19.0.99.json')), format='json')
    
    template_az = pd.read_csv((os.path.join(cwd ,'models/azurify_template.txt')), sep='\t')
    df = df[template_az.columns]
    df = df.infer_objects(copy=False).fillna('-999')

    pred = azurify.predict(data=df)[:,0]
    prob = azurify.predict_proba(X=df)
    return(pred, prob)

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
    print("\nVersion 0.9.9")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Azurify classifies the pathogencity of small variants based on clinical training labels.")
    parser.add_argument("-l", "--litvar", action="store_true", required=False, help="Should we ping LitVar to pull publications. Takes 1 second per ping")
    parser.add_argument("-i", "--input_file", metavar="FILE_PATH", required=False, type=str, help="Specify the input file.")
    parser.add_argument("-o", "--output_filename", metavar="OUTPUT_FILENAME", required=False,type=str, help="Specify the output filename.")
    parser.add_argument("-g", "--geneom_build", metavar="GENOME_BUILD", required=False, type=str, help="Specify the input genome, i.e hg19 or hg38.")
    parser.add_argument("-s", "--snpeff_jar_path", metavar="SNPEFF_JAR_PATH", required=True, type=str, help="Specify the path to the snpeff jar file.")
    parser.add_argument("-d", "--no_drug_targets", metavar="Drug Targets", required=False, type=str, help="Omit the use of drug targets in the model.")  
    # Parse the command-line arguments
    args = parser.parse_args()

    # Check for required arguments
    if not args.input_file or not args.output_filename:
        parser.error("Input file and Output file are required. Ex. Usage: python azurify.py -i /path/input.tsv -o /path/output.ts")
    return args

def main():
    print_ascii()
    progress_bar.update(1)

    #grab args and files
    args = parse_arguments()
    out_file = args.output_filename
    uin = args.input_file
    udf = pd.read_csv(uin, sep='\t',low_memory=False)
    progress_bar.update(5)

    #get output directory for intermediate files
    directory = os.path.dirname(out_file)
    filename = os.path.basename(out_file)
    name, ext = os.path.splitext(filename)

    #convert to hg19 if required
    if args.geneom_build == 'hg38':
        o19 = f"{name}{"hg19"}{ext}"
        uin = convert38(uin, o19)
        progress_bar.update(7)

    #run snpeff
    if args.snpeff_jar_path:
        so = f"{name}{"snpeff"}{ext}"
        az_in = run_snpeff(args.snpeff_jar_path, o19, so)
        progress_bar.update(9)

    
    udf = pd.read_csv(az_in, sep='\t',low_memory=False)


    #split the df into columns used/not used by model
    in_cols=['CHROM','POS','REF','ALT','FAF','GENE','PCHANGE','EFFECT', 'EXON_Rank']
    adf, xdf = split_df(udf, in_cols)

    #add features
    ddf = add_domain(adf)
    progress_bar.update(10)
    kdf = add_kegg(ddf)
    progress_bar.update(10)
    ldf = add_litvar(kdf)
    progress_bar.update(10)

    #add keyed feeatures
    mdf = merge_keys(ldf)


    #run model
    if args.no_drug_targets:
        preds, prob = run_azurify(mdf, False)
    else:
        preds, prob = run_azurify(mdf, True)

    mdf['Pathogenicity'] = preds
    mdf['BP'] = prob[:,0]
    mdf['PP'] = prob[:,1]
    mdf['LBP'] = prob[:,2]
    mdf['LPP'] = prob[:,3]
    mdf['VP'] = prob[:,4]

    final = pd.concat([mdf, xdf], axis=1)

    final = final.drop(['KEY'], axis=1)

    #replace to standard ACMG nomenclature
    final['Pathogenicity'] = final['Pathogenicity'].astype(str).str.replace('Disease Associated', 'Pathogenic')
    final['Pathogenicity'] = final['Pathogenicity'].astype(str).str.replace('VOUS', 'VUS')
    final['Pathogenicity'] = final['Pathogenicity'].astype(str).str.replace('Probably DA', 'Likely Pathogenic')

    final.to_csv(out_file, sep="\t", index=False)
    progress_bar.update(19)
    progress_bar.close()

    print("Classifications Complete, results written to: " + out_file)

if __name__ == '__main__':
    main()