import os
import argparse
import subprocess
import glob
from collections import OrderedDict
from tqdm import tqdm
import pandas as pd
from catboost import CatBoostClassifier
from liftover import get_lifter


#global get working directory
cwd = os.getcwd()



def split_df(df, in_cols):
    included_df = df[in_cols]
    excluded_columns = [col for col in df.columns if col not in in_cols]
    excluded_df = df[excluded_columns]
    return included_df, excluded_df

def convert38(input_file, output_file):
    converter = get_lifter('hg38', 'hg19', one_based=True)
    with open(input_file, 'r') as f:
        lines = f.readlines()[1:]

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

    return output_file

def run_snpeff(snpeff_jar_path, input_file, output_file): 
    snpeff_command = f"java -jar {snpeff_jar_path} -canon hg19 {input_file} > {output_file}"
    
    subprocess.run(snpeff_command, shell=True, check=True)
    
    return output_file

def format_annotations(input_file, output_file):
    with open(input_file, 'r') as anno_in, open(output_file, 'w') as anno_out:
        annod = OrderedDict.fromkeys(['CHROM', 'POS', 'REF', 'ALT', 'FAF', 'GENE', 'PCHANGE', 'EFFECT', 'EXON_Rank'])
        anno_out.write('\t'.join(annod.keys()) + '\n')

        for line in anno_in:
            if line.startswith('##') or line.startswith('#') or line.startswith('CHROM'):
                continue
            values = line.strip().split('\t')
            annod['CHROM'] = values[0]
            annod['POS'] = values[1]
            annod['REF'] = values[3]
            annod['ALT'] = values[4]
            if values[4] == '*': #snpeff won't annotate these so move on
                continue
            annod['FAF'] = values[5]

            sample_ann = values[7]
            try:
                anns = sample_ann.split('=')[1].split('|')
            except:
                print(sample_ann)
                continue
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

   # df.to_csv('/home/ec2-user/Azurify/results/domain.tsv',sep='\t', index=False)

    return(df)

def add_litvar(df):
    litvar_file_path = os.path.join(cwd, 'utils/litvar_pchanges.txt')
    litvar_df = pd.read_csv(litvar_file_path, sep='\t', header=None, names=['PCHANGE', 'PMID_COUNT'], dtype={'PCHANGE': str, 'PMID_COUNT': int})

    merged_df = pd.merge(df, litvar_df, on="PCHANGE", how='left')
   # df.to_csv('/home/ec2-user/Azurify/results/litvar.tsv',sep='\t', index=False)
    return(merged_df)

def add_kegg(df):
    kegg_file_path = os.path.join(cwd, 'utils/kegg_pathways.txt')
    kegg_df = pd.read_csv(kegg_file_path, sep='\t', header=None, names=['GENE'], dtype={'PCHANGE': str})

    kegg_df['IN_KEGG'] = 1
    plus_kegg = pd.merge(df,kegg_df[['GENE', 'IN_KEGG']],on='GENE', how='left')
    plus_kegg["KEGG"] = plus_kegg['IN_KEGG'].isnull()*1
    plus_kegg = plus_kegg.drop(['IN_KEGG'], axis=1)
   # plus_kegg.to_csv('/home/ec2-user/Azurify/results/kegg.tsv',sep='\t', index=False)
    return(plus_kegg)

def add_mvp(df):
    utilp = os.path.join(cwd, 'utils')
    df['KEY'] = df['CHROM'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']

    for filename in os.listdir(utilp):
        if filename.startswith('mvp'):
            files = glob.glob(os.path.join(utilp, "mvp*"))
            mvp = [pd.read_csv(f, sep="\t", low_memory=False) for f in files]
            k = pd.concat(mvp,ignore_index=True)
            df = df.merge(k, on='KEY', how='left')
            break
   # df.to_csv('/home/ec2-user/Azurify/results/mvp.tsv',sep='\t', index=False)
    return(df)


def merge_keys(df):
    utilp = os.path.join(cwd, 'utils')
    df['KEY'] = df['CHROM'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']

    for filename in os.listdir(utilp):
        if filename.endswith('_key.tsv.gz'):
            k = pd.read_csv((os.path.join(utilp, filename)), sep='\t',low_memory=False)
            df = df.merge(k, on='KEY', how='left')
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
    parser.add_argument("-i", "--input_file", metavar="FILE_PATH", required=False, type=str, help="Specify the input file.")
    parser.add_argument("-o", "--output_directory", metavar="OUTPUT_DIRECTORY", required=False,type=str, help="Specify where the output and intermediate files should be written.")
    parser.add_argument("-g", "--geneom_build", metavar="GENOME_BUILD", required=False, type=str, help="Specify the input genome, i.e hg19 or hg38.")
    parser.add_argument("-s", "--snpeff_jar_path", metavar="SNPEFF_JAR_PATH", required=True, type=str, help="Specify the path to the snpeff jar file.")
    parser.add_argument("-d", "--no_drug_targets", metavar="Drug Targets", required=False, type=str, help="Omit the use of drug targets in the model.")  
    # Parse the command-line arguments
    args = parser.parse_args()

    # Check for required arguments
    if not args.input_file or not args.output_directory:
        parser.error("Input file and Output Directory are required. Ex. Usage: python azurify.py -i /path/input.tsv -o /path/output/")
    return args

def main():
    print_ascii()

    #global progress bar
    progress_bar = tqdm(total=100, desc="Progress", unit="iter")

    #grab args and files
    args = parse_arguments()
    out_dir = args.output_directory
    uin = args.input_file
    udf = pd.read_csv(uin, sep='\t',low_memory=False)
    progress_bar.update(5)


    #set up directory and naming
    os.makedirs(out_dir, exist_ok=True)
    prefix = uin.split('.')[0]


    #split the df into columns used/not used by model
    in_cols=['CHROM','START','STOP','REF','ALT','VAF']
    adf, xdf = split_df(udf, in_cols)

    #convert to hg19 if required
    if args.geneom_build == 'hg38':
        o19 = os.path.join(out_dir, (prefix + ".hg19.txt"))
        uin = convert38(uin, o19)
        progress_bar.update(5)

    #run snpeff
    so = os.path.join(out_dir, (prefix + ".snpeffhg19.txt"))
    format_in = run_snpeff(args.snpeff_jar_path, uin, so)
    progress_bar.update(10)
     
    #format output
    format_out = os.path.join(out_dir, (prefix + ".snpeff.hg19.az_format.tsv"))
    az_in = format_annotations(format_in, format_out)

    
    ddf = pd.read_csv(az_in, sep='\t',low_memory=False)


    #split the df into columns used/not used by model
   # in_cols=['CHROM','POS','REF','ALT','FAF','GENE','PCHANGE','EFFECT', 'EXON_Rank']
   # adf, xdf = split_df(udf, in_cols)

    #add features
    ddf = add_domain(ddf)
    progress_bar.update(10)
    kdf = add_kegg(ddf)
    progress_bar.update(10)
    ldf = add_litvar(kdf)
    progress_bar.update(10)
    mvf =add_mvp(ldf)

    #add key feeatures
    mdf = merge_keys(mvf)
    progress_bar.update(20)

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
    mdf = mdf.round({'BP': 5,'PP': 5,'LBP': 5, 'LPP': 5, 'VP': 5})

    progress_bar.update(20)
    final = pd.concat([mdf, xdf], axis=1)

    final = final.drop(['KEY'], axis=1)

    #replace to standard nomenclature
    final = final.rename(columns={'FAF': 'VAF'})
    final['Pathogenicity'] = final['Pathogenicity'].astype(str).str.replace('Disease Associated', 'Pathogenic')
    final['Pathogenicity'] = final['Pathogenicity'].astype(str).str.replace('VOUS', 'VUS')
    final['Pathogenicity'] = final['Pathogenicity'].astype(str).str.replace('Probably DA', 'Likely Pathogenic')

    out_file = os.path.join(out_dir, (prefix + ".azurify.tsv"))
    final.to_csv(out_file, sep="\t", index=False)
    progress_bar.update(10)
    progress_bar.close()

    print("Classifications Complete, results written to: " + out_file)

if __name__ == '__main__':
    main()
