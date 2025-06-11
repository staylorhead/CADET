#!/usr/bin/env python

#########################################################
# import dependencies 
#########################################################
import functools
import os
import subprocess
import sys
import traceback
from numpy import linalg
from functools import reduce
from io import StringIO
from scipy.stats.distributions import chi2
from scipy.stats import norm
import pandas as pd
import numpy as np
import getopt
import pysam
import glob
import resource
import multiprocessing
import shutil
from time import time

#########################################################
# load child functions
#########################################################

# return human readable elapsed time string
def format_elapsed_time(time_secs):
    val = abs(int(time_secs))
    day = val // (3600*24)
    hour = val % (3600*24) // 3600
    mins = val % 3600 // 60
    secs = val % 60
    res = '%02d:%02d:%02d:%02d' % (day, hour, mins, secs)
    if int(time_secs) < 0:
        res = "-%s" % res
    return res


# Decrease memory by downcasting 'CHROM' column to integer, integer and float columns to minimum size that will not lose info
def optimize_cols(df: pd.DataFrame):
    if 'CHROM' in df.columns:
        df['CHROM'] = df['CHROM'].astype(str).astype(int)
    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')
    floats = df.select_dtypes(include=['float64']).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast='float')
    return df


# determine indices of file cols to read in, dtype of each col
def get_cols_dtype(cols):
    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'A1': object,
        'A2': object,
        'SNP': object,
        'ES': np.float64,
        'Z': np.float64,
        'P': np.float64,
        'Beta': np.float64,
        'N': np.int64,
        'GeneEnd': np.int64,
        'GeneName': object,
        'GeneStart': np.int64,
        'snpID': object,
        'TargetID': object,
        'bp': object}
    out_dtype_dict = {x: dtype_dict[x] for x in cols}
    return out_dtype_dict

def tabix(in_file, format='bgzip', tmp_file=None):
    """ Tabix a txt or txt.gz file with option to keep the original file or not """
    # if the file is not compressed, firt compress and save to a temporary file
    if format == 'txt':
        try:
            pysam.tabix_compress(in_file, tmp_file, force=True)
        except:
            print('Compress failed for ' + in_file)
            return None
    else:
        tmp_file = in_file
    try:
        out_file = pysam.tabix_index(tmp_file, force=True,
                                     seq_col=0, start_col=1, end_col=1,
                                     line_skip=1)
    except:
        print('Tabix failed for ' + in_file)
        return None
    return out_file    

# check if a path exists. If not, create the path
def check_path(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def create_file_title(out_cols, out_dir, out_file):
    out_path = os.path.join(out_dir, out_file)
    pd.DataFrame(columns=out_cols).to_csv(
        out_path,
        sep='\t',
        index=None,
        header=True,
        mode='w')

def get_snpIDs(df: pd.DataFrame, flip=False):
    chroms = df['CHROM'].astype('str').values
    pos = df['POS'].astype('str').values
    A1 = df['A1'].values
    A2 = df['A2'].values
    if flip:
        return ['_'.join(i) for i in zip(chroms, pos, A1, A2)]
    else:
        return ['_'.join(i) for i in zip(chroms, pos, A2, A1)]

def read_format_ref_bim(ref_dir, ref_file):
    target_bim = os.path.join(ref_dir, ref_file)
    col_names = ["CHROM", "SNP", "bp", "POS", "A1", "A2"]
    dtypes = get_cols_dtype(col_names)
    ref_chunks = pd.read_csv(target_bim, sep='\t',
                             low_memory=False,
                             header=None,
                             names=col_names,
                             iterator=True,
                             chunksize=1000,
                             dtype=dtypes)
    target_ref = pd.concat([chunk for chunk in ref_chunks]).reset_index(drop=True)
    if len(target_ref) == 0:
        return None
    # format snp IDs in the reference bim file
    target_ref['snpID'] = get_snpIDs(target_ref, flip=False)
    target_ref[['CHROM', 'snpID', 'bp', 'POS', 'A1', 'A2']].to_csv(
        target_bim,
        sep='\t',
        index=None,
        header=None,
        mode='w')
    return target_ref

def read_in_clumped(clumped_file, chrom):
    with open(clumped_file) as file_in:
        lines = []
        snps = []
        for line in file_in:
            lines.append(line)
            res = line.split()
            if len(res) != 0:
                if res[0] == chrom:
                    snps.append(res[2])
    return snps

def handle_flip_snps(df_ref, df_1, statistics):
    """Handle flipped snp IDs"""
    # filter out non-matching snpID rows in df_ref and df_1
    df_1 = df_1[np.any(df_1[['snpID', 'snpIDflip']].isin(df_ref.snpID.values), axis=1)].reset_index(drop=True)
    # if snpID is not in df_ref.snpID, assumed flipped; if flipped, flip Zscore sign
    df_1['flip'] = np.where(df_1.snpID.isin(df_ref.snpID.values), 1, -1)
    if not np.all(df_1['flip'] == 1):
            idx = (df_1['flip'] == -1)
            df_1.loc[idx, ['snpID']] = df_1.loc[idx, ['snpIDflip']].values
            df_1.loc[idx, ['A1', 'A2']] = df_1.loc[idx, ['A2', 'A1']].values
            df_1[statistics] = df_1['flip'] * df_1[statistics]
    return df_1


def match_snp_ID_double(df_ref, df_1):
    df_ref['snpID'] = get_snpIDs(df_ref, flip=False)
    df_1['snpID'] = get_snpIDs(df_1, flip=False)
    df_1['snpIDflip'] = get_snpIDs(df_1, flip=True)
    # find overlapped snpIDs
    overlap_list = np.intersect1d(df_ref['snpID'], df_1[['snpID', 'snpIDflip']])
    if overlap_list.size:
        df_ref = df_ref[df_ref.snpID.isin(overlap_list)]
        df_1 = handle_flip_snps(df_ref, df_1, 'Beta')
    else:
        return None, None, None
    return df_ref, df_1, overlap_list

def output_data(out_df, out_cols, out_path, write_mode, out_header=None):
    out_df[out_cols].to_csv(
        out_path,
        sep='\t',
        index=None,
        header=out_header,
        mode=write_mode)

def call_PLINK_clump(bim_path, r2, pvalue_path, work_dir, window=1000, p=1, snp_field="snpID", p_field="P"):
    cmd = ["plink --bfile " + bim_path+" --clump-p1 " + str(p) + " --clump-r2 " + str(r2) +
           " --clump-kb " + str(window)+" --clump " + pvalue_path +
           " --clump-snp-field " + snp_field + " --clump-field " + p_field +
           " --keep-allele-order --out "+bim_path]
    # perform LD-clumping
    try:
        proc = subprocess.check_call(cmd,
                                     stdout=subprocess.PIPE,
                                     cwd=work_dir,
                                     shell=True)
        print('Done LD Clumping.')
    except subprocess.CalledProcessError:
        print('LD Clumping Failed. \n')
        return None

def filter_df_rows(df, filter_by, filter_cols):
    df = df[np.any(df[filter_cols].isin(filter_by), axis=1)].reset_index(drop=True)
    return df

def save_results(model, out_df, out_dir):
    out_df[['CHROM', 'POS', 'A1', 'A2', 'TargetID', model]].to_csv(
        os.path.join(out_dir, model + '.txt'),
        sep='\t',
        index=None,
        header=None,
        mode='a')
    print('Finish ' + model + '.')

def save_results_PT(p, out_df, out_dir):
    out_df[['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'Beta']].to_csv(
        os.path.join(out_dir, 'P' + str(p) + '.txt'),
        sep='\t',
        index=None,
        header=None,
        mode='a')
    print('Finish P+T with p-value < ' + str(p) + '.')


def pos_def_matrix(mat):
    """ convert the input matrix to the cloest positive definite matrix"""
    # Make sure the ld is positive definite matrix
    _, s, v = linalg.svd(mat)
    h = np.dot(v.T, np.dot(np.diag(s), v))
    mat_pos_def = (mat+h)/2
    return mat_pos_def

def format_save_results(work_dir, out_dir, model, sst_df):
    raw_out_dir = os.path.join(work_dir, model + '.txt')
    if model == 'SDPR':
        if not os.path.exists(raw_out_dir):
            print('SDPR failed. Please check the input arguments for SDPR')
            return None
        else:
            raw_names = ["snpID", "A1", model]
    elif model == 'lassosum':
        if not os.path.exists(raw_out_dir):
            print('lassosum failed. Please check the input arguments for lassosum')
            return None
        else:
            raw_names = ['CHROM', 'POS', 'A1', 'A2', model]
    raw_chunks = pd.read_csv(raw_out_dir, sep='\t',
                             low_memory=False,
                             header=0,
                             names=raw_names,
                             iterator=True,
                             chunksize=1000)
    raw_weights = pd.concat([chunk for chunk in raw_chunks]).reset_index(drop=True)
    if model == 'SDPR':
        out_df = raw_weights.merge(sst_df, left_on=['snpID', 'A1'],
                                   right_on=['snpID', 'A1'], how="inner")
    elif model == 'lassosum':
        out_df = raw_weights.merge(sst_df, on=['CHROM', 'POS', 'A1', 'A2'],
                                   how="inner")
    # remove the raw outputs
    os.remove(raw_out_dir)
    save_results(model, out_df, out_dir)


def lassosum_cmd(chrom, bim_dir, sst_dir, out_dir, lassosum_path, lassosum_LD_block):
    cmd = ['Rscript ' + lassosum_path +
           ' --bim_file=' + bim_dir +
           ' --sst_file=' + sst_dir +
           ' --out_path=' + out_dir +
           ' --chr=' + str(chrom) +
           ' --LDblocks=' + lassosum_LD_block]
    return cmd

def read_sst_by_chunks(tabix_out, sst, target):
    if not tabix_out:
        return None
    if sst == 'eQTL weights':
        col_names = ["CHROM", "POS", "A1", "A2", "TargetID", "ES"]
    elif sst == 'GWAS':
        col_names = ["CHROM", "POS", "A1", "A2", "Z"]
    elif sst == 'eQTL sst':
        col_names = ["CHROM", "POS", "A1", "A2", "Beta","P", "TargetID", "N"]
    dtypes = get_cols_dtype(col_names)
    chunks = pd.read_csv(StringIO(tabix_out.decode('utf-8')), sep='\t',
                         low_memory=False,
                         header=None,
                         names=col_names,
                         iterator=True,
                         chunksize=1000,
                         dtype=dtypes)
    if sst == 'GWAS':
        target_df = pd.concat([chunk for chunk in chunks]).drop_duplicates(["CHROM", "POS", "A1", "A2"]).reset_index(drop=True)
    else:
        target_df = pd.concat([chunk[chunk.TargetID == target] for chunk in chunks]).drop_duplicates(["CHROM", "POS", "A1", "A2"]).reset_index(drop=True)
    if target_df.empty:
        return None
    target_df = optimize_cols(target_df)
    return target_df

def read_sst(sst_file, sst, target, chrom, start_pos, end_pos):
    # call tabix to extract estimated eQTL effect sizes for target gene
    tabix_out = call_tabix(sst_file, chrom, start_pos, end_pos)
    # read in estimated eQTL effect sizes
    target_df = read_sst_by_chunks(tabix_out=tabix_out,
                                   sst=sst,
                                   target=target)
    return target_df


def call_tabix(path, chrom, start, end):
    chrom = str(chrom)
    proc = subprocess.Popen(
        ["tabix "+path+" "+chrom+":"+start+"-"+end],
        shell=True,
        stdout=subprocess.PIPE)
    proc_out = bytearray()
    # process while subprocesses running
    while proc.poll() is None:
        line = proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += line
    # get any remaining lines
    for line in proc.stdout:
        proc_out += line
    return proc_out

def match_snp_ID_impute(df_ref, df_1):
    df_ref['snpID'] = get_snpIDs(df_ref, flip=False)
    df_1['snpID'] = get_snpIDs(df_1, flip=False)
    df_1['snpIDflip'] = get_snpIDs(df_1, flip=True)
    # find overlapped snpIDs
    overlap_list = np.intersect1d(df_ref['snpID'], df_1[['snpID', 'snpIDflip']])
    if overlap_list.size:
        df_ref = df_ref[df_ref.snpID.isin(overlap_list)]
        df_1 = handle_flip_snps(df_ref, df_1, 'ES')
    else:
        return None, None, None
    return df_ref, df_1, overlap_list

def call_PLINK_extract(bim_path, out_path, target, chrom, start_pos, end_pos):
    # save the range of the gene
    range = os.path.join(out_path, 'range.txt')
    with open(range, 'w') as ff:
        ff.write('%s\t%s\t%s\t%s\n' % (chrom, start_pos, end_pos, target))
    # extract the genotype data for this range
    out_geno = os.path.join(out_path, target)
    cmd = ["plink --bfile "+bim_path+" --keep-allele-order --extract range " + range + " --make-bed --out " + out_geno]
    try:
        proc = subprocess.check_call(cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True)
    except subprocess.CalledProcessError:
        print('There is no genotype reference data.')
        return None
    return True

def prepare(target, target_anno, chrom, window,
            geno_dir, out_dir, sst_dir, clump_r2):
    ################# PLINK Binary Files #####################
    print('...Preparing Input Files...')
    start = str(max(int(target_anno.GeneStart) - window, 0))
    end = str(int(target_anno.GeneEnd) + window)
    # set output path of the target gene
    target_dir = os.path.join(out_dir, target)
    check_path(target_dir)
    # generate command to call PLINK to extract the binary file for the target gene
    extract_proc = call_PLINK_extract(bim_path=geno_dir,
                                          out_path=target_dir,
                                          target=target,
                                          chrom=chrom,
                                          start_pos=start,
                                          end_pos=end)
    if not extract_proc:
        print('Remove temporary files. \n')
        shutil.rmtree(target_dir)
        return None, None, None
    ################# Read in eQTL summary statistics #####################
    target_sst = read_sst(sst_file=sst_dir,
                              sst='eQTL sst',
                              target=target,
                              chrom=chrom,
                              start_pos=start,
                              end_pos=end)
    if target_sst is None:
        print('There is no estimated eQTL summary statistics')
        print('Remove temporary files. \n')
        shutil.rmtree(target_dir)
        return None, None, None
    ################# Read in SNPs in LD reference #####################
    # read in snps in LD reference panel
    target_ref = read_format_ref_bim(ref_dir=target_dir,
                                         ref_file=target + '.bim')
    if target_ref is None:
        print('There is no reference bim file.')
        return None
    ########### Match SNPs in eQTL reference and LD reference #########
    # print('Check overlapped SNPs between eQTL reference and LD reference...')
    target_ref, target_sst, snp_overlap = match_snp_ID_double(df_ref=target_ref,
                                                                  df_1=target_sst)
    if not snp_overlap.size:
        print('No overlapping test eQTLs')
        print('Remove temporary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None, None, None
    ################ Calculate median sample size ####################
    # print('*Calculate median sample size of eQTLs...*')
    #median_N = np.nanmedian(target_sst['N'])
    ################# LD clumping #############################
    # print('*Perform LD clumping...*')
    # generate summary statistics of p-value to perform LD-clumping
    # chi.sf returns 0 for large Z scores
    # we divided it by 5 here to shrink Zscore to prevent p-values = 0
    #target_sst['P'] = chi2.sf(np.power(target_sst['Z']/5, 2), 1)
    output_data(out_df=target_sst,
                    out_cols=['snpID', 'A1', 'A2', 'P'],
                    out_path=os.path.join(target_dir, target + '.pvalue'),
                    write_mode='w',
                    out_header=True)
    # use PLINK to perform LD-clumping 
    if clump_r2 < 1:
        call_PLINK_clump(target, clump_r2, target + '.pvalue', target_dir, window/1000)
        # read in remaining eQTLs after LD-clumping
        target_clumped = os.path.join(target_dir, target + '.clumped')
        clumped_snp = read_in_clumped(clumped_file=target_clumped,
                                      chrom=chrom)
        # filter summary statistics by clumping results
        target_sst = filter_df_rows(df=target_sst,
                                    filter_by=clumped_snp,
                                    filter_cols=['snpID'])
        print('Number of SNPs remaining after LD clump in ref:' + ' ' + str(len(target_sst.index)))
    else:    
        print('No LD clumping in ref performed.')
        print('Number of SNPs:' + ' ' + str(len(target_sst.index)))    
    ####### Prepare Inputs for Imputation Models ########
    # print('*Start prepare input summary statistics...*')
    # Prepare Zscore input for SDPR
    # print('Done generating Zscore.')
    #ots.output_data(out_df=target_sst,
    #                out_cols=['snpID', 'A1', 'A2', 'Z'],
    #                out_path=os.path.join(target_dir, target+'_Zscore.txt'),
    #                write_mode='w',
    #                out_header=['SNP', "A1", "A2", "Z"])
    # Prepare the standardized beta input for lassosum and P+T
    #target_sst['Beta'] = target_sst['Z']/np.sqrt(median_N)
    output_data(out_df=target_sst,
                    out_cols=['snpID', 'A1', 'A2', 'Beta'],
                    out_path=os.path.join(target_dir, target+'_beta.txt'),
                    write_mode='w',
                    out_header=['SNP', "A1", "A2", "Beta"])
    # print('Done generating standardized beta.')
    print('Done prepare inputs.')
    return target_dir, target_sst

#########################################################
# parse argument
#########################################################

def parse_param():
    long_opts_list = ['anno_file=', 'geno_dir=','out_dir=',
                      'r2=',
                      'window=', 'thread=', 'models=',
                      'pt=',
                      'seed=',
                      'sst_file=', 'lassosum_LD_block=', 'script_dir=',
                      'help']
    param_dict = {'anno_file': None, 
                  'geno_dir': None, 'out_dir': None, 
                  'r2': 0.99, 'window': 1000000, 'thread': 1, 'models': None,
                  'pt': [0.001, 0.05],
                  'seed': None, 'sst_file': None, 'lassosum_LD_block': None,'script_dir': None}
    print('\n')
    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)
        except:
            print('Option not recognized.')
            print('Use --help for usage information.\n')
            sys.exit(2)
        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--anno_file":
                param_dict['anno_file'] = arg
            elif opt == "--geno_dir":
                param_dict['geno_dir'] = arg
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--r2":
                param_dict['r2'] = float(arg)
            elif opt == "--window":
                param_dict['window'] = int(arg)
            elif opt == "--models":
                param_dict['models'] = arg.split(',')
            elif opt == "--pt":
                param_dict['pt'] = arg.split(',')
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
            elif opt == "--seed":
                param_dict['seed'] = int(arg)
            elif opt == "--sst_file":
                param_dict['sst_file'] = arg
            elif opt == "--lassosum_LD_block":
                param_dict['lassosum_LD_block'] = arg
            elif opt == "--script_dir":
                param_dict['script_dir'] = arg    
    else:
        print(__doc__)
        sys.exit(0)
    if param_dict['script_dir'] is None:
        print('* Please specify the directory to CADET directory --script_dir\n')
        sys.exit(2)
    elif param_dict['anno_file'] is None:
        print('* Please specify the directory to the gene annotation file using --anno_dir\n')
        sys.exit(2)
    elif param_dict['geno_dir'] is None:
        print('* Please specify the directory to the binary file of LD reference panel --geno_dir\n')
        sys.exit(2)
    elif param_dict['out_dir'] is None:
        print('* Please specify the output directory\n')
        sys.exit(2)
    elif param_dict['models'] is None:
        print('* Please specify the imputation models --models\n')
        sys.exit(2)
    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))
    print('\n')
    return param_dict

param_dict = parse_param()

# Create directory for output files
out_dir = param_dict['out_dir']
check_path(out_dir)

# Create output files and load required tools
out_cols = ['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES', 'N']
for model in param_dict['models']:
    if model == 'PT':
        for p in param_dict['pt']:
            create_file_title(out_cols, out_dir, 'P' + str(p) + '.txt')
    elif model == 'lassosum':
        create_file_title(out_cols, out_dir, model + '.txt')
        lassosum_path = os.path.join(param_dict['script_dir'], 'PRSmodels/lassosum.R')
    else:
        print('Please specify models among: PT,lassosum')

# Create directory for temporary files
tmp_dir = os.path.join(out_dir, 'tmp')
check_path(tmp_dir)
print('Tabixing eQTL summary statistics...')
tabix_sst = tabix(in_file=param_dict['sst_file'], format='txt',
                      tmp_file=os.path.join(tmp_dir, 'sst.txt.gz'))

if not tabix_sst:
    sys.exit(2)

col_names = ['CHROM', 'GeneEnd', 'GeneStart', 'TargetID']
dtypes = get_cols_dtype(col_names)
anno = pd.read_csv(
        param_dict['anno_file'],
        sep='\t',
        header=0,
        usecols=col_names,
        dtype=dtypes)
target_vec = anno['TargetID']
n_targets = target_vec.size

#####  threaded defn ###
def thread_process(num):
    target = target_vec[num]
    print('num=' + str(num) + ' TargetID=' + target)
    target_anno = anno.iloc[[num]]
    param_dict['chrom'] = str(target_anno['CHROM'].iloc[0]) # set chrom to that of gene under study
    # prepare data
    target_dir, target_sst = prepare(target=target,
                                                    target_anno=target_anno,
                                                    chrom=param_dict['chrom'],
                                                    window=param_dict['window'],
                                                    geno_dir=param_dict['geno_dir'],
                                                    out_dir=out_dir,
                                                    sst_dir=tabix_sst,
                                                    clump_r2=param_dict['r2'])
    if (target_dir is None):
        return None
    ################# Start Imputation Models #############################
    print('...Training eQTL weights...')
    # Clumping and thresholding P+T #
    if 'PT' in param_dict['models']:
        # print('*Start P+T...*')
        #target_sst['P'] = chi2.sf(np.power(target_sst['Z'], 2), 1)
        for p in param_dict['pt']:
            PT_out = target_sst[target_sst.P < float(p)]
            save_results_PT(p=p, out_df=PT_out, out_dir=out_dir)
    # lassosum #
    if 'lassosum' in param_dict['models']:
        # print("*Start lassosum...*")
        try:
            lassosum_arg = lassosum_cmd(lassosum_LD_block=param_dict['lassosum_LD_block'], 
                                        chrom=param_dict['chrom'], 
                                        bim_dir=target, 
                                        sst_dir=target + '_beta.txt', 
                                        out_dir='lassosum.txt',
                                        lassosum_path=lassosum_path)
            proc = subprocess.check_call(lassosum_arg,
                                         stdout=subprocess.PIPE,
                                         cwd=target_dir,
                                         shell=True)
            # save lassosum results
            format_save_results(work_dir=target_dir,
                                    out_dir=out_dir,
                                    model='lassosum',
                                    sst_df=target_sst)
        except subprocess.CalledProcessError:
                print('lassosum failed for TargetID: ' + target)
    ############################ Clean temporary files #########################
    shutil.rmtree(target_dir)
    print('Done. \n')

print('Starting model training')
# time calculation
start_time = time()

############################################################
if __name__ == '__main__':
    print('Starting training eQTL weights for ' + str(n_targets) + ' target genes.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Removing temporary files.')
    shutil.rmtree(os.path.join(out_dir, 'tmp'))

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = format_elapsed_time(elapsed_sec)
print('Total computation time (DD:HH:MM:SS): ' + elapsed_time)

# peak memory usage
print('Peak memory usage in kilobytes:')
mem=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(mem)


