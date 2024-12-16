import pandas as pd
import numpy as np
from statsmodels.multivariate.manova import MANOVA
import statsmodels.api as sm
# import dask.dataframe as dd
import multiprocessing
from tqdm import tqdm
import argparse
import time
import os
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
from scipy.stats import ks_2samp

def group_func(name, temp_df, method,coverage_cutoff=20,balance=False):
    try:
        if temp_df.shape[0] >= coverage_cutoff * 2:
            label_list = temp_df["Label"].values
            temp_df.drop(["Label",'Position'], axis=1, inplace=True)

            sample = temp_df[label_list == 'Sample']
            control = temp_df[label_list == 'Control']
            sample_cov = sample.shape[0]
            control_cov = control.shape[0]
            if balance == True:
                if len(sample) > len(control):
                    sample = sample.sample(n=len(control), random_state=42)
                else:
                    control = control.sample(n=len(sample), random_state=42)
                temp_df = pd.concat([sample, control], axis=0)
                ones_array = np.ones(len(control))

                # 创建一个长度为 n 的全为 0 的 NumPy 数组
                zeros_array = np.zeros(len(control))

                # 将两个数组拼接在一起
                label_list = np.concatenate((ones_array, zeros_array))
            else:
                label_list = np.where(label_list == 'Sample', 1, np.where(label_list == 'Control', 0, label_list))
                label_list = label_list.astype(int)
            # new_df['Label'] = label_list
            # new_df.columns = ['PC1', 'PC2', "Group"]

            if sample.shape[0] >= coverage_cutoff and control.shape[0] >= coverage_cutoff:
                mean_differ = sample['Mean'].mean() - control['Mean'].mean()
                median_differ = sample['Median'].mean() - control['Median'].mean()
                dwell_differ = sample['Dwell time'].mean() - control['Dwell time'].mean()
                std_differ = sample['STD'].mean() - control['STD'].mean()
                if method == "lr":
                    new_df = sm.add_constant(temp_df)  # 增加模型的常數，使更為符合回歸模型
                    model = sm.OLS(label_list, new_df)  # OLS回歸
                    results = model.fit()
                    pvalue = results.pvalues.iloc[1]
                    del new_df
                elif method=='manova':
                    temp_df['Group'] = label_list
                    temp_df.columns=['Mean', 'STD', 'Median', 'dwell_time', 'Group']
                    manova = MANOVA.from_formula('Mean + STD + Median + dwell_time ~ Group', data=temp_df)
                    # 执行多元方差分析
                    results = manova.mv_test()
                    pvalue = results.summary().tables[3].iloc[0, 5]

                elif method == 'ks':
                    statistic, pvalue = ks_2samp(sample['Mean'].values, control['Mean'].values,method='asymp')
                else:
                    raise ValueError('Method must be either "lr" or "manova" or "ks"')
                line = name.split(":")
                line.extend([sample_cov,control_cov, pvalue, mean_differ,median_differ,dwell_differ,std_differ])
                line = [str(element) for element in line]
                del control, sample
                return line
                # pbar.update(1)
        return None
    except Exception as e:
        print(1)

def update(*a):
    if a[0] is not None:
        output_file.write(','.join(a[0]) + '\n')
    pbar.update(1)


def list_files_in_directory(directory_path):
    result_list=[]
    files_and_dirs = os.listdir(directory_path)
    for item in files_and_dirs:
        if '.csv' in item and '_' in item:
            result_list.append(item)
    return result_list

def detect_chunks(input_path, control_path):
    input_files = set(list_files_in_directory(input_path))
    control_files = set(list_files_in_directory(control_path))
    chunks_name = input_files.intersection(control_files)
    if len(chunks_name) == 0:
        raise Exception("No same file name found in both folders")
    sorted_list_desc = sorted(chunks_name, reverse=False)
    return sorted_list_desc
def merge_df(input_df, control_df):
    df_Sample = pd.read_csv(input_df, header=None)
    df_Sample['Label'] = 'Sample'

    df_Control = pd.read_csv(control_df, header=None)
    df_Control['Label'] = 'Control'

    df = pd.concat([df_Control, df_Sample], axis=0)
    df.columns = ['Read ID', 'Chrom', 'Position','Strand','Mean', 'STD', 'Median', 'Dwell time', 'Label']
    df['Position'] = df['Chrom']+':'+df['Position'].astype(str)+":"+df['Strand']
    df.drop(['Read ID', 'Chrom', 'Strand'], axis=1, inplace=True)
    del df_Sample, df_Control
    return df
def main(input_path, control_path, output_path, cpu, balance,method,coverage_cutoff):
    global pbar, output_file
    chunks_name = detect_chunks(input_path, control_path)
    timestamp = int(time.time())
    tmp_output_name='nanosundial_'+str(timestamp)+'.csv'
    output_file = open(tmp_output_name, 'w+')
    with tqdm(total=len(chunks_name),unit='chunk') as pbar_chunk:
        for chunk_name in chunks_name:
            input_file = input_path +'/'+ chunk_name
            control_file = control_path+ '/' + chunk_name
            df = merge_df(input_file, control_file)
            df_groups = df.groupby('Position',sort=True)
            pbar = tqdm(total=df_groups.ngroups, position=0, leave=True)
            with multiprocessing.Pool(processes=cpu) as pool:
                for key,temp_df in df_groups:
                    pool.apply_async(group_func, args=(key, temp_df, method, coverage_cutoff,balance), callback=update)
                pool.close()
                pool.join()
            pbar.close()
            pbar_chunk.update(1)
            del df,df_groups
    output_file.close()
    print('Start to run FDR correction')

    df = pd.read_csv(tmp_output_name)
    df.columns = ['chrom', 'position', 'strand', 'sample_coverage', 'control_coverage', 'pvalue', 'mean_differ','median_differ',
                  'dwell_differ','std_differ']
    rejected, df['adj_p'], _, _ = multipletests(df['pvalue'].values, method='fdr_bh')
    if output_path.split('.')[-1] != 'csv':
        output_file=output_path+'.csv'
    df.sort_values(by=['chrom', 'position'], inplace=True)
    df.to_csv(output_path,index=False)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",required=True,help="Input file path")

    parser.add_argument('-c', "--control",required=True,help="Input file path")
    parser.add_argument('-o', "--output", default='nanosundial_result.csv',
                        help="output file")
    parser.add_argument('-t', "--cpu", default=4, type=int, help="Process numbers")
    parser.add_argument("--balance", action='store_true', help="Turn on the balance mode")
    parser.add_argument("--method", choices=['manova', 'lr', 'ks'], default='manova',
                        help="statistical algorithms provided(default:manova)")
    parser.add_argument("--minimum_coverage", default=20, type=int,
                        help="coverage cutoff of this comparison")
    args = parser.parse_args()

    # args.subsample = 1
    start_time = time.time()
    args.balance = True
    main(args.input, args.control, args.output, args.cpu, args.balance, args.method,args.minimum_coverage)

    end_time = time.time()
    run_time = end_time - start_time
    # Convert to minutes and hours
    minutes = int(run_time // 60)  # Integer division to get the whole minutes
    seconds = int(run_time % 60)  # Remainder to get the remaining seconds

    hours = int(minutes // 60)  # Integer division to get the whole hours
    minutes %= 60  # Modulo operation to get the remaining minutes

    # Print the running time
    print("Program running time:", hours, "hours", minutes, "minutes", seconds, "seconds")