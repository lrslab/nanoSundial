import os
import time
from collections import OrderedDict
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from scipy import stats
import copy

kmer_model_size ={
    'r9+RNA': 5,
    'r9+DNA': 6,
    'r10+DNA': 9,
    'rna004+RNA': 9
}

base_shift_dict ={
    'r9RNA+': -1,
    'r9RNA-': -3,
    'r9DNA+': -2,
    'r9DNA-': -3,
    'r10DNA+': -6,
    'r10DNA-': -2,
    'rna004RNA+': -3,
    'rna004RNA-': -5,
}


def run_cmd(cmd):
    try:
        # 执行命令并捕获输出
        print(cmd)
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
        # 处理命令输出
        # print(output)
    except subprocess.CalledProcessError as e:
        # 处理命令执行错误
        print("Command execution failed:", cmd)
        print("Reason:", e.output)
        raise RuntimeError('There are some errors in the cmd as below, please check your env ')
def read_fasta_to_dic(filename):
    """
    function used to parser small fasta
    still effective for genome level file
    """
    fa_dic = OrderedDict()

    with open(filename, "r") as f:
        for n, line in enumerate(f.readlines()):
            if line.startswith(">"):
                if n > 0:
                    fa_dic[short_name] = "".join(seq_l)  # store previous one

                full_name = line.strip().replace(">", "")
                short_name = full_name.split(" ")[0]
                seq_l = []
            else:  # collect the seq lines
                if len(line) > 8:  # min for fasta file is usually larger than 8
                    seq_line1 = line.strip()
                    seq_l.append(seq_line1)

        fa_dic[short_name] = "".join(seq_l)  # store the last one
    return fa_dic


def identify_file_path(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError("File do not exist! Please check your path : " + file_path)


def generate_bam_file(fastq_file, reference, cpu,output,subsample_ratio=1):


    bam_file = output+'/' + '.'.join(fastq_file.split('.')[:-1]) + '_aligned.bam'

    if not os.path.exists(bam_file):
        cmds = 'minimap2 -ax map-ont -t ' + cpu + ' --MD --secondary=no ' + reference + ' ' + fastq_file + ' | samtools view -hbS -F ' + str(
            3332) + '  - | samtools sort -@ ' + cpu + ' -o ' + bam_file
        print('Start to alignment ...')
        run_cmd(cmds)
        print('bam file is saved in ' + bam_file)
    else:
        print(bam_file + ' existed. Will skip the minimap2 ... ')

    if subsample_ratio < 1:
        new_bam = output+'/' +'.'.join(bam_file.split('.')[:-1]) + '_sub.bam'
        if not os.path.exists(new_bam):
            cmds = "samtools view -hbS -s " +str(subsample_ratio) +' ' + bam_file +' > ' + new_bam
            run_cmd(cmds)
        else:
            print(bam_file + ' existed. Will skip the subsample ... ')
        bam_file = new_bam

    if not os.path.exists(bam_file+'.bai'):
        cmds = 'samtools index ' + bam_file
        run_cmd(cmds)

    new_fastq_file = '.'.join(bam_file.split('.')[:-1]) + '.fastq'
    if not os.path.exists(new_fastq_file):
        cmds = 'samtools bam2fq ' + bam_file + ' > '+ new_fastq_file
        run_cmd(cmds)
    return new_fastq_file, bam_file

def generate_paf_file_eventalign(fastq_file, blow5_file,bam_file,fasta_file,pore,rna,cpu,output):
    paf_file =  '.'.join(fastq_file.split('.')[:-1]) + '_ev.paf'
    if not os.path.exists(paf_file):
        if not os.path.exists(blow5_file+'.idx'):
            cmds = 'slow5tools index ' + blow5_file
            run_cmd(cmds)

        cmds = 'f5c index --slow5 ' +blow5_file+' '+ fastq_file
        run_cmd(cmds)

        cmds = 'f5c eventalign -r '+ fastq_file +" -g "+fasta_file+ ' --slow5 ' + blow5_file + ' --pore '+ pore+' -b ' + bam_file +' -c --min-mapq 0' + ' -t ' + str(cpu)
        if rna:
            cmds =cmds +' --rna'
        cmds =cmds +' > '+ paf_file
        # print(cmds)
        print('Start to eventalign ...')
        run_cmd(cmds)
        print('Generated paf file : ' + paf_file)
    else:
        print(paf_file + ' existed. Will skip the f5c eventalign ... ')
    return paf_file

def generate_paf_file_resquiggle(fastq_file, blow5_file,pore,rna,cpu):
    paf_file =  '.'.join(fastq_file.split('.')[:-1]) + '_re.paf'
    if not os.path.exists(paf_file):
        if not os.path.exists(blow5_file+'.idx'):
            cmds = 'slow5tools index ' + blow5_file
            run_cmd(cmds)

        cmds = 'f5c index --slow5 ' +blow5_file+' '+ fastq_file
        run_cmd(cmds)

        cmds = 'f5c resquiggle -c ' + fastq_file + ' ' + blow5_file + ' --pore ' + pore + ' -o ' + paf_file +' -t '+ str(cpu)
        if rna:
            cmds =cmds +' --rna'
        print('Start to f5c resquiggle ...')
        run_cmd(cmds)
        print('Generated paf file : ' + paf_file)
    else:
        print(paf_file + ' existed. Will skip the f5c resquiggle ... ')
    return paf_file

def prepare_move_table_file(bam_file, reference, cpu,sig_move_offset,kmer_length):
    paf_file = '.'.join(bam_file.split('.')[:-1]) + '.paf'
    fastq_file = '.'.join(bam_file.split('.')[:-1]) + '.fastq'
    cmds = 'samtools index ' + bam_file
    run_cmd(cmds)
    print('Start to generate paf file  ...')
    cmds = 'squigualiser reform --sig_move_offset '+sig_move_offset+' --kmer_length '+kmer_length+' -c --bam ' + bam_file +' -o ' + paf_file
    print(cmds)
    run_cmd(cmds)
    print("Start alignment ...")
    cmds = 'samtools bam2fq '+bam_file+' >' + fastq_file
    run_cmd(cmds)
    aligned_fastq,aligned_bam = generate_bam_file(fastq_file, reference, cpu)
    return aligned_bam,paf_file


def run_samtools(fastq_file, location, reference, result_path, group, cpu):
    _,bam_file = generate_bam_file(fastq_file, reference, cpu)
    cmds = 'samtools mpileup ' + bam_file + ' -r ' + location + ' --no-output-ins --no-output-del -B -Q 0 -f ' + reference + ' -o ' + result_path + 'temp.txt'
    run_cmd(cmds)
    temp_file = pd.read_csv(result_path + 'temp.txt', sep='\t', header=None)
    temp_file['Group'] = group
    run_cmd('rm ' + result_path + 'temp.txt')
    return temp_file


def build_out_path(results_path):
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    else:
        print("Output file existed! It will be overwrite or add content after 5 secs ...")
        time.sleep(5)
        print("Continue ...")

def reverse_fasta(ref):
    base_dict = {'A': 'T', "T": "A", "C": "G", "G": "C"}
    ref = np.array(list(ref.upper()))
    for index, element in enumerate(ref):
        ref[index] = base_dict[element]
    reverse_fasta = list(reversed(''.join(ref)))
    return ''.join(reverse_fasta)

def extract_kmer_feature(df_input, kmer, position):
    df = copy.deepcopy(df_input)
    def is_odd(number):
        if number % 2 == 1:
            return True
        else:
            return False
    if not is_odd(kmer) or kmer <= 0:
        raise Exception("The kmer should be an odd number and greater than zero.")

    kmer_size = (kmer-1)//2
    # df.loc[:, 'Dwell time'] = np.log10(df['Dwell time'])
    # df.loc[:, 'Dwell time'] = stats.zscore(df['Dwell time'])
    df = df[(df['Position'] >= position-kmer_size) & (df['Position'] <= position+kmer_size)]
    # df = df.sort_values(['Read ID', 'Position'], ascending=True)
    # df = df.reset_index(drop=True)
    df.loc[:, 'Dwell time'] = np.log10(df['Dwell time'])
    #df['Dwell time'] = np.log10(df['Dwell time'])
    # df.loc[:, 'Dwell time'] = stats.zscore(df['Dwell time'])
    # df.loc[:, 'STD'] = stats.zscore(df['STD'])
    grouped_df = df.groupby('Read ID')
    # df = df.apply(stats.zscore, axis=0)

    result_list=[]
    label_list=[]
    for key,temp in grouped_df:
        item = temp[['Mean','STD','Median','Dwell time']].values
        item = item.reshape(-1,).tolist()
        if len(item) < kmer * 4:
            continue
        if len(item) > kmer * 4:
            print(1)
        result_list.append(item)
        label_list.append([temp['Group'].values[0]])
    feature_matrix = pd.DataFrame(result_list)
    label = pd.DataFrame(label_list)
    return feature_matrix, label


def process_group(read_df, shift_size):
    # fill nan row
    read_df.reset_index(drop=True, inplace=True)
    positions = read_df['Position'].astype(int)
    missing_positions = read_df[positions.diff() > 1].index

    # 构建空行 DataFrame
    empty_df = pd.DataFrame(index=missing_positions, columns=read_df.columns)

    # 合并空行和原始数据
    merged_df = pd.concat([read_df, empty_df]).sort_index()

    data_df = merged_df[['Mean', 'STD', 'Median', 'Dwell time']]
    previous_row = data_df.shift(shift_size)
    next_row = data_df.shift(-shift_size)

    # 将当前行、前一行和后一行的数据合并成新的一行
    merged_data = pd.concat([previous_row, data_df, next_row], axis=1)
    merged_data.columns = range(merged_data.shape[1])
    merged_data['Position'] = read_df['Chrom'] + ":" + read_df['Position'].astype(str) + ':' + read_df['Strand']
    # 删除空的列并重塑成新的一行
    merged_data = merged_data.dropna(axis=0).reset_index(drop=True)
    return merged_data.values.tolist()

def extract_feature_parallel(df, kmer_size, label, cpu):
    def is_odd(number):
        if number % 2 == 1:
            return True
        else:
            return False

    if not is_odd(kmer_size) or kmer_size < 0:
        raise Exception("The kmer should be an odd number and greater than zero.")

    window_size = (kmer_size - 1) // 2
    df_group = df.groupby('Read ID')

    results_list=[]
    pool = multiprocessing.Pool(processes=cpu)
    pbar = tqdm(total=df_group.ngroups, position=0, leave=True,unit='reads')

    def update(*a):
        pbar.update(1)

    for key, read_df in df_group:
        result_per_read = pool.apply_async(process_group, args=(read_df,window_size,),callback=update)
        results_list.append(result_per_read)
    pool.close()
    pool.join()

    final_df = []
    for item in results_list:
        item = item.get()
        if item is not None:
            final_df.extend(item)

    pbar.close()
    print("Merging the result table ...")
    final_df = pd.DataFrame(final_df)
    columns = list(range(final_df.shape[1]-1))
    columns.append('Position')
    final_df.columns = columns
    final_df['Label'] = label
    return final_df

def extract_feature(df, kmer_size, label):
    def is_odd(number):
        if number % 2 == 1:
            return True
        else:
            return False

    if not is_odd(kmer_size) or kmer_size < 0:
        raise Exception("The kmer should be an odd number and greater than zero.")

    shift_size = (kmer_size - 1) // 2
    df_group = df.groupby('Read ID')

    pbar = tqdm(total= df_group.ngroups, position=0, leave=True)

    final_df = pd.DataFrame()
    for key, read_df in df_group:
        # fill nan row
        read_df.reset_index(drop=True, inplace=True)
        positions = read_df['Position'].astype(int)
        missing_positions = read_df[positions.diff() > 1].index

        # 构建空行 DataFrame
        empty_df = pd.DataFrame(index=missing_positions, columns=read_df.columns)

        # 合并空行和原始数据
        merged_df = pd.concat([read_df, empty_df]).sort_index()

        data_df = merged_df[['Mean', 'STD', 'Median', 'Dwell time']]
        previous_row = data_df.shift(shift_size)
        next_row = data_df.shift(-shift_size)

        # 将当前行、前一行和后一行的数据合并成新的一行
        merged_data = pd.concat([previous_row, data_df, next_row], axis=1)
        merged_data.columns = range(merged_data.shape[1])
        merged_data['Position'] = read_df['Chrom'] + ":" + read_df['Position'] + ':' + read_df['Strand']
        # 删除空的列并重塑成新的一行
        merged_data = merged_data.dropna(axis=0).reset_index(drop=True)
        final_df = pd.concat([final_df, merged_data], axis=0)
        pbar.update(1)
    pbar.close()
    final_df['Label'] = label
    return final_df

def save_fasta_dict(fasta_dict, path):
    f = open(path, 'w+')
    for key, value in fasta_dict.items():
        f.write('>' + key + '\n')
        line = len(value) // 80 + 1
        for i in range(0, line):
            f.write(value[i * 80:(i + 1) * 80] + '\n')
    f.close()

def calculate_MANOVA_result(position,df,length_size,subsample_num=500,windows_len=10,kmer=3):
    print("Start to run the MANOVA analysis on target region ...")
    from statsmodels.multivariate.manova import MANOVA
    from sklearn.decomposition import PCA
    # import umap
    kmer_size =(kmer-1)//2
    methylation_list=list(range( position - length_size + kmer_size ,position + length_size +1-kmer_size))
    result_list=[]
    for item in methylation_list:
        # subsample the reads
        control = df[df['Group']=='Control']
        # if control.shape[0] > subsample_num*(2*windows_len+1):
        #     control=control.iloc[0:subsample_num*(2*windows_len+1), :]
        sample = df[df['Group'] == 'Sample']
        # if sample.shape[0] > control.shape[0] * 2:
        #     sample = sample.iloc[0:control.shape[0],:]
        df = pd.concat([control,sample],axis=0).reset_index(drop=True)
        feature,label = extract_kmer_feature(df,kmer,item)
        pca = PCA(n_components=2,whiten=True)
        try:
            new_df = pd.DataFrame(pca.fit_transform(feature))
        except Exception as e:
            print(1)
        # reducer = umap.UMAP(n_components=2)  # Create a UMAP object with 2 dimensions
        # new_df = reducer.fit_transform(feature)
        new_df = pd.concat([pd.DataFrame(new_df),label], axis =1)
        new_df.columns=['PC1','PC2','Group']
        if np.sum(new_df['Group']=='Sample') > 10 and np.sum(new_df['Group']=='Control') > 10:
            manova = MANOVA.from_formula('PC1 + PC2 ~ Group', data=new_df)
            # 执行多元方差分析
            results = manova.mv_test()
            pvalue = results.summary().tables[3].iloc[0,5]
            mean_differ = feature[label[0]=='Sample'][4].median() - feature[label[0]=='Control'][4].median()
            result_list.append([item,pvalue,mean_differ])
        else:
            result_list.append([item, None])
    new_df = pd.DataFrame(result_list)
    new_df[1] = np.log10(new_df[1]) * (-1)
    new_df.columns = ['Position', 'P value(-log10)','Norm_differ']
    return new_df

# fasta=read_fasta_to_dic("../example/data/23S_rRNA.fasta")
# for key,value in fasta.items():
#     fasta[key]=reverse_fasta(value)
# save_fasta_dict(fasta,'../example/data/23S_rRNA_re.fasta')
