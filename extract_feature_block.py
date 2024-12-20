import pandas as pd
import re
import numpy as np
from cem_utils import generate_paf_file_eventalign,identify_file_path,generate_bam_file,build_out_path,base_shift_dict,read_fasta_to_dic
import pyslow5
from tqdm import tqdm
from normalization import normalize_signal_with_lim
import argparse
import multiprocessing
import time

def create_output_file(fasta,output_path,win_size=100000):
    fasta =read_fasta_to_dic(fasta)
    global result_dict
    result_dict={}
    for key,value in fasta.items():
        file_number= len(value)//win_size
        for item in range(0,file_number+1):
            file_name = key+'_'+str(item)
            raw_output_file = open(output_path + '/' + file_name+'.csv', 'w+')
            result_dict[file_name] = raw_output_file

def write_output_file(chrom,pos_list,total_feature_per_reads,win_size=100000):
    # read_id = 'd857a5e2-e055-4372-bea7-e2c3517ac13b'
    # if total_feature_per_reads[0][0] ==read_id:
    #     print(1)
    # else:
    #     return
    start,end=pos_list
    if start//win_size == end//win_size:
        file_number = start//win_size
        file_name = chrom+'_'+str(file_number)
        raw_output_file = result_dict[file_name]
        for line in total_feature_per_reads:
            line_string = ','.join(line) + '\n'
            raw_output_file.write(line_string)
    else:
        for line in total_feature_per_reads:
            file_number = int(line[2]) // win_size
            file_name = chrom + '_' + str(file_number)
            raw_output_file = result_dict[file_name]
            line_string = ','.join(line) + '\n'
            raw_output_file.write(line_string)

def extract_feature(line, base_shift, norm=True, nucleotide_type='DNA',pore='r9',win_size=100000,flip_win=10):

    read_id = line[0]
    # if read_id == 'edba3d28-323d-47ce-94c7-fa45d3751990':
    #     print(1)
    # tackle moves tag
    moves_string = line[14]
    moves_string = re.sub('ss:Z:', '', moves_string)
    moves_string = re.sub('D', 'D,', moves_string)
    moves_string = re.sub('I', 'I,', moves_string)
    # print(moves_string)
    moves = re.split(r',+', moves_string)
    moves = moves[:-1]
    # extract index and generate event_length and event_start
    insertion = 0
    event_length = []
    for i, item in enumerate(moves):
        if 'D' in item:
            deletion = int(item[:-1])
            for i in range(deletion):
                event_length.append(0)
        elif 'I' in item:
            if i == 0:
                continue
            else:
                return None
        elif '=' in item:
            return None
        else:
            event_length.append(int(item))
    # build event_length from move table
    read_event = s5.get_read(read_id, aux=["read_number", "start_mux"], pA=True)
    start_index = line[2]
    end_index = line[3]
    event_length = np.array(event_length)
    strand = line[4]

    try:
        assert end_index - start_index == np.sum(event_length)
        assert read_event['len_raw_signal'] == line[1]
        assert event_length.shape[0] == line[10]
        assert abs(line[8] - line[7]) == line[10]
    except Exception:
        print("Warning: 1 read's length of signal is not equal between blow5 and paf")
        return None

    signal = read_event['signal']
    signal = signal[start_index:end_index]

    if norm:
        signal,shift,scale = normalize_signal_with_lim(signal)

    event_starts = event_length.cumsum()
    event_starts = np.insert(event_starts, 0, 0)[:-1]

    # shift
    if base_shift == 'auto':
        base_shift = base_shift_dict[pore + nucleotide_type + strand]

    # base shift
    ref_start = np.min([line[7], line[8]])
    ref_end = np.max([line[7], line[8]])
    start_pos = ref_start - base_shift
    end_pos = ref_end - base_shift

    total_feature_per_reads = []
    for x in range(0,line[10]):
        # to filter the first and the last part of one read
        if x < flip_win or x >= line[10]-flip_win:
            continue
        element = signal[event_starts[x]:event_starts[x] + event_length[x]]
        if event_length[x] == 0:
            continue
        if (nucleotide_type == 'RNA' and line[4]=='+') or (nucleotide_type == 'DNA' and line[4]=='-'):
            position = end_pos - x
        else:
            position = start_pos + x
        temp = [read_id,line[5],str(position),line[4],str(np.mean(element).round(4)), str(np.std(element).round(4)), str(np.median(element).round(4)), str(event_length[x])]
        total_feature_per_reads.append(temp)

    if len(total_feature_per_reads) > 0:
        write_output_file(line[5],(start_pos,end_pos),total_feature_per_reads,win_size)
        # temp_df = pd.DataFrame(total_feature_per_reads)
        # temp_df.columns=['Read ID','Chrom','Position','Strand','Mean', 'STD', 'Median', 'Dwell time']
        # convert_per_read(temp_df,1)


def update(*a):
    pbar.update(1)

def read_blow5(fastq_file,slow5_file, reference, pore, subsample_ratio,output,base_shift, norm=True, cpu=4, rna=False,win_size=100000):
    global s5, pbar,raw_output_file,kmer_output_file
    identify_file_path(fastq_file)
    identify_file_path(slow5_file)
    build_out_path(output)
    new_fastq,bam_file = generate_bam_file(fastq_file, reference, str(cpu),output,subsample_ratio)
    paf_file = generate_paf_file_eventalign(new_fastq, slow5_file,bam_file,reference,pore,rna,cpu,output)
    print('Parsing the paf file ...')
    df = pd.read_csv(paf_file, sep='\t', header=None)
    # df = df.tail(100000)
    # paf_dict = df.set_index(df.columns[0]).to_dict(orient='index')
    print("There are "+str(df.shape[0])+' reads did the f5c eventalign')

    if rna:
        nucleotide_type = 'RNA'
    else:
        nucleotide_type = 'DNA'
    # kmer = kmer_model_size[pore+'+'+nucleotide_type]
    print('Parsing the blow5 file ...')
    s5 = pyslow5.Open(slow5_file, 'r')

    pbar = tqdm(total=df.shape[0], position=0, leave=True,  unit='read')
    create_output_file(reference,output,win_size)

    with multiprocessing.Pool(processes=cpu) as pool:
        for idx, row in df.iterrows():
            pool.apply_async(extract_feature, args=(row,base_shift, norm, nucleotide_type,pore,win_size), callback=update)
        pool.close()
        pool.join()
    pbar.close()
    for key,value in result_dict.items():
        value.close()
    # kmer_output_file.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq",help="input fastq file")
    parser.add_argument("-b", "--blow5", help="input blow5 file")
    parser.add_argument('-o', "--output", default='nanoSundial_feature', help="output file path")
    parser.add_argument('-r', "--ref",  help="reference path (fasta file)")
    parser.add_argument('-t', "--cpu", default=1,type=int, help="Process numbers")
    parser.add_argument("--pore", default='rna004', choices=['r9', 'r10','rna004'], help="flowcell")
    parser.add_argument("--subsample", default=1, type=float, help="subsample ratio (0-1)")
    parser.add_argument("--win_size", default=100000, type=int, help="windows size")
    parser.add_argument("--flip", default=10, type=int, help="length of base of start or end")
    parser.add_argument('--rna', action='store_true', help='Turn on the RNA mode')
    parser.add_argument('--base_shift', choices=['auto', '0', '-1', '-2', '-3', '-4', '-5', '-6', '-7', '-8'], default='auto',
                            help='base shift option')
    args = parser.parse_args()
    pore = args.pore
    if args.rna:
        nucleotide_type = "RNA"
    else:
        nucleotide_type = 'DNA'
    args.norm = True
    # args.subsample = 1
    start_time = time.time()
    base_shift = args.base_shift

    read_blow5(args.fastq,args.blow5, args.ref, args.pore, args.subsample,args.output,base_shift=base_shift,norm=args.norm, cpu=args.cpu,rna=args.rna,win_size=args.win_size)
    end_time = time.time()
    run_time = end_time - start_time
    # Convert to minutes and hours
    minutes = int(run_time // 60)  # Integer division to get the whole minutes
    seconds = int(run_time % 60)  # Remainder to get the remaining seconds

    hours = int(minutes // 60)  # Integer division to get the whole hours
    minutes %= 60  # Modulo operation to get the remaining minutes

    # Print the running time
    print("Program running time:", hours, "hours", minutes, "minutes", seconds, "seconds")