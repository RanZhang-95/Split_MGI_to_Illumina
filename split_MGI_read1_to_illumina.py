
#This function parse through the sample indexes file and returns a dictionary in a format like: {'Sample_A_CITE': ['ATTACTTC', 'TGCGAACT', 'GCATTCGG', 'CAGCGGAA'], ......}
def get_sample_barcode(sample_barcode_file):
    import pandas as pd
    sample_barcode_df=pd.read_csv(sample_barcode_file,sep="\t",header=[0],index_col=0)
    sample_barcode_df=sample_barcode_df.drop(sample_barcode_df.columns[[5]],axis=1) #removethe last column since it's all nan
    sample_barcode_df = sample_barcode_df.drop(["10X Index"],axis=1)
    samples=sample_barcode_df.index.tolist()
    sample_barcode_dict={}
    for sample in samples:
        sample_row=sample_barcode_df.loc[[sample]]
        sample_barcode_list=sample_row.values.flatten().tolist()
        sample=sample.replace(" ", "_")
        sample_barcode_dict[sample]=sample_barcode_list

    return (sample_barcode_dict)


#This function parse throught the record_barcode file, which is a output file from <split_read1_MGI.py>.
#This file records the sequences with its sample indexes and its illumina header from parsing through read2 file.
#As output, this function returns a dictionary that records which sequence to split and header for read1 file.
#The dictionary has format like:{tile:x_pos:y_pos': ['Sample_B_CITE', 'MGISEQ-2000:1:V300037074:3:0:1:1 1:N:0:ACCCGACG\n'],...}
#The key of the dictionary is the tile number which should be the same for the corresponding read1 and 2 sequences, the values are sample and illumina header for read1
def get_paired_header(header_from_read2_file):
    read1_header_dict={}
    with open(header_from_read2_file) as header_from_read2:
        for line in header_from_read2:
            splited_forread_line = (line.split(" "))

            splited_forread_line_front = splited_forread_line[0]
            read2_pos_info = splited_forread_line_front.split(":")
            read2_pos = read2_pos_info[4] + ":" + read2_pos_info[5] + ":" + read2_pos_info[6]
            splited_forread_line_back = splited_forread_line[1]
            barcode_info = splited_forread_line_back.split(":")
            barcode = barcode_info[-1]
            barcode = barcode.replace("\n", "")
            splited_forread_line_back = splited_forread_line_back.replace("2", "1")
            illumina_read1 = splited_forread_line_front + " " + splited_forread_line_back



            read1_info=[]

            for sample in sample_barcode_dict:
                if barcode in sample_barcode_dict[sample]:
                    read1_sample=sample
                    read1_info.append(read1_sample)
                    read1_info.append(illumina_read1)
                    read1_header_dict[read2_pos]=read1_info


    return(read1_header_dict)

#
#
#
#
#

def split_fastq_for_sample_barcodes(path_to_splitted_fastq, read1_file):
    from Bio.Seq import Seq
    from Bio import SeqIO
    from gzip import open as gzopen


    for record in SeqIO.parse(gzopen(read1_file, "rt"), format="fastq"):
        BGI_header = record.id
        tile_fastq = int(BGI_header[20:-2])
        tile_fastq = str(tile_fastq)
        x_pos = int(BGI_header[13:16])
        x_pos=str(x_pos)
        y_pos = int(BGI_header[17:20])
        y_pos=str(y_pos)
        read1_pos=tile_fastq+":"+x_pos+":"+y_pos
        try:
            read1_info= read1_header_dict[read1_pos]
            sample=read1_info[0]
            read1_header=read1_info[1]
            record.description = record.description.replace(record.id, "")
            record.id=read1_header
            splitted_fastq_file=path_to_splitted_fastq + sample + "_read1.fq"
            with open(splitted_fastq_file, "a") as output_handle:
                SeqIO.write(record, output_handle, "fastq")
        except:
            continue





#
# record_barcode_file="/Users/ranzhang/Documents/Garvan/Projects/MGI/record_barcode.txt"
# read1_fastq="/directflow/GWCCGPipeline/projects/sequencing/MGI/L04/V300037074_L04_read_1.fq.gz"
# read_position_for_record(record_barcode_file,read1_fastq)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-sample_barcode_file', type=str,help="please input sample barcode file")
    parser.add_argument('-path_to_splitted_fastq', type=str,help="please input target path")
    parser.add_argument('-fastq_file', type=str,help="please input text file")
    parser.add_argument('-read2_header', type=str, help="please input text file")
    args = parser.parse_args()
    sample_barcodes=args.sample_barcode_file
    path=args.path_to_splitted_fastq
    fastq=args.fastq_file
    read2_header_file=args.read2_header
    print('successfully getting arguments:',sample_barcodes,path,fastq, read2_header_file)
    sample_barcode_dict=get_sample_barcode(sample_barcodes)
    read1_header_dict=get_paired_header(read2_header_file)
    split_fastq_for_sample_barcodes(path, fastq)





