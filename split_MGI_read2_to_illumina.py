#This python script takes a sample barcode text file and read2 fastq.gz file, convert the  MGI header of the fastq file and split records into corresponding sample files.
#Noted that the sample barcodes are located at the beggining of the reverse complement read2 sequence
#
#
#
#
#
#

# #This function takes a sample barcode text file and returns a sample barcode dictionary
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



def split_fastq_for_sample_barcodes(read2_file,path_to_splitted_fastq,dict):
    from Bio.Seq import Seq
    from Bio import SeqIO
    from gzip import open as gzopen

    samples=list(dict.keys())
    sample_total=[]
    fastqs=[]
    for i in samples:
        sample_total.append(dict[i])
        fastqs.append(path_to_splitted_fastq+str(i)+"_read2.fq")

    
    record_barcode_file=path_to_splitted_fastq+"record_barcode.txt"
    no_barcode=path_to_splitted_fastq+"no_barcode.fq"

    for record in SeqIO.parse(gzopen(read2_file,"rt"), format="fastq"):
        sequence=str(record.seq)
        barcode_in_fastq = Seq(sequence[:])
        barcode_in_fastq = barcode_in_fastq.reverse_complement()
        barcode_in_fastq=barcode_in_fastq[0:8]
        findornot = False
        for i in range(len(sample_total)):
            if barcode_in_fastq in sample_total[i]:
                findornot=True
                BGI_header = record.id
                flowcell_id = BGI_header[0:10]
                lane = BGI_header[11:12]
                tile = int(BGI_header[20:-2])
                x_pos = int(BGI_header[13:16])
                y_pos = int(BGI_header[17:20])
                tile=str(tile)
                x_pos=str(x_pos)
                y_pos=str(y_pos)
                read_no = BGI_header[-1]
                sample_barcode =str(barcode_in_fastq)
                illumina_header = "MGISEQ-2000:1:" + flowcell_id + ":" + lane + ":" + tile + ":" + x_pos + ":" + y_pos + " " + read_no + ":N:0:" + sample_barcode
                record.description=record.description.replace(record.id, "")
                record.id=illumina_header
                with open(fastqs[i], "a") as output_handle:
                    SeqIO.write(record, output_handle, "fastq")
                with open(record_barcode_file,"a") as record_barcode:
                    record_barcode.write(illumina_header+"\n")



        if findornot==False:
            print(barcode_in_fastq)
            print(record)
            print(type(record))
            with open(no_barcode, "a") as recording_no_barcode:
                print(type(record))
                SeqIO.write(record, recording_no_barcode, "fastq")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-sample_barcode_file', type=str,help="please input sample barcode file")
    parser.add_argument('-path_to_splitted_fastq', type=str,help="please input target path")
    parser.add_argument('-fastq_file', type=str,help="please input text file")
    args = parser.parse_args()
    sample_barcodes=args.sample_barcode_file
    path=args.path_to_splitted_fastq
    fastq=args.fastq_file
    print('successfully getting arguments:',sample_barcodes,path,fastq)
    dict=get_sample_barcode(sample_barcodes)
    split_fastq_for_sample_barcodes(fastq,path,dict)




    #python split_samples_L03_change_header_soft.py -sample_barcode_file /directflow/SCCGGroupShare/projects/ranzha/MGI/MGI_sequencing_sample_indexes.txt --path_to_splitted_fastq /directflow/SCCGGroupShare/projects/ranzha/MGI/L03_2_s/ -fastq_file /home/ranzha/MGI/V300037074_L03_read_2.fq.gz