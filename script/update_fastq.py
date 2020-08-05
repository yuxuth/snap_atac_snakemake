import gzip
import itertools

# def from_file_to_barcode_list(filename):
#     '''
#     Load the barcode whitelist into the memory
#     '''
#    
#     # f = open(filename,"r")
#     f = gzip.open(filename,'rt')
#     barcodes = []
#     while(True):
#         line = f.readline().rstrip('\n')
#         if not line:
#             break
#         barcodes.append(line)
#     f.close()
#     print()
#     return set(barcodes)


def newbarcode_white_list_dic(barcode_map):
    barcode_dic = {}
    with open(barcode_map,'r') as f:
        for line in f:
            (mismatch_barcode,white_list_barcode) = line.rstrip('\n').strip().split(" ")
            barcode_dic[mismatch_barcode] = white_list_barcode
    return barcode_dic

# def update_fastq(r1,r2,out_r1,out_r2, barcode_dic ): ## process two files
def update_fastq(r1,out_r1, barcode_dic ): ## process one files
    """"
    modify the barcode
    """

    # fastq_reader1 = HTSeq.FastqReader(r1)
    # fastq_reader2 = HTSeq.FastqReader(r2)
    # 
    # output1 = gzip.open(out_r1, "wt" )
    # output2 = gzip.open(out_r2, "wt" )

    """
    a = zip(fastq_reader1,fastq_reader2 )
    read1,read2 = next(a)
    """

    # for read1,read2 in zip(fastq_reader1,fastq_reader2 ):
    #     cur_barcode = read1.name.split(':')[0]
    #     if  cur_barcode in barcode_dic:
    #         new_read_name = (barcode_dic[cur_barcode]+ read1.name[32:]) ## 32bp cell barcode
    #         read1.name = new_read_name
    #         read2.name = new_read_name
    #     read1.write_to_fastq_file(output1)
    #     read2.write_to_fastq_file(output2)
    

    # for read1,read2 in zip(fastq_reader1,fastq_reader2 ):
    #     try:
    #         # cur_barcode = read1.name.split(':')[0]
    #         # if  cur_barcode in barcode_dic:
    #         #     new_read_name = (barcode_dic[cur_barcode]+ read1.name[32:]) ## 32bp cell barcode
    #             # read1.name = new_read_name
    #             # read2.name = new_read_name
    #         # read1.write_to_fastq_file(output1)
    #         # read2.write_to_fastq_file(output2)
    #     except StopIteration:
    #         pass
    # r1,r2,out_r1,out_r2, barcode_set
    
    f_r1 = gzip.open(r1, 'rt')
    # f_r2 = gzip.open(r2, 'rt')
#     f_out_r1 = gzip.open(out_r1, 'wt') ## using pigz for zip
    f_out_r1 = open(out_r1, 'w') 
    # f_out_r2 = gzip.open(out_r2, 'wt')

    while True:
        cur_r1_name = f_r1.readline().strip()[1:]
        cur_r1_read = f_r1.readline().strip()
        cur_r1_plus = f_r1.readline().strip()
        cur_r1_qual = f_r1.readline().strip()
        # 
        # cur_r2_name = f_r2.readline().strip()[1:]
        # cur_r2_read = f_r2.readline().strip()
        # cur_r2_plus = f_r2.readline().strip()
        # cur_r2_qual = f_r2.readline().strip()
    
        # if cur_r1_name == "" or cur_r2_name == "" : break
        if cur_r1_name == "" : break
        
        # if not (cur_r1_name.split()[0] == cur_r2_name.split()[0] ): sys.exit("error: read name does not match")
    
        cur_barcode = (cur_r1_name[:32].upper())
        if  cur_barcode in barcode_dic:
            cur_r1_name = (barcode_dic[cur_barcode]+ cur_r1_name[32:]) ## 32bp cell barcode
            # cur_r2_name = (barcode_dic[cur_barcode]+ cur_r2_name[32:])
        
        f_out_r1.write('@' + cur_r1_name +"\n")
        f_out_r1.write(cur_r1_read+"\n")
        f_out_r1.write(cur_r1_plus+"\n")
        f_out_r1.write(cur_r1_qual+"\n")     
    
        # f_out_r2.write('@' + cur_r1_name +"\n")
        # f_out_r2.write(cur_r1_read+"\n")
        # f_out_r2.write(cur_r1_plus+"\n")
        # f_out_r2.write(cur_r1_qual+"\n")    


    f_r1.close()
    # f_r2.close()
    f_out_r1.close()
    # f_out_r2.close()

r1 = str(snakemake.input[0])
# r2 = str(snakemake.input['r2'])
barcode_dic = newbarcode_white_list_dic(snakemake.input[1])
out_r1 = snakemake.output[0]
# out_r2 = snakemake.output['r2']
# barcode_set = from_file_to_barcode_list('barcode_whitelist.txt.gz') ## provide in the snakefile current location

# print(type(r1))

# update_fastq(r1,r2,out_r1,out_r2,barcode_dic )
update_fastq(r1,out_r1,barcode_dic )

