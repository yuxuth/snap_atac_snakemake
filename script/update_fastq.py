import HTSeq
import gzip
import itertools

def from_file_to_barcode_list(filename):
    '''
    Load the barcode whitelist into the memory
    '''
   
    # f = open(filename,"r")
    f = gzip.open(filename,'rt')
    barcodes = []
    while(True):
        line = f.readline().rstrip('\n')
        if not line:
            break
        barcodes.append(line)
    f.close()
    print()
    return set(barcodes)


def newbarcode_white_list_dic(barcode_map):
	barcode_dic = {}
	with open(barcode_map,'r') as f:
	    for line in f:
	        (mismatch_barcode,white_list_barcode) = line.rstrip('\n').trip().split(" ")
	        barcode_dic[barcode_dic] = white_list_barcode
    return barcode_dic

def update_fastq(r1,r2,barcode_dic,out_r1,out_r2, barcode_set ): ## process two files
    """
    modify the barcode
    """
    fastq_reader1 = HTSeq.FastqReader(r1)
    fastq_reader2 = HTSeq.FastqReader(r2)

    output1 = gzip.open(out_r1, "wt" )
    output2 = gzip.open(out_r2, "wt" )
    
    """
    a = zip(fastq_reader1,fastq_reader2 )
    read1,read2 = next(a)
    """

    for read1,read2 in zip(fastq_reader1,fastq_reader2 ):
    	cur_barcode = read1.name.split(':')[0]
    	if not cur_barcode in barcode_set:
    		new_read_name = (barcode_dic[cur_barcode]+ read1.name[32:]) ## 32bp cell barcode
    		read1.name = new_read_name
    		read2.name = new_read_name
        read1.write_to_fastq_file(output1)
        read2.write_to_fastq_file(output2 )
    
    output1.close()
    output2.close()


r1 = snakemake.input['r1']
r2 = snakemake.input['r2']
barcode_dic = newbarcode_white_list_dic(snakemake.input['barcode_info'])
out_r1 = snakemake.output['r1']
out_r2 = snakemake.output['r2']
barcode_set = from_file_to_barcode_list('barcode_whitelist.txt.gz') ## provide in the snakefile current location

update_fastq(r1,r2,barcode_dic,out_r1,out_r2, barcode_set )

