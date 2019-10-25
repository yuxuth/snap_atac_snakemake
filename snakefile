shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")


FILES = json.load(open('./samples.json'))
SAMPLES = sorted(FILES.keys())

TARGETS = []
final_snap_file = expand("02_snap/{sample}_2nd.snap", sample = SAMPLES)
TARGETS.extend(final_snap_file)

rule all:
	input: TARGETS

baw_index='/home/gb148/almanlab/genome_ref/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa'
genome_size='../mm10.chrom.sizes'
genome_name='mm10'

## start from the demutiplexed fastqs
# 1. using the snaptools to finish the mapping, extract the barcode information 
# 2. examine the barcode collision to correct the fastq file 
# 3. modify the fastq file based on the barcode hash 
# 4. rerun the snaptools to create the final snapfile
# 4.1 doublet selection to remove the potential doublet
# 5. using snapatac to generate the bin-guided clustering 
# 6. export the peak for each cluster 
# 7. update the peak matrix in the snapfile. 
# 8. recluster everything based on the peak matrix
# 9. update the cluster information


rule snaptools_1st_align:
	input :
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		("01_bam/{sample}_1st.bam") ## temp bam file
	threads: 4
	shell:
		"""
		snaptools align-paired-end  \
		  --input-reference={baw_index}   \
		  --input-fastq1={input.r1} \
		  --input-fastq2={input.r2} \
		  --output-bam={output} \
		  --aligner=bwa \
		  --read-fastq-command=zcat  \
		  --min-cov=0 \
		  --num-threads={threads}  \
		  --if-sort=True  \
		  --tmp-folder=./ \
		  --overwrite=TRUE
	    """		

rule snaptools_1st_snap:
	input :
		"01_bam/{sample}_1st.bam"
	output:
		("02_snap/{sample}_1st.snap")
	shell:
		"""
		snaptools snap-pre  \
		  --input-file={input}  \
		  --output-snap={output}  \
		  --genome-name={genome_name}  \
		  --genome-size={genome_size}  \
		  --min-mapq=30  \
		  --min-flen=0  \
		  --max-flen=1000  \
		  --keep-chrm=TRUE  \
		  --keep-single=TRUE  \
		  --keep-secondary=False  \
		  --overwrite=True  \
		  --max-num=1000000  \
		  --min-cov=0  \
		  --verbose=False
	    """	

rule barcode_dump:
	input: 
		"02_snap/{sample}_1st.snap"
	output:
		"03_barcode_info/{sample}.barcode_info"
	shell:
		"""
		snaptools dump-barcode \
			--snap-file={input} \
			--output-file={output} \
		"""

rule barcode_count:
	input: "03_barcode_info/{sample}.barcode_info"
	output: "03_barcode_info/{sample}.barcode_count"
	shell:
		" awk '{{print $3, $2}}' {input}  > {output}"

rule barcode_info_update:
	input: count = "03_barcode_info/{sample}.barcode_count"
	output: 
			sum = "03_barcode_info/{sample}.barcode_final_summary",
			map = "03_barcode_info/{sample}.barcode_final_map",
			log = "03_barcode_info/{sample}.barcode_log",
	script:
		"script/update_barcode.py"

rule update_fastq_file_r1:
	input :
		lambda wildcards: FILES[wildcards.sample]['R1'],
		"03_barcode_info/{sample}.barcode_final_map"
	output :
		"updated/{sample}_L001_R1_001.fastq.gz"
	script:
		"script/update_fastq.py"


rule update_fastq_file_r2:
	input :
		lambda wildcards: FILES[wildcards.sample]['R2'],
		"03_barcode_info/{sample}.barcode_final_map"
	output :
		"updated/{sample}_L001_R2_001.fastq.gz"
	script:
		"script/update_fastq.py"		



rule snaptools_2nd_align:
	input :
		r1 = "updated/{sample}_L001_R1_001.fastq.gz",
		r2 = "updated/{sample}_L001_R2_001.fastq.gz"
	output:
		("01_bam/{sample}_2nd.bam") ## temp bam file
	threads: 12
	shell:
		"""
		snaptools align-paired-end  \
		  --input-reference={baw_index}   \
		  --input-fastq1={input.r1} \
		  --input-fastq2={input.r2} \
		  --output-bam={output} \
		  --aligner=bwa \
		  --read-fastq-command=zcat  \
		  --min-cov=0 \
		  --num-threads={threads}  \
		  --if-sort=True  \
		  --tmp-folder=./ \
		  --overwrite=TRUE
	    """		
rule snaptools_2nd_snap_pre:
	input :
		"01_bam/{sample}_2nd.bam"
	output:
		("02_snap/{sample}_2nd.snap") ## temp bam file
	shell:
		"""
		snaptools snap-pre  \
		  --input-file={input}  \
		  --output-snap={output}  \
		  --genome-name={genome_name}  \
		  --genome-size={genome_size}  \
		  --min-mapq=30  \
		  --min-flen=0  \
		  --max-flen=1000  \
		  --keep-chrm=TRUE  \
		  --keep-single=TRUE  \
		  --keep-secondary=False  \
		  --overwrite=True  \
		  --max-num=1000000  \
		  --min-cov=0  \
		  --verbose=False
	    """	

