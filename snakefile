shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")


FILES = json.load(open('./samples.json'))
SAMPLES = sorted(FILES.keys())

TARGETS = []
final_snap_file = expand("02_snap/{sample}_2nd.snap", sample = SAMPLES)
# final_snap_file_bm_log = expand("02_snap/{sample}_add_bm.log", sample = SAMPLES)


final_fragment_file = expand("04_fragment/{sample}.bed.gz", sample = SAMPLES)


TARGETS.extend(final_snap_file)
# TARGETS.extend(final_snap_file_bm_log)

TARGETS.extend(final_fragment_file)
localrules: all, barcode_count

rule all:
	input: TARGETS

# baw_index='/home/gb148/almanlab/genome_ref/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa'
#baw_index='/home/gb148/almanlab/genome_ref/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'
#genome_size='./hg38.chrom.sizes'
#genome_name='hg38'

configfile: "config.yaml"
    
baw_index= config['baw_index']
genome_size= config['genome_size']
genome_name= config['genome_name']

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


rule 1st_align:
	input :
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		temp("01_bam/{sample}_1st.bam") ## temp bam file
	threads: 24
	shell:
		"""
		module load BWA
		module load samtools
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
		  --tmp-folder=./tmp \
		  --overwrite=TRUE
	    """		

rule 1st_snap:
	input :
		"01_bam/{sample}_1st.bam"
	output:
		temp("02_snap/{sample}_1st.snap")
	shell:
		"""
		module load BWA
		module load samtools
		snaptools snap-pre  \
		  --input-file={input}  \
		  --output-snap={output}  \
		  --genome-name={genome_name}  \
		  --genome-size={genome_size}  \
		  --min-mapq=10  \
		  --min-flen=0  \
		  --max-flen=1000  \
		  --keep-chrm=TRUE  \
		  --keep-single=TRUE  \
		  --keep-secondary=False  \
		  --overwrite=True  \
		  --max-num=10000000  \
		  --min-cov=100  \
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
		" sed '1d' {input}  > {output}"

rule update_bd:
	input: "03_barcode_info/{sample}.barcode_count"
	output: 
			sum = "03_barcode_info/{sample}.barcode_final_summary",
			map = "03_barcode_info/{sample}.barcode_final_map",
			log = "03_barcode_info/{sample}.barcode_log",
	script:
		"script/update_barcode.py"

rule r1:
	input :
		lambda wildcards: FILES[wildcards.sample]['R1'],
		"03_barcode_info/{sample}.barcode_final_map"
	output :
		"updated/{sample}_L001_R1_001.fastq"
	script:
		"script/update_fastq.py"
		
rule r1_zip:
	input :
		"updated/{sample}_L001_R1_001.fastq"
	output :
		"updated/{sample}_L001_R1_001.fastq.gz"
	threads: 8
	shell:
		"pigz -p {threads} {input}"


rule r2:
	input :
		lambda wildcards: FILES[wildcards.sample]['R2'],
		"03_barcode_info/{sample}.barcode_final_map"
	output :
		"updated/{sample}_L001_R2_001.fastq"
	script:
		"script/update_fastq.py"		

rule r2_zip:
	input :
		"updated/{sample}_L001_R2_001.fastq"
	output :
		"updated/{sample}_L001_R2_001.fastq.gz"
	threads: 8
	shell:
		"pigz -p {threads} {input}"


rule 2nd_align:
	input :
		r1 = "updated/{sample}_L001_R1_001.fastq.gz",
		r2 = "updated/{sample}_L001_R2_001.fastq.gz"
	output:
		("01_bam/{sample}_2nd.bam") 
	threads: 24
	shell:
		"""
		module load BWA
		module load samtools
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
		  --tmp-folder=./tmp \
		  --overwrite=TRUE
	    """		
	    
rule 2nd_snap_pre:
	input :
		"01_bam/{sample}_2nd.bam"
	output:
		("02_snap/{sample}_2nd.snap") 
	benchmark:
		"benchmarks/{sample}.2nd_align.txt"
	shell:
		"""
		module load BWA
		module load samtools
		snaptools snap-pre  \
		  --input-file={input}  \
		  --output-snap={output}  \
		  --genome-name={genome_name}  \
		  --genome-size={genome_size}  \
		  --min-mapq=10  \
		  --min-flen=0  \
		  --max-flen=1000  \
		  --keep-chrm=TRUE  \
		  --keep-single=TRUE  \
		  --keep-secondary=False  \
		  --overwrite=True  \
		  --max-num=10000000  \
		  --min-cov=100  \
		  --verbose=False
	    """	

rule snaptools_2nd_snap_add_bm:
	input :
		("02_snap/{sample}_2nd.snap") 
	output:
		"02_snap/{sample}_add_bm.log"
	shell:
		"""
		snaptools snap-add-bmat  \
			--snap-file={input}  \
			--bin-size-list 5000 > {output}
		"""

rule fg_dump:
	input :
		("02_snap/{sample}_2nd.snap") 
	output:
		"04_fragment/{sample}.bed.gz"
	benchmark:
		"benchmarks/{sample}.fragment_dump.txt"
	log:
		"logs/{sample}_dump_fragment.log"
	shell:
		"""
		snaptools dump-fragment --snap-file {input} \
		--output-file {output}  --buffer-size 10000 \
		--tmp-folder ./tmp &> {log}
		"""


		
