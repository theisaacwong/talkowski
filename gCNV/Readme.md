This is a repository of small scripts intended to automate the batching, clustering, and basic QC of gCNV. Requires java 14. 
Please use the R markdown (.Rmd) file for the full pipeline. 
<pre>
        Version: 2.22

        java -jar gCNV_helper.jar [Command] [required argument(s)] {optional arguement(s)}

    addGnomadAnnotations [input] [output] [gencode_gtf] [gene_column_name]
	add gnomAD gene annotation scheme column for a given column of genes
	[input] - input file 
	[output] - output file
	[gencode_gtf] - gencode gtf file
	[gene_column_name] - name of gene column

    annotateWithGenes [mode] [gencode_gtf] [input] [output] [svtype_column_name]
	annotate a bed file with genes
	[mode] - either "any" (at least 1 bp overlap) or "strict" (10/75% DEL/DUP exon space overlap) 
	[gencode_gtf] - path to gencode gtf file
	[input] - input bed file
	[output] - output path
	[svtype_column_name] - name of column containing DEL/DUP info

    bedcluster [input] [output] {prefix} {meta_file_output}
	cluster SV calls based on reciprocal overlap
	[input] - input file
	[output] - output file
	{prefix} - optional - prefix for variant names 
	{meta_file_output} - optional - clustering traceback file 

    calculateFrequency [input] [variant_column_name] [output]
	[input]
	[variant_column_name]
	[output]

    condenseBedtoolsIntersect [input] [output] [columnsToHashOn] [columnsToMerge] [columnsToKeep]
	basically bedtoolc intersect -c but with more options. Merge successive lines based on certain columns and merge the values of specified fields
	[input] - input file
	[output] - output file
	[columnsToHashOn] - comma separated integer list to identify which rows to merge, eg "1,2,3,4"
	[columnsToMerge] - comma separated integer list to identify which column values to merge into comma separated lists, eg "5,6,7"
	[columnsToKeep] - comma separated integer list to identify which columns to keep the first instance of
	
    convertToEnsemble [input] [output] [geneColumnName] [gencode_gtf]
	convert gene names to ensemble IDs
	[input] - input file
	[output] - output file
	[geneColumnName] - name of column containing gene names
	[gencode_gtf] - path to gtf file
	
    convertVCFsToBEDFormat [directory] [output] [prefix_regex] [suffix_regex]
	convert gCNV VCF output files to bed format
	[directory] - directory where VCF files are located in, requires trailing "/"
	[output] - output path for final BED file
	[prefix_regex] - prefix to trim from file name, eg 'genotyped-segments-'
	[suffix_regex] - suffix used to identify VCF files, used also to trim from file name. eg '.vcf'

    countExons [input] [gene_column_name] [gencode_gtf] [output] 
	count the number of protein coding exons an SV spans
	[input] - BED file input
	[gene_column_name] - name of column containing genes
	[gencode_gtf] - gencode gtf 
	[output] - output path
	
    defragment [input] [output] [intervals]
	defragment CNV calls
	[input] - input file
	[output] - output path
	[intervals] - BED format interval list
	
    downloadSegmentsVCFs [manifest] [directory] [column_name]
	Download VCF file outputs from the main gCNV algorithm
	[manifest] - full path to sample_set_entity.tsv file
	[directory] - Directory to download files to, w/ trailing '/'
	[column_name] - name of manifest column. eg 'segment_vcfs'

    filter [input] [output]
	Apply filter metrics to CNVs
	[input] -  input
	[output] - output

    getBarcodeCounts [manifest] [directory] [column_name]
	Download read count files specified in entity file
	[manifest] - full path to sample_set_entity.tsv file
	[directory] - Directory to download files to, w/ trailing '/'
	[column_name] - name of manifest column. eg 'segment_vcfs
	
    getBEDtrack [directory] [output] [suffix]
	Read in all read count files and generate a bed track file
	[directory] - directory where read count files are located, searches recursively
	[output]- output path 
	[suffix] - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'
	
    getCountsMatrix [directory] [output] [suffix]
	Read in all read count files and generate a matrix file
	[directory] - directory where read count files are located, searches recursively
	[output]- output path 
	[suffix] - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'

    getCountsMatrixBuffered [directory] [output] [suffix] [buffer_size]
	Read in all read count files and generate a matrix file
	[directory] - directory where read count files are located, searches recursively
	[output]- output path 
	[suffix] - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'
	[buffer_size] - number of lines to store in memory for each thread before writing

    getPerSampleMetrics [input] [column_name] [output]
	Sum integer columns by classification of a column
	[input] - input file
	[column_name] - name of column to factor by, eg "sample" or "name"
	[output] - output path

    jointCallVCF [input] [output] [variant_column] [sample_column] [sep] [svtype_column]
	convert gcnv output to joint call format
	[input] - input 
	[output]- output
	[variant_column] - name of variant column
	[sample_column] - name of sample column
	[sep] - delimiter to split samples by when writing output file, eg ","
	[svtype_column] - name of the svtype column, eg "svtype"
	
    labelMedianQS [input] [variantColumn] [qsColumn] [output]
	Add a column for the median QS value
	[input] - input
	[variantColumn] - column name to aggregate by
	[qsColumn]  - name of column to calculate median with
	[output] - output path

    simplify [input] [output] [column_to_aggregate_by] [columns_to_aggregate] [column_to_split_by]
	generate one row per unique value of a specified column, merging unique values in other columns
	[input] - input file
	[output] - output path
	[column_to_aggregate_by] - name of column to aggregate by, eg "sample"
	[columns_to_aggregate] - comma separated list of column numbers to aggregate, starting from 0
	[column_to_split_by] - split aggregation into categories, eg "svtype"
	
    subsetAnnotations [input] [main_output] [metadata_output] [gene_column_name] [gene_list_1] ... [gene_list_n]
	add a column for each gene list containing gene hits. write a per sample metadata file
	[input] - input
	[main_output] - output path
	[metadata_output]- output path
	[gene_column_name] - name of gene column to subset
	[gene_list_1] ... [gene_list_n] - space separated list of gene lists

    transposeTSV [input] [output]
	transpose a tsv file
	[input]
	[output]

    validateSubsetAnnotations [gencode_gtf] [gene_list_1] ... [gene_list_n]
	check that genes from list to be subsetted were originally labeled for in the source gtf file
	[gencode_gtf] - gencode gtf file
	[gene_list_1] ... [gene_list_n] - space separated list of gene lists
</pre>
