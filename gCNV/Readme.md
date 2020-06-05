This is a repository of small scripts intended to automate the batching, clustering, and basic QC of gCNV. 
Please use the R markdown (.Rmd) file for the full pipeline. 

        Version: 2.17

        java -jar gCNV_helper.jar [Command] [required argument(s)] {optional arguement(s)}

        getBarcodeCounts [entityPath] [working-directory] {counts-field-name}
                Download the read count files specified in the entity file
                [entityPath] - Full path to the entity file. eg: '/home/gCNV/sample_set_entity.tsv'
                [working-directory] - Directory to download files to. eg '/home/gCNV/'
                {counts-field-name} - optional - The name of the column in the entity file containing the counts paths. eg 'output_counts_barcode'

        getCountsMatrix [sourceFolder] [OUTPUT_PATH] {regex}
                Read in all read count files and generate a matrix file
                [sourceFolder] - Directory where read count files are located in, files can be in sub-directories.
                [OUTPUT_PATH] -  The full output path where the matrix file will be written to
                {regex} - optional - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'

        getCountsMatrixBuffered [sourceFolder] [OUTPUT_PATH] [regex] {buffer-size}
                Read in all read count files and generate a matrix file
                [sourceFolder] - Directory where read count files are located in, files can be in sub-directories.
                [OUTPUT_PATH] -  The full output path where the matrix file will be written to
                [regex] - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'
                {buffer-size} - number of lines to store in memory for each thread before writing

        downloadSegmentsVCFs [entityPath] [working-directory] {column-name}
                Download the VCF file outputs from running the main gCNV algorithm
                [entityPath] - Full path to the entity file. eg: '/home/gCNV/sample_set_entity.tsv'
                [working-directory] - Directory to download files to. eg '/home/gCNV/'
                {column-name} - optional - The name of the column in the entity file containing the counts paths. eg 'segments_vcfs'

        convertVCFsToBEDFormat [working-directory] [output-path] {prefix-regex} {suffix-regex}
                Convert the gCNV VCF output files to BED format for svtk bedlcuster input
                [working-directory] - Directory where VCF files are located in, files can be in sub-directories.
                [ouput-path] - The output path for the final consolidated BED file
                {prefix-regex} - prefix to trim from file name, eg 'genotyped-segments-'
                {suffix-regex} - suffix used to identify VCF files, used also to trim from file name. eg '.vcf'

        svtkMatch [svtk_input] [svtk_output] [output_path]
                Match up the gCNV meta data with the svtk bedcluster meta data and write to file
                [svtk_input] - The BED file that was given to svtk bedcluster
                [svtk_output] - The output file from svtk bedcluster
                [output_path] - The full path to write the output file to.
                        
        getPerSampleMetrics [input_path] [column_name] [output_path]
                Sum integer columns by classification of a column
                [input_path] - gcnv output file
                [column_name] - the name of column to factor by
                [output_path] - The full path to write the output file to.

        labelMedianQS [input_path] [variantColumn] [qsColumn] [output_path]
                Add a column for the median QS value
                [input_path] - gcnv output file
                [variantColumn] - column name to aggregate by
                [qsColumn] - the name of column to calculate median with
                [output_path] - The full path to write the output file to.

        jointCallVCF [input_path] [output_path] [variant_column] [sample_column] {sep} {svtpye_column}
                convert gcnv output to joint call format
                [input_path] - gcnv output file
                [output_path] - The full path to write the output file to.
                [variant_column] - The name of the variant column
                [sample_column] - The name of the sample column
                {sep} - seperator to split samples by when writing output file. default ','
                {svtype_column} - The name of the svtype column, default 'svtype'
                	
        annotateWithGenes [mode] [gcnv_input_path] [annotation_input_path] [output_path]
		Annotate a bed file with overlapping intervals from a second bed file
		[mode] - either 'strict' or 'any'; strict requires 10%/75% exon space overlap for DEL/DUP, any requires at least 1bp overlap
		[gcnv_input_path] - gcnv file path
		[annotation_input_path] - any: a bed file with columns: chr, start, end, name; strict: gencode gtf file uncompressed
		[output_path] - The full path to write the output file to.
		
        condenseBedtoolsIntersect [input_path] [output_path] [columns_to_hash_on] [columns_to_merge] [columns_to_keep]
		Condense the output of bedtools intersect by adding a column for annotations instead of having each row be a unique annotation
		[input_path] - input path
		[output_path] - output path
		[columns_to_hash_on] - comma separated columns to identify rows to merge. eg 1,2,3 
		[columns_to_merge] - comma separated columns to merge into annotation column. eg 7,9
		[columns_to_keep] - comma separated columns to keep, eg 4,5,6
		
        defragment [input_path] [output_path] 
		defragment gcnv calls
		[input_path] - input path
		[output_path] - output path

        validateSubsetAnnotations [gtfFile] [annotationSubsets] ... [annotationSubsets] 
		check to see if the genes you will be subsetting were contained in the original annotation file
		[gtfFile] - gencode gtfFile
		[annotationSubsets] - space separated list of gene lists

        subsetAnnotations [gcnvInput] [output] [sourceColumnName] [annotationSubsets] ... [annotationSubsets]
		add a column for each gene list for annotated genes contained in said list
		[gcnvInput] - input file
		[output] - output path
		[sourceColumnName] - column name of annotated genes
		[annotationSubsets] - space separated list of gene lists

