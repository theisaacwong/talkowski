
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
                [svtk_input]
                [svtk_output]
                [output_path]
