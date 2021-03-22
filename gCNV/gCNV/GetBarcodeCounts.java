package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class GetBarcodeCounts extends gCNVHelperTool {

	public static String[] inputs = new String[] {COLUMN_NAME, DIRECTORY, MANIFEST_PATH};

	public GetBarcodeCounts(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException, InterruptedException {
		String entityPath = args.get(MANIFEST_PATH);
		String wd = args.get(DIRECTORY);
		String output_counts_barcode_regex = args.get(COLUMN_NAME);
		this.getBarcodeCounts(entityPath, wd, output_counts_barcode_regex);
	}
	

	/**
	 * Downloads the barcode counts file from terra/firecloud
	 * 
	 * @param entityPath                  - full path to the entity file, eg
	 *                                    "/documents/sample_set_entity.tsv"
	 * @param wd                          - the directory to download the counts
	 *                                    files to
	 * @param output_counts_barcode_regex - the name of the column containing the
	 *                                    counts files, eg "output_counts_barcode"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getBarcodeCounts(String entityPath, String wd, String output_counts_barcode_regex)
			throws IOException, InterruptedException {
		DataFrame entityDF = new DataFrame(entityPath, true, "\\t", "@");
		System.out.println(entityDF.nrow());
		System.out.println(entityDF.columnMapping.toString());
		for (int i = 0; i < entityDF.nrow(); i++) {
			String individualToDownload = entityDF.get("entity:sample_set_id", i);
			String pathInGoogleBucket = entityDF.get(output_counts_barcode_regex, i).replaceAll("\\[|\\]|\"", "")
					.replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");

			if (files.length == 1 && (files[0].trim().equals("") || files[0].trim().equals("0"))) {
				continue;
			}

			String temp_suffix = "_temp_files_to_download.txt";

			File file = new File(wd + individualToDownload + temp_suffix);
			BufferedWriter output = new BufferedWriter(new FileWriter(file));
			for (String filePath : files) {
				output.write(filePath + "\n");
			}
			output.close();

			String pathToDownloadTo = wd + individualToDownload + "/";
			File dirToDownloadTo = new File(pathToDownloadTo);
			if (!dirToDownloadTo.exists()) {
				dirToDownloadTo.mkdir();
			}

			String wdUnix = wd.replace("C:", "/mnt/c");
			String pathToDownloadToUnix = wdUnix + individualToDownload;

			downloadFiles(wdUnix + individualToDownload + temp_suffix, pathToDownloadToUnix);

		}
	}
	
}
