package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class DownloadFilteredIntervals extends gCNVHelperTool {

	public static String[] inputs = new String[] {MANIFEST_PATH, DIRECTORY, COLUMN_NAME, OUTPUT_PATH};

	public DownloadFilteredIntervals(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException, InterruptedException {
		String entityPath = args.get(MANIFEST_PATH);
		String wd = args.get(DIRECTORY);
		String filtered_intervals_columnName = args.get(COLUMN_NAME);
		String outputPath = args.get(OUTPUT_PATH);
		this.downloadFilteredIntervals(entityPath, wd, filtered_intervals_columnName, outputPath);
	}
	

	/**
	 * @param entityPath               - full path to the entity file, eg
	 *                                 "sample_set_entity.tsv"
	 * @param wd                       - working directory, where the files should
	 *                                 be downloaded to
	 * @param segments_vcfs_columnName - the name of the column in the entity file
	 *                                 containing the segment_vcfs, eg
	 *                                 "segments_vcfs"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void downloadFilteredIntervals(String entityPath, String wd, String segments_vcfs_columnName, String outputPath)
			throws IOException, InterruptedException {
		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
		
		
		for (int i = 0; i < sampleSetEntity.nrow(); i++) {
			String individualToDownload = sampleSetEntity.get("entity:sample_set_id", i);
			String pathInGoogleBucket = sampleSetEntity.get(segments_vcfs_columnName, i).replaceAll("\\[|\\]|\"", "")
					.replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");

			// if(pathInGoogleBucket.trim().equals("") || (files.length==1 &&
			// files[0].trim().equals("")) || files.length==0) {continue;}
			if (files.length == 1 && (files[0].trim().equals("") || files[0].trim().equals("0"))) {
				continue;
			}

			String temp_suffix = "_files_to_download.txt";

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

			// //this is only tested on windows, and it works so far. ITS VERY DELICATE.
			String[] command;

			// unzip vcfs, working as of 12/3/19
			if (System.getProperty("os.name").contains("Windows")) {
				String[] dosCommand = { "bash", "-c", "'gunzip", "-k", pathToDownloadToUnix + "/*'" };
				command = dosCommand;
			} else {
				String[] unixCommand = { "gunzip", "-k" + pathToDownloadToUnix + "/*'" };
				command = unixCommand;
			}
			System.out.println(String.join(" ", command));
			try {
				new ProcessBuilder(command).inheritIO().start().waitFor();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

}
