package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

public class GetCountsMatrix extends gCNVHelperTool {

	public static String[] inputs = new String[] {DIRECTORY, OUTPUT_PATH, COUNTS_REGEX};

	public GetCountsMatrix(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException, InterruptedException {
		String sourceFolder = args.get(DIRECTORY);
		String OUTPUT = args.get(OUTPUT_PATH);
		String countsRegex = args.get(COUNTS_REGEX);
		this.getCountsMatrix(sourceFolder, OUTPUT, countsRegex);
	}
	

	/**
	 * This is an improvement upon the previous version as it now more accurately
	 * forces GC to remove the DF after its been read and only keep the important
	 * field Reads in barcode counts files in parallel, and writes a counts matrix
	 * 
	 * @param sourceFolder - directory where files are located in, files can be in
	 *                     sub-directories.
	 * @param OUTPUT_PATH  - the full output path where the matrix file will be
	 *                     written to
	 * @param countsRegex  - the regex suffix to identify counts files
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrix(String sourceFolder, String OUTPUT_PATH, String countsRegex)
			throws IOException, InterruptedException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();

		System.out.print("finding files.\t");
		Path p = Paths.get(sourceFolder);
		final int maxDepth = 10;
		Stream<Path> matches = Files.find(p, maxDepth,
				(path, basicFileAttributes) -> String.valueOf(path).endsWith(countsRegex));
		matches.filter(s -> s.getFileName().toString().contains(countsRegex)).forEach(barcodeCountsPaths::add);
		matches.close();
		for (Path fp : barcodeCountsPaths) {
			barcodeCountsFiles.add(fp.toAbsolutePath().toString());
			sampleNames.add(fp.getFileName().toString().replaceAll(countsRegex, ""));
		}
		System.out.println("found " + sampleNames.size() + " files");

		System.out.println("reading files. \t");

		// toRead functions as a synchronized queue so each thread knows which file to
		// read next
		// each dataframe is then mapped to a unique index number so that the arraylist
		// of dataframes
		// can be assembled again in order even though each one is read out of order
		List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
		Map<Integer, String> doneReading2 = new ConcurrentHashMap<>();
		for (int i = 0; i < sampleNames.size(); i++) {
			toRead.add(i);
		}
		int N_THREADS = Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		int totalNFiles = toRead.size();
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {
					while (toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading2.put(currentFile, String.join("\t", countsDF.getColumn("COUNT")));
							progressPercentage(totalNFiles - toRead.size(), totalNFiles,
									barcodeCountsFiles.get(currentFile));
							countsDF = null; // java gc is a fickle mistress
						} catch (IOException e) {
							System.out.println();
							e.printStackTrace();
							System.out.println("error with file: " + currentFile);
						}
					}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

		ArrayList<String> countsArrayList = new ArrayList<>();
		for (int i = 0; i < doneReading2.size(); i++) {
			countsArrayList.add(doneReading2.get(i));
		}

		System.out.println("done reading files");

		// 'labels' is the exon labels, / column names
		System.out.print("generating counts matrix. \t");
		ArrayList<String> labels = new ArrayList<>();
		DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
		for (int i = 0; i < countsDF.nrow(); i++) {
			labels.add(countsDF.get("CONTIG", i) + "_" + countsDF.get("START", i) + "_" + countsDF.get("END", i));
		}
		System.out.println("done generating counts matrix");

		// for sanity checking, number should all appropriately match
		System.out.println("sampleNames.size()\t" + sampleNames.size());
		System.out.println("countsArrayList.size()\t" + countsArrayList.size());
		System.out.println("labels.size()\t" + labels.size());
		System.out.println("writing counts matrix");

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		for (int i = 0; i < labels.size(); i++) {
			output.write(labels.get(i));
			if (i != labels.size() - 1) {
				output.write("\t");
			}
		}
		output.write("\n");

		if (countsArrayList.size() != sampleNames.size()) {
			System.out.println("error 51");
		}

		// written in such a way as to be readable in R, adding rownames for the sample
		// ID
		for (int i = 0; i < sampleNames.size(); i++) {
			StringBuilder line = new StringBuilder();
			line.append(sampleNames.get(i) + "\t");
			line.append(countsArrayList.get(i));
			if (i != sampleNames.size() - 1) {
				line.append("\n");
			}
			output.write(line.toString());
		}
		output.close();
	}

	
}
