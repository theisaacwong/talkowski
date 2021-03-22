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
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

public class GetCountsMatrixBuffered extends gCNVHelperTool {

	public GetCountsMatrixBuffered(String[] args) {
		super(args);
	}

	@Override
	public void run() throws IOException, InterruptedException {
		String sourceFolder = args[1];
		String OUTPUT_PATH = args[2];
		String countsRegex = args[3];
		int BUFFER_SIZE = Integer.parseInt(args[4]);
		this.getCountsMatrixBuffered(sourceFolder, OUTPUT_PATH, countsRegex, BUFFER_SIZE);
	}
	

	/**
	 * write buffer after every n files is read per thread, uses less memory as less
	 * count data is held in memory at any given time however, uses more I/O, which
	 * will be slower overall
	 * 
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @param countsRegex
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrixBuffered(String sourceFolder, String OUTPUT_PATH, String countsRegex, int BUFFER_SIZE)
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
		List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
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
					File file = new File(OUTPUT_PATH + "_" + Thread.currentThread().getId());
					BufferedWriter output = null;
					try {
						output = new BufferedWriter(new FileWriter(file));
					} catch (IOException e1) {
						e1.printStackTrace();
					}

					HashMap<String, String> doneReading = new HashMap<>(); // map of sample name to line string

					int bufferCounter = 0;

					ArrayList<String> labels = new ArrayList<>();
					DataFrame labelsDF = null;
					try {
						labelsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
					} catch (IOException e1) {
						e1.printStackTrace();
					}
					for (int i = 0; i < labelsDF.nrow(); i++) {
						labels.add(labelsDF.get("CONTIG", i) + "_" + labelsDF.get("START", i) + "_"
								+ labelsDF.get("END", i));
					}
					labelsDF = null;

					for (int i = 0; i < labels.size(); i++) {
						try {
							output.write(labels.get(i));
						} catch (IOException e) {
							e.printStackTrace();
						}
						if (i != labels.size() - 1) {
							try {
								output.write("\t");
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
					}
					try {
						output.write("\n");
					} catch (IOException e1) {
						e1.printStackTrace();
					}

					while (toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading.put(sampleNames.get(currentFile),
									String.join("\t", countsDF.getColumn("COUNT")));
							progressPercentage(totalNFiles - toRead.size(), totalNFiles,
									barcodeCountsFiles.get(currentFile));
							countsDF = null; // java gc is a fickle mistress

							if (bufferCounter >= BUFFER_SIZE) {
								StringBuilder line = new StringBuilder();
								for (String snKey : doneReading.keySet()) {
									line.append(snKey); // sample name key
									line.append("\t");
									line.append(doneReading.get(snKey));
									line.append("\n");
								}
								output.write(line.toString());
								doneReading = new HashMap<>();
								bufferCounter = 0;
							} else {
								bufferCounter++;
							}

						} catch (IOException e) {
							e.printStackTrace();
						}
					}

					// flush remaining buffer
					StringBuilder line = new StringBuilder();
					for (String snKey : doneReading.keySet()) {
						line.append(snKey); // sample name key
						line.append("\t");
						line.append(doneReading.get(snKey));
						line.append("\n");
					}
					try {
						output.write(line.toString());
					} catch (IOException e2) {
						e2.printStackTrace();
					}
					doneReading = new HashMap<>();
					bufferCounter = 0;

					try {
						output.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

	}
	
}
