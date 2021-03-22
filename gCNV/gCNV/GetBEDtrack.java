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

public class GetBEDtrack extends gCNVHelperTool {

	public GetBEDtrack(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException, InterruptedException {
		String sourceFolder = args[1];
		String OUTPUT_PATH = args[2];
		String countsRegex = args[3];
		this.getBEDtrack(sourceFolder, OUTPUT_PATH, countsRegex);
	}
	

	public void getBEDtrack(String sourceFolder, String OUTPUT_PATH, String countsRegex) throws IOException, InterruptedException {
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
		Map<Integer, ArrayList<String>> doneReading2 = new ConcurrentHashMap<>();
		for (int i = 0; i < sampleNames.size(); i++) {
			toRead.add(i);
		}
		int N_THREADS = Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		int totalNFiles = toRead.size();
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@SuppressWarnings("unchecked")
				@Override
				public void run() {
					while (toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							//doneReading2.put(currentFile, countsDF.getColumn(countsDF.fieldNames[3]));
							doneReading2.put(currentFile, (ArrayList<String>) countsDF.getColumn(countsDF.fieldNames[3]).clone()); // Shallow copy should work right
							progressPercentage(totalNFiles - toRead.size(), totalNFiles,
									barcodeCountsFiles.get(currentFile));
							countsDF = null; // java gc is a fickle mistress
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

		//		ArrayList<String> countsArrayList = new ArrayList<>();
		//		for (int i = 0; i < doneReading2.size(); i++) {
		//			countsArrayList.add(doneReading2.get(i));
		//		}

		System.out.println("done reading files");

		// 'labels' is the exon labels, / column names
		System.out.print("generating counts matrix. \t");
		ArrayList<String> chr = new ArrayList<>();
		ArrayList<String> start = new ArrayList<>();
		ArrayList<String> end = new ArrayList<>();
		DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
		for (int i = 0; i < countsDF.nrow(); i++) {
			chr.add(countsDF.get(i, 0));
			start.add(countsDF.get(i, 1));
			end.add(countsDF.get(i, 2));
		}
		System.out.println("done generating counts matrix");

		// for sanity checking, number should all appropriately match
		System.out.println("sampleNames.size()\t" + sampleNames.size());
		System.out.println("countsArrayList.size()\t" + doneReading2.get(0).size());
		System.out.println("labels.size()\t" + chr.size());
		System.out.println("writing counts matrix");

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		// write the header
		StringBuilder header = new StringBuilder();
		header.append("#chr\tstart\tend\t");
		for(int i = 0; i < sampleNames.size()-1; i++) {
			header.append(sampleNames.get(i));
			header.append("\t");
		}
		header.append(sampleNames.get(sampleNames.size()-1));
		header.append("\n");
		output.write(header.toString());


		for(int r = 0; r < chr.size(); r++) {
			StringBuilder line = new StringBuilder();
			line.append(chr.get(r));
			line.append("\t");
			line.append(start.get(r));
			line.append("\t");
			line.append(end.get(r));
			line.append("\t");
			for(int c = 0; c < sampleNames.size()-1; c++) {
				line.append(doneReading2.get(c).get(r));
				line.append("\t");
			}
			line.append(doneReading2.get(sampleNames.size()-1).get(r));
			line.append("\n");
			output.write(line.toString());
		}

		output.close();
	}

}
