package data_management;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class DupCallerPipe {
	public final static int N_THREADS = Runtime.getRuntime().availableProcessors();
	
	public static void main(String[] args) throws NumberFormatException, IOException, InterruptedException {
		if(args.length == 4 && args[0].equals("splitFiles")) {
			splitFiles(args[1], Integer.parseInt(args[2]), args[3]);
		} else if(args.length == 2 && args[0].equals("calcMD5")) {
			calcMD5(args[1]);
		} else if(args.length == 4 && args[0].equals("mergeFiles")){
			mergeFiles(args[1], Integer.parseInt(args[2]), args[3]);
		} else {
			System.out.println("Usage:");
			System.out.println("java -jar DupCallerPipe splitFiles [INPUT_PATH] [nFiles] [BASE_PATH]");
			System.out.println("java -jar DupCallerPipe calcMD5 [INPUT_PATH]");
			System.out.println("java -jar DupCallerPipe mergeFiles [BASE_PATH] [nFiles] [OUTPUT_PATH]");
		}
	}
	
	public static void sortFiles(String OUTPUT_PATH, String INPUT_PATH) {
		
		
		
	}
	
	public static void mergeFiles(String BASE_PATH, int n, String OUTPUT_PATH) throws IOException {
		ArrayList<ArrayList<String>> files = new ArrayList<>();
		HashMap<Long, HashSet<String>> sizeToMD5 = new HashMap<>();
		HashMap<String, ArrayList<Integer[]>> MD5ToTupleFileAndLineNumber = new HashMap<>(); //key: asize, value: al of tuples of file index and line number
		
		for(int i = 0; i < n; i++) {
			files.add(new ArrayList<>());
			FileInputStream inputStream = new FileInputStream(BASE_PATH + "_" + (i+1));;
			Scanner sc = new Scanner(inputStream, "UTF-8");
			System.out.println("Reading File: " + BASE_PATH + "_" + (i+1));
			String line = "";
			int lineNumber = 0;
			while (sc.hasNextLine()) {
				try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
				
				String[] linee = line.split("\t");
				if(linee.length < 2) {
					continue;
				}
				
				files.get(i).add(line);
				
				Long asize = Long.parseLong(linee[1]);
				String md5 = linee[5];
				
				if(!sizeToMD5.containsKey(asize)) {
					sizeToMD5.put(asize, new HashSet<>());
				}
				sizeToMD5.get(asize).add(md5);
				
				if(!MD5ToTupleFileAndLineNumber.containsKey(md5)) {
					MD5ToTupleFileAndLineNumber.put(md5, new ArrayList<>());
				}
				MD5ToTupleFileAndLineNumber.get(md5).add(new Integer[] {i, lineNumber});
				
				lineNumber++;
			}
			sc.close();			
		}
		
		ArrayList<Long> sortedSizes = new ArrayList<>(sizeToMD5.keySet());
		Collections.sort(sortedSizes, Collections.reverseOrder());
		
		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));
		output.write("file\tasize\tdsize\thuman_readable\tfile_extension\tmd5sum\tnotes\n");		
		
		for(int i = 0; i < sortedSizes.size(); i++) {
			Long currSize = sortedSizes.get(i);
			for(String currMD5 : sizeToMD5.get(currSize)) {
				if(MD5ToTupleFileAndLineNumber.get(currMD5).size() < 1) {
					continue;
				}
				for(int k = 0; k < MD5ToTupleFileAndLineNumber.get(currMD5).size(); k++) {
					int fileNumber = MD5ToTupleFileAndLineNumber.get(currMD5).get(k)[0];
					int lineNumber = MD5ToTupleFileAndLineNumber.get(currMD5).get(k)[1];
					output.write(files.get(fileNumber).get(lineNumber) + "\n");
				}
			}
		}
		output.close();
	}
	
	
	public static void splitFiles(String INPUT_PATH, int nFiles, String BASE_PATH) throws IOException {
		ArrayList<BufferedWriter> outputs = new ArrayList<>();
		for(int i = 1; i <= nFiles; i++) {
			outputs.add(new BufferedWriter(new FileWriter(new File(BASE_PATH + "_" + i))));
		}
		int currWriter = 0;
		
		FileInputStream inputStream = new FileInputStream(INPUT_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		sc.nextLine(); // skip the header line
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			outputs.get(currWriter).write(line + "\n");
			currWriter = (currWriter + 1) % nFiles;
		}
		sc.close();
		for(int i = 0; i < outputs.size(); i++) {
			outputs.get(i).close();
		}
	}
	
	
	public static void calcMD5(String INPUT_PATH) throws FileNotFoundException, InterruptedException {
		List<String> lines = Collections.synchronizedList(new ArrayList<String>());
		
		FileInputStream inputStream = null;
		Scanner sc = null;
		inputStream = new FileInputStream(INPUT_PATH);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			lines.add(line);
		}
		sc.close();
		
		System.out.println("n threads: "  + N_THREADS + "\t" + INPUT_PATH);
		
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {
					System.out.println("Thread: " + Thread.currentThread().getId());
					while(lines.size() > 0) {
						String currentLine = lines.remove(0);
						String currentFile = currentLine.split("\t")[0];
						
						String md5val = "na";
						try {
							Runtime r = Runtime.getRuntime();
							Process p = r.exec("md5sum " + currentFile); // linux
							p.waitFor();
							BufferedReader b = new BufferedReader(new InputStreamReader(p.getInputStream()));
							md5val = b.readLine().split(" ")[0];
							b.close();
						} catch (Exception e) {
							e.printStackTrace();
						}
						
						System.out.println(currentLine + "\t" + md5val);
					}
				}
			});
		}
		System.out.println("heap size: " + getHeapSize());
		// Initiates an orderly shutdown in which previously submitted tasks are executed, but no new tasks will be accepted. Invocation
		// has no additional effect if already shut down.
		// This method does not wait for previously submitted tasks to complete execution. Use awaitTermination to do that.
		exServer.shutdown();
		
		// Blocks until all tasks have completed execution after a shutdown request, or the timeout occurs, or the current thread is
		// interrupted, whichever happens first.
		exServer.awaitTermination(1, TimeUnit.DAYS);
	}
	
	/**
	 * 
	 * @return returns the current JVM heap size in human readable format
	 */
	public static String getHeapSize() {
		return humanReadableByteCount(Runtime.getRuntime().totalMemory());
	}
	
	/**
	 * 
	 * @param bytes
	 * @return human readable version of bytes
	 */
	public static String humanReadableByteCount(long bytes) {
		if (bytes < 1024) return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp-1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}
	
	
	
	
}
