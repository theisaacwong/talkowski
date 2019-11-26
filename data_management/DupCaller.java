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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * 
 * @author Isaac Wong
 * November 26, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski Lab
 * 
 * This is designed to read in a list of files, with a full file path, and then calculate md5sum on all files larger than some value
 *
 */


public class DupCaller {

	public final static int N_THREADS = Runtime.getRuntime().availableProcessors();
	public Map<Long, ArrayList<String>> sizeToPath;
	public Map<String, Long> pathToSize;
	public Map<String, ArrayList<String>> md5ToPath;
	public Map<Long, HashSet<String>> sizeToMd5;
	public List<String> filePaths;
	
	public DupCaller() {
		md5ToPath = new ConcurrentHashMap<>();
		sizeToMd5 = new ConcurrentHashMap<>();
		sizeToPath = new ConcurrentHashMap<>();
		pathToSize = new ConcurrentHashMap<>();
		filePaths = Collections.synchronizedList(new ArrayList<String>()); 
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		System.out.println("Duplicate File Caller version 0.14.8");
		System.out.println("Java version: " + System.getProperty("java.version"));
		System.out.println("N_THREADS: " + N_THREADS);
		if(args[0].equals("-h") || args[0].equals("--help")) {
			System.out.println("java -jar DupCaller [input_file] [output_file]");
		} else {
			String INPUT_PATH = args[0];
			String OUTPUT_PATH = args[1];
			
			DupCaller dupCaller = new DupCaller();
			dupCaller.run(INPUT_PATH, OUTPUT_PATH);
		}

	}
	
	public void run(String INPUT_PATH, String OUTPUT_PATH) throws IOException, InterruptedException {
		long t1 = System.currentTimeMillis();
		System.out.println("step\tlap_time(s)\telapsed_time(s)\theap_size");
		
		System.out.println("Reading in file\t0\t0\t" + getHeapSize());
		this.parseFile(INPUT_PATH);
		long t2 = System.currentTimeMillis();
		System.out.println("Done reading file\t" + (t2 - t1)/1000F + "\t" + (t2 - t1)/1000F + "\t" + getHeapSize());

		System.out.println("Calculating md5 sums");
		this.calcMD5();
		long t3 = System.currentTimeMillis();
		System.out.println("Done calculating\t" + (t3 - t2)/1000F + "\t" + (t3 - t1)/1000F + "\t" + getHeapSize());
		
		System.out.println("writing md5 sums");
		this.writeFile(OUTPUT_PATH);
		long t4 = System.currentTimeMillis();
		System.out.println("Done calculating\t" + (t4 - t3)/1000F + "\t" + (t4 - t1)/1000F + "\t" + getHeapSize());
		
	}
	
	public void writeFile(String OUTPUT_PATH) throws IOException {
		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));
		output.write("file_path\tasize\tmd5\n");
		
		ArrayList<Long> sortIndex = new ArrayList<>(sizeToMd5.keySet());
		Collections.sort(sortIndex, Collections.reverseOrder());
		
		//System.out.println(toString(md5ToPath));
		
		//System.out.println("sortIndex.toString() "+ sortIndex.toString());
		
		for(Long l : sortIndex) {
			String hSize = humanReadableByteCount(l);			
			for(String md5 : sizeToMd5.get(l)) {			
				if(md5ToPath.get(md5).size() > 1) {
					for(String filePath : md5ToPath.get(md5)) {
						output.write(filePath + "\t" + hSize + "\t" + md5 + "\n");
					}					
				}
			}
		}
		output.close();
	}
	
	public void calcMD5() throws IOException, InterruptedException {
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {

					while(filePaths.size() > 0) {
						String currentFile = filePaths.remove(0);
						//System.out.println(Thread.currentThread().getId() + "\t" + currentFile);
						
						String md5val = "na";
						try {
							Runtime r = Runtime.getRuntime();
							Process p = r.exec("md5sum " + currentFile);
							//Process p = r.exec("bash -c 'md5sum " + currentFile + "'");
							p.waitFor();
							BufferedReader b = new BufferedReader(new InputStreamReader(p.getInputStream()));
							md5val = b.readLine().split(" ")[0];
							b.close();
						} catch (Exception e) {
							e.printStackTrace();
						}
						
						if(!md5ToPath.containsKey(md5val))
							md5ToPath.put(md5val, new ArrayList<String>());
						md5ToPath.get(md5val).add(currentFile);
						
						Long asize = pathToSize.get(currentFile);
						if(!sizeToMd5.containsKey(asize))
							sizeToMd5.put(asize, new HashSet<String>());
						sizeToMd5.get(asize).add(md5val);
					}
				}
			});
		}
		// Initiates an orderly shutdown in which previously submitted tasks are executed, but no new tasks will be accepted. Invocation
		// has no additional effect if already shut down.
		// This method does not wait for previously submitted tasks to complete execution. Use awaitTermination to do that.
		exServer.shutdown();
		
		// Blocks until all tasks have completed execution after a shutdown request, or the timeout occurs, or the current thread is
		// interrupted, whichever happens first.
		exServer.awaitTermination(6, TimeUnit.DAYS);
	}
	
	public void parseFile(String INPUT_PATH) throws FileNotFoundException {
		FileInputStream inputStream = null;
		Scanner sc = null;
		inputStream = new FileInputStream(INPUT_PATH);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		sc.nextLine();
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			String filePath = line.split("\t")[0];
			Long asize = Long.parseLong(line.split("\t")[1]);

			filePaths.add(filePath);
			
			pathToSize.put(filePath, asize);
			
			if(!sizeToPath.containsKey(asize)) 
				sizeToPath.put(asize, new ArrayList<>());
			sizeToPath.get(asize).add(filePath);
		}
		sc.close();
	}
	
	public static String getHeapSize() {
		return humanReadableByteCount(Runtime.getRuntime().totalMemory());
	}
	
	public static String humanReadableByteCount(long bytes) {
		if (bytes < 1024) return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp-1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}
	
	public String toString(Map<String, ArrayList<String>> m) {
		StringBuilder s = new StringBuilder();
		s.append(m.size());
		s.append("\n");
		for(String k : m.keySet() ) {
			s.append(k);
			s.append("\n");
			for(String al : m.get(k)) {
				s.append("\t");
				s.append(al);
				s.append("\n");
			}
		}
		return s.toString();
	}
	
	
	
	
	
	
}
