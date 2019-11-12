package data_management;
/**
 * Isaac Wong
 * November 5, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski Lab
 * 
 * This is designed to parse the output of ncdu (https://dev.yorhel.nl/ncdu) and list files with their full files paths and file size
 * Also output will be a sorted file so that large files can easily be identified (and deleted/archived)
 * 
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Stack;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NCDU_Parser {

	public static void main(String[] args) throws IOException {
		System.out.println("Java version " + System.getProperty("java.version"));
		if(args[0].equals("-h") || args[0].equals("--help")) {
			System.out.println("java NCDU_Parser [ncdu_file] [output_file]");
		} else {
			ncduParser(args[0], args[1]);
		}
	}


	/*
	 * https://programming.guide/java/formatting-byte-size-to-human-readable-format.html
	 * https://stackoverflow.com/questions/3758606/how-to-convert-byte-size-into-human-readable-format-in-java
	 */
	public static String humanReadableByteCount(long bytes) {
	    if (bytes < 1024) return bytes + " B";
	    int exp = (int) (Math.log(bytes) / 6.907755278982137);
	    char pre = "KMGTPE".charAt(exp-1);
	    return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}
	public static String humanReadableByteCount(String bytes_string) {
		long bytes = Long.parseLong(bytes_string);
	    if (bytes < 1024) return bytes + " B";
	    int exp = (int) (Math.log(bytes) / 6.907755278982137);
	    char pre = "KMGTPE".charAt(exp-1);
	    return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}
	
	
	
	public static void ncduParser(String INPUT_PATH, String OUTPUT_BASE) throws IOException {
		long start = System.currentTimeMillis();
		System.out.println("Reading file: ");

		BufferedWriter output = null;
		File file = new File(OUTPUT_BASE + ".tsv");
		output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(INPUT_PATH);
		sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		sc.nextLine();

		Stack<String> filePaths = new Stack<String>();

		output.write("file\tasize\tdsize\thuman_readable\tfile_extension\n");

		HashMap<Long, String> lineMap = new HashMap<Long, String>();
		HashMap<Long, ArrayList<Long>> sortMap = new HashMap<Long, ArrayList<Long>>();
		long lineNumber = 1;

		Pattern ap = Pattern.compile("(?<=asize\":)\\d+");
		Pattern dp = Pattern.compile("(?<=dsize\":)\\d+");
		Pattern np = Pattern.compile("(?<=name\":\").+?(?=\")");
		Pattern ep = Pattern.compile("(?<=\\.).*");
		int count = 0;

		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			Matcher am = ap.matcher(line); 
			Matcher dm = dp.matcher(line); 
			Matcher nm = np.matcher(line);  nm.find();

			if(line.contains("[")) {	// new folder opening
				filePaths.push(nm.group());
				count = line.length() - line.replaceAll("]", "").length();	// folder(s) are closed
				for(int i = 0; i < count; i++) {
					if(filePaths.size() > 0)
						filePaths.pop();
				}
				continue;
			} 
			
			if(am.find()==false || dm.find()==false) {
				count = line.length() - line.replaceAll("]", "").length();
				for(int i = 0; i < count; i++) {
					if(filePaths.size() > 0)
						filePaths.pop();
				}
				continue;
			}

			String asize = am.group(0);
			String dsize = dm.group(0);
			String name  = nm.group(0);
			String fileExtension = "NA";
			
			Matcher em = ep.matcher(name);
			if(em.find()) {
				fileExtension = em.group(0);
			} 
			
			StringBuilder sb = new StringBuilder();
			for (Object s : filePaths.toArray()) {
			    sb.append(s);
			    sb.append("/");
			}
			sb.append(name);
			sb.append("\t");
			sb.append(asize);
			sb.append("\t");
			sb.append(dsize);
			sb.append("\t");
			sb.append(humanReadableByteCount(asize));
			sb.append("\t");
			sb.append(fileExtension);
			sb.append("\n");
			String currPath = sb.toString();
			output.write(currPath);
			

			count = line.length() - line.replaceAll("]", "").length();
			for(int i = 0; i < count; i++) {
				if(filePaths.size() > 0)
					filePaths.pop();
			}

			lineMap.put(lineNumber, currPath);
			Long aSize = Long.parseLong(asize);
			if(sortMap.containsKey(aSize)) {
				sortMap.get(aSize).add(lineNumber);
			} else {
				sortMap.put(aSize, new ArrayList<Long>());
				sortMap.get(aSize).add(lineNumber);
			}
			lineNumber++;
		}
		inputStream.close();
		sc.close();
		output.close();

		long doneReading = System.currentTimeMillis();
		System.out.println("done reading: ");
		System.out.println("elapsed time: " + (doneReading - start) / 1000F + " seconds");
		System.out.println("now sorting, n = " + sortMap.keySet().size());

		output = null;
		file = new File(OUTPUT_BASE + "_sorted.tsv");
		output = new BufferedWriter(new FileWriter(file));

		output.write("file\tasize\tdsize\thuman_readable\n");
		ArrayList<Long> sortIndex = new ArrayList<>(sortMap.keySet());
		Collections.sort(sortIndex, Collections.reverseOrder());

		long doneSorting = System.currentTimeMillis();
		System.out.println("elapsed time: " + (doneSorting - doneReading) / 1000F + " seconds");
		System.out.println("done sorting, now writing");

		for(Long si : sortIndex) {
			for(Long ln : sortMap.get(si)) {
				output.write(lineMap.get(ln));
			}
		}
		output.close();
		System.out.println("done writing");
		long end = System.currentTimeMillis();
		System.out.println("total elapsed time: " + (end - start) / 1000F + " seconds");

	}

}

