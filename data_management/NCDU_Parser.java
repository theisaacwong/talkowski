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
import java.io.FileNotFoundException;
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

	public static HashMap<String, String> errorMessages = new HashMap<String, String>();
	public static final String HEADER = "file\tasize\tdsize\tuser\thuman_readable\tfile_extension\n";

	public static void main(String[] args) throws IOException {
		//args = new String[] {"C:/Users/iwong/Documents/MGH/temp/ncdu.txt", "C:/Users/iwong/Documents/MGH/temp/ncdu_"};
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

	public static String getFileName(String str) {
		String file = str.split("\t")[0];
		String[] filee = file.split("/");
		int len = filee.length;
		return filee[len - 1];
	}
	
	public static HashMap<String, String> getUIDtoNameMap(String pathToPasswd) throws FileNotFoundException{
		HashMap<String, String> uidToName = new HashMap<>();
		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToPasswd);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			String[] linee = line.split(":");
			uidToName.put(linee[2], linee[0]);
		}
		sc.close();
		return uidToName;
	}

	public static void ncduParser(String INPUT_PATH, String OUTPUT_BASE) throws IOException {
		long start = System.currentTimeMillis();
		System.out.println("Reading file: ");

		BufferedWriter output = null;
		File file = new File(OUTPUT_BASE + ".tsv");
		output = new BufferedWriter(new FileWriter(file));

		BufferedWriter err = null;
		File errfile = new File("err.txt");
		err = new BufferedWriter(new FileWriter(errfile));

		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(INPUT_PATH);
		sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		sc.nextLine();

		Stack<String> filePaths = new Stack<String>();

		output.write(HEADER);

		HashMap<Long, String> lineMap = new HashMap<Long, String>();
		HashMap<Long, ArrayList<Long>> sortMap = new HashMap<Long, ArrayList<Long>>();
		long lineNumber = 1;

		Pattern ap = Pattern.compile("(?<=asize\":)\\d+");
		Pattern dp = Pattern.compile("(?<=dsize\":)\\d+");
		Pattern np = Pattern.compile("(?<=name\":\").+?(?=\")");
		Pattern ep = Pattern.compile("(?<=\\.).*");
		Pattern up = Pattern.compile("(?<=uid\":)\\d+");
		
		int temp_1 = -1;

		var uidToName = getUIDtoNameMap("/etc/passwd");
		
		/*
		 * For each line
		 * 		if line starts with '['
		 * 			push "name" to the filePath stack
		 * 			while line has ending ']'
		 * 				remove the ending ']'
		 * 				pop top dir name from stack
		 * 		else
		 * 			add "name" to the lineEntry to write
		 * 		 	while line has ending ']'
		 * 				remove the ending ']'
		 * 				pop top dir name from stack
		 */
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.charAt(0) == '[') {
				Matcher nm = np.matcher(line);  
				if(nm.find()) {
					filePaths.push(nm.group());
				} else { // error
					System.out.println("error 38472");
				}
				while(line.lastIndexOf("]") > line.lastIndexOf("}")) {
					temp_1 = line.lastIndexOf("]");
					line = line.substring(0, temp_1) + line.substring(temp_1 + 1);
					if(!filePaths.isEmpty()) {
						filePaths.pop();	
					}
				}				
			} else {
				Matcher am = ap.matcher(line); String asize = am.find() ? am.group() : "-1";
				Matcher dm = dp.matcher(line); String dsize = dm.find() ? dm.group() : "-1";
				Matcher nm = np.matcher(line); String name  = nm.find() ? nm.group() : "NA";
				Matcher em = ep.matcher(name); String fileExtension = em.find() ? em.group() : "NA";
				Matcher um = up.matcher(line); String uid = um.find() ? um.group() : "NA";
				
				uid = uidToName.containsKey(uid) ? uidToName.get(uid) : uid;
				
				
				StringBuilder sb = new StringBuilder();
				for (Object s : filePaths.toArray()) {
					sb.append(s);
					sb.append("/");
				}
				sb.append(name); 	sb.append("\t");
				sb.append(asize); 	sb.append("\t");
				sb.append(dsize);	sb.append("\t");
				sb.append(uid);		sb.append("\t");
				sb.append(humanReadableByteCount(asize));	sb.append("\t");
				sb.append(fileExtension);	sb.append("\n");
				String currPath = sb.toString();
				output.write(currPath);

				if(currPath.contains("talkowski") == false) {
					err.write(lineNumber + "\t" + line + "\t" + currPath + "\n");
					System.out.println(lineNumber + "\t" + line + "\t" + currPath + "\n");
				}

				while(line.lastIndexOf("]") > line.lastIndexOf("}")) {
					temp_1 = line.lastIndexOf("]");
					line = line.substring(0, temp_1) + line.substring(temp_1 + 1);
					if(!filePaths.isEmpty()) {
						filePaths.pop();	
					}
				}

				Long aSize = Long.parseLong(asize);
				lineMap.put(lineNumber, currPath); 
				if(sortMap.containsKey(aSize)) { 
					sortMap.get(aSize).add(lineNumber); 
				} else {
					sortMap.put(aSize, new ArrayList<Long>());
					sortMap.get(aSize).add(lineNumber); 
				}
				lineNumber++;					
			}
		}
		inputStream.close();
		sc.close();
		output.close();
		err.close();

		long doneReading = System.currentTimeMillis();
		System.out.println("done reading: ");
		System.out.println("elapsed time: " + (doneReading - start) / 1000F + " seconds");
		System.out.println("now sorting, n = " + sortMap.keySet().size());

		output = null;
		file = new File(OUTPUT_BASE + "_sorted.tsv");
		output = new BufferedWriter(new FileWriter(file));

		output.write(HEADER);
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


		output = null;
		file = new File(OUTPUT_BASE + "_sorted_duplicates.tsv");
		output = new BufferedWriter(new FileWriter(file));
		output.write(HEADER);

		// for each file size bucket/bin observed
		for(Long si : sortIndex) {
			// if that bucket is greater than two, there might be a duplicate	- this check is unnecessary, but saves time
			if(sortMap.get(si).size() > 1) {
				// create a hashmap where keys:baseFileName and values:lineEntryToWrite
				var fileNames = new HashMap<String, ArrayList<String>>();
				// for each line in the current file bin size, add the file path to the hashmap
				for(Long ln : sortMap.get(si)) {
					String fn = getFileName(lineMap.get(ln));
					// basic, checks to see if file is already observed or not
					if(fileNames.containsKey(fn)) {
						fileNames.get(fn).add(lineMap.get(ln));
					} else {
						fileNames.put(fn, new ArrayList<String>());
						fileNames.get(fn).add(lineMap.get(ln));
					}
				}
				// for all value arraylists, write their contents if size>=2
				for(var fileDups : fileNames.values()) {
					if(fileDups.size() > 1) {
						for(String lineEntry : fileDups) {
							output.write(lineEntry);
						}
					}
				}
			}
		}
		output.close();
		System.out.println("done writing");

	}

}

