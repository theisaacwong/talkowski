package data_management;


/**
 * @author Isaac Wong <iwong1@mgh.harvard.edu>
 * @version 2.1 - January 6, 2020
 * @since November 5, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski Lab
 * 
 * Data Management Toolkit
 * This is a small set of tools to manage large directories. 
 * Designed for use in conjunction with ncdu (https://dev.yorhel.nl/ncdu) 
 * 
 */

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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Stack;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DATA_MANAGEMENT {

	public final static int N_THREADS = Runtime.getRuntime().availableProcessors();
	public static final String ARCHIVE = "archive";
	public static final String IMPORTANT = "important";
	public static final String USER_FIELD_NAME = "user";
	public static final String ASIZE_FIELD_NAME = "asize";
	public static final String HEADER = "file\tasize\tdsize\tuser\thuman_readable\tfile_extension\tnotes\n"; // generic header
	public static final String MD5_HEADER = "file\tasize\tdsize\tuser\thuman_readable\tfile_extensiont\tmd5sum\tnotes\n"; // header of files with MD5 info

	public static void main(String[] args) throws NumberFormatException, IOException, InterruptedException {
		if(args.length == 4 && args[0].equalsIgnoreCase("splitFiles")) {
			splitFiles(args[1], Integer.parseInt(args[2]), args[3]);
		} else if(args.length == 2 && args[0].equalsIgnoreCase("calcMD5")) {
			calcMD5(args[1]);
		} else if(args.length == 4 && args[0].equalsIgnoreCase("mergeFiles")){
			mergeFiles(args[1], Integer.parseInt(args[2]), args[3]);
		} else if(args.length == 3 && args[0].equalsIgnoreCase("addOwner")){
			addOwner(args[1], args[2]);
		} else if(args.length == 4 && args[0].equalsIgnoreCase("shame")) {
			shame(args[0], args[1], args[2], args[3]);
		} else if(args.length == 3 && args[0].equalsIgnoreCase("ncduParser")){
			ncduParser(args[1], args[2]);
		}else {
			System.out.println("Usage:");
			System.out.println("java -jar DMTK.jar splitFiles [INPUT_PATH] [nFiles] [BASE_PATH]");
			System.out.println("java -jar DMTK.jar calcMD5 [INPUT_PATH]");
			System.out.println("java -jar DMTK.jar mergeFiles [BASE_PATH] [nFiles] [OUTPUT_PATH]");
			System.out.println("java -jar DMTK.jar shame [BLACKLIST] [USER_DIRECTORIES] [NCDU_PARSE_FILE] [OUTPUT_PATH]");
			System.out.println("java -jar DMTK.jar ncduParser [ncdu_file] [output_file]");
		}
	}


	/***************************************************************************************
	 ***************************************************************************************
	 ************************** DUPLICATE CALLER PIPELINE **********************************
	 ***************************************************************************************
	 ***************************************************************************************
	 */


	/**
	 * This is designed to read in multiple files, sort them, and write a single file to output
	 * @since v1.0
	 * @param BASE_PATH - The base file name string, files will be names basename_1, basename_2, basename_3, etc
	 * @param n - The number of files to merge
	 * @param OUTPUT_PATH - the output path to write to 
	 * @throws IOException
	 */
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
		//output.write("file\tasize\tdsize\thuman_readable\tfile_extension\tmd5sum\tnotes\n");		
		output.write(MD5_HEADER);

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


	/**
	 * Reads in a file and splits it into n lines, line 1 is written to file 1, line 2 to file 2, line 3 to file 3, etc. 
	 * line x is written to file x%n
	 * @since v1.0
	 * @param INPUT_PATH - file to split
	 * @param nFiles - number of files to split into 
	 * @param BASE_PATH - base name of shard files, eg basename_1, basename_2, basename_3
	 * @throws IOException
	 */
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

	/**
	 * @deprecated, too slow, I/O is awful
	 * @param INPUT_PATH
	 * @param OUTPUT_PATH
	 * @throws IOException
	 */
	public static void addOwner(String INPUT_PATH, String OUTPUT_PATH) throws IOException {

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = null;
		Scanner sc = null;
		inputStream = new FileInputStream(INPUT_PATH);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			String currentFile = line.split("\t")[0];
			String owner = "NaN";

			try {
				Runtime r = Runtime.getRuntime();
				String[] cmd = new String[] {"ls", "-l", currentFile};
				Process p = r.exec(cmd); // linux
				p.waitFor();
				BufferedReader b = new BufferedReader(new InputStreamReader(p.getInputStream()));
				owner = b.readLine().split(" ")[2];
				b.close();
			} catch (Exception e) {
				System.out.println(line);
				e.printStackTrace();
			}


			output.write(line + "\t" + owner + "\n");

		}
		sc.close();
		output.close();

	}

	/**
	 * Reads in a file where the first column is the full file path, 
	 * Writes the md5sum value of the file
	 * @since v1.0
	 * @param INPUT_PATH - path to a shard file
	 * @throws FileNotFoundException
	 * @throws InterruptedException
	 */
	public static void calcMD5(String INPUT_PATH) throws FileNotFoundException, InterruptedException {
		List<String> lines = Collections.synchronizedList(new ArrayList<String>());

		FileInputStream inputStream = new FileInputStream(INPUT_PATH);;
		Scanner sc = new Scanner(inputStream, "UTF-8");;
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

						String md5val = "NaN";
						try {
							Runtime r = Runtime.getRuntime();
							Process p = r.exec("md5sum " + currentFile); // linux
							p.waitFor();
							BufferedReader b = new BufferedReader(new InputStreamReader(p.getInputStream()));
							md5val = b.readLine().split(" ")[0];
							b.close();
						} catch (Exception e) {
							System.out.println("Eror: " + currentLine);
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






	/***************************************************************************************
	 ***************************************************************************************
	 *********************************** SHAME FILE*****************************************
	 ***************************************************************************************
	 ***************************************************************************************
	 */


	/**
	 * Returns a map where
	 * 		key: file, full path
	 * 		value: classification, eg who is the owner, isit archive, important, etc
	 * @since v1.0
	 * @param pathToBlacklist - full path to blacklist
	 * @return - a mapping of filepath to user classification
	 * @throws FileNotFoundException
	 */
	public static HashMap<String, String> getBlacklist(String pathToBlacklist) throws FileNotFoundException {
		HashMap<String, String> fileToClassification = new HashMap<>();
		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToBlacklist);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			String[] linee = line.split("\\t");
			fileToClassification.put(linee[1], linee[0]);
		}
		sc.close();
		return fileToClassification;
	}

	/**
	 * Returns a map where
	 * 		key: directory path
	 * 		value: owner string, eg 'iw068'
	 * @since v1.0
	 * @param pathToUserDirList
	 * @return a mapping of directory to owner
	 * @throws FileNotFoundException
	 */
	public static HashMap<String, String> getUserDirectories(String pathToUserDirList) throws FileNotFoundException{
		HashMap<String, String> directoryToUser = new HashMap<>();
		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToUserDirList);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.charAt(0) == '#') {continue;}

			String[] linee = line.split("\\t");
			directoryToUser.put(linee[0], linee[1]);
		}
		sc.close();
		return directoryToUser;		
	}

	/**
	 * returns an map where 
	 * 		key: user name, eg 'iw068'
	 * 		value: 0L
	 * @since v1.0
	 * @param pathToNcduParsedFile
	 * @return a mapping of usernames to a size (initially set to 0L)
	 * @throws FileNotFoundException
	 */
	public static HashMap<String, Long> createUserToSize(String pathToNcduParsedFile) throws FileNotFoundException{
		HashMap<String, Long> userToSize = new HashMap<>();
		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToNcduParsedFile);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		line = sc.nextLine();
		String[] fields = line.split("\\t");
		int userFieldIndex = 3;
		for(int i = 0; i < fields.length; i++) {
			if(fields[i].equalsIgnoreCase(USER_FIELD_NAME)) {
				userFieldIndex = i;
			}
		}
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			String[] linee = line.split("\\t");
			userToSize.put(linee[userFieldIndex], 0L);
		}
		sc.close();
		return userToSize;		
	}

	/**
	 * looks for files in user directories and adds those file sizes to user's size score
	 * @since v1.0
	 * @param userToSize - HashMap<String, Long> key: user name 		value: bytes(long)
	 * @param directoryToUser - HashMap<String, String> key: user's personal directory		value: user name 
	 * @param blacklist - HashMap<String, String> key: file path		value: classification
	 * @param pathToNcduParseFile
	 * @throws FileNotFoundException
	 */
	public static void measureUserDirectories(HashMap<String, Long> userToSize, HashMap<String, String> directoryToUser, HashMap<String, String> blacklist, String pathToNcduParseFile) throws FileNotFoundException {

		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToNcduParseFile);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		line = sc.nextLine();
		String[] fields = line.split("\\t");
		//		int userFieldIndex = 3;
		int asizeFieldIndex = 1;
		for(int i = 0; i < fields.length; i++) {
			if(fields[i].equalsIgnoreCase(USER_FIELD_NAME)) {
				//				userFieldIndex = i;
			}
			if(fields[i].equalsIgnoreCase(ASIZE_FIELD_NAME)) {
				asizeFieldIndex = i;
			}
		}
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			String[] linee = line.split("\\t");
			Long asize = Long.parseLong(linee[asizeFieldIndex]);
			//			String user = linee[userFieldIndex];
			String filePath = linee[0];

			/*
			 * if the current line contains a file in a user directory
			 * 		if the current line contains a file not in the blacklist	
			 * 			add that file to user size
			 * 
			 */

			String currentUser = "";
			String directoryName = filePath.split("/")[3]; 
			if(directoryToUser.containsKey(directoryName)) {
				if(blacklist.containsKey(filePath) && (blacklist.get(filePath).equals(ARCHIVE) || blacklist.get(filePath).equals(IMPORTANT)) ) {
					continue;
				} 
				currentUser = directoryToUser.get(directoryName);
				if(userToSize.containsKey(currentUser)==false) userToSize.put(currentUser, 0L);
				userToSize.put(currentUser, userToSize.get(currentUser) + asize);
			} 

		}
		sc.close();
	}

	/**
	 * looks for files in shared project directories and adds those file sizes to user's size score
	 * @since v1.0
	 * @param userToSize - HashMap<String, Long> key: user name 		value: bytes(long)
	 * @param directoryToUser - HashMap<String, String> key: user's personal directory		value: user name 
	 * @param blacklist - HashMap<String, String> key: file path		value: classification
	 * @param pathToNcduParseFile
	 * @throws FileNotFoundException
	 */
	public static void measureProjectDirectories(HashMap<String, Long> userToSize, HashMap<String, String> directoryToUser, HashMap<String, String> blacklist, String pathToNcduParseFile) throws FileNotFoundException {

		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToNcduParseFile);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		line = sc.nextLine();
		String[] fields = line.split("\\t");
		int userFieldIndex = 3;
		int asizeFieldIndex = 1;
		for(int i = 0; i < fields.length; i++) {
			if(fields[i].equalsIgnoreCase(USER_FIELD_NAME)) {
				userFieldIndex = i;
			}
			if(fields[i].equalsIgnoreCase(ASIZE_FIELD_NAME)) {
				asizeFieldIndex = i;
			}
		}
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			String[] linee = line.split("\\t");
			Long asize = Long.parseLong(linee[asizeFieldIndex]);
			String user = linee[userFieldIndex];
			String filePath = linee[0];

			/*
			 * if the current line contains a file in a user directory
			 * 		if the current line contains a file not in the blacklist	
			 * 			add that file to user size
			 * 
			 */

			String currentUser = user;
			String directoryName = filePath.split("/")[3]; 
			if(directoryToUser.containsKey(directoryName) == false) { // only change from measureUserDirectories
				if(blacklist.containsKey(filePath) && (blacklist.get(filePath).equals(ARCHIVE) || blacklist.get(filePath).equals(IMPORTANT)) ) {
					continue;
				} 
				if(userToSize.containsKey(currentUser)==false) userToSize.put(currentUser, 0L);
				userToSize.put(currentUser, userToSize.get(currentUser) + asize);
			} 

		}
		sc.close();
	}



	/**
	 * Creates a file of users and their total file size 
	 * @param pathToBlackList
	 * @param pathToUserDirList
	 * @param pathToNcduParseFile
	 * @param outputPath
	 * @throws IOException
	 */
	public static void shame(String pathToBlackList, String pathToUserDirList, String pathToNcduParseFile, String outputPath) throws IOException {
		HashMap<String, String> blacklist = getBlacklist(pathToBlackList);
		HashMap<String, String> userDirectories = getUserDirectories(pathToUserDirList);
		HashMap<String, Long> userToSize = createUserToSize(pathToNcduParseFile);

		// i could do this in one pass instead of two for efficiency, however, I'd prefer this to be more modular, and I/O speed in reading files is slow but not that slow
		measureUserDirectories(userToSize, userDirectories, blacklist, pathToNcduParseFile);
		measureProjectDirectories(userToSize, userDirectories, blacklist, pathToNcduParseFile);

		// sort users based on size

		// print sorted users

		BufferedWriter output = null;
		File file = new File(outputPath);
		output = new BufferedWriter(new FileWriter(file));


		// https://www.geeksforgeeks.org/sorting-a-hashmap-according-to-values/
		Map<String, Long> hm1 = sortByValue(userToSize); 

		// print the sorted hashmap 
		for (Map.Entry<String, Long> en : hm1.entrySet()) { 
			//System.out.println("Key = " + en.getKey() +  ", Value = " + en.getValue()); 
			output.write(en.getKey() + "\t" + humanReadableByteCount(en.getValue()) + "\n");
		} 

		output.close();
	}




	/***************************************************************************************
	 ***************************************************************************************
	 ******************************** NCDU PARSER ******************************************
	 ***************************************************************************************
	 ***************************************************************************************
	 */


	/**
	 * Generates a map of user ID #'s to user name strings, eg 4672135 to iw068  
	 * @since v1.0
	 * @param pathToPasswd, typically /etc/passwd .older users who have been removed wil not show up in here
	 * @return a map of strings where key: UID value: user name
	 * @throws FileNotFoundException
	 */
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

	/**
	 * Parsers the JSON output file from NCDU
	 * @since v 1.0
	 * @param INPUT_PATH - The full file path to the ncdu JSON output file
	 * @param OUTPUT_BASE - The full output path to write the info file
	 * @throws IOException
	 */
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



	/***************************************************************************************
	 ***************************************************************************************
	 ******************************** UTILITIES ********************************************
	 ***************************************************************************************
	 ***************************************************************************************
	 */

	/**
	 * 
	 * @return returns the current JVM heap size in human readable format
	 */
	public static String getHeapSize() {
		return humanReadableByteCount(Runtime.getRuntime().totalMemory());
	}

	// function to sort hashmap by values 
	// https://www.geeksforgeeks.org/sorting-a-hashmap-according-to-values/
	public static HashMap<String, Long> sortByValue(HashMap<String, Long> hm) 
	{ 
		// Create a list from elements of HashMap 
		List<Map.Entry<String, Long> > list = 
				new LinkedList<Map.Entry<String, Long> >(hm.entrySet()); 

		// Sort the list 
		Collections.sort(list, new Comparator<Map.Entry<String, Long> >() { 
			public int compare(Map.Entry<String, Long> o1,  
					Map.Entry<String, Long> o2) 
			{ 
				return (o1.getValue()).compareTo(o2.getValue()); 
			} 
		}); 

		// put data from sorted list to hashmap  
		HashMap<String, Long> temp = new LinkedHashMap<String, Long>(); 
		for (Map.Entry<String, Long> aa : list) { 
			temp.put(aa.getKey(), aa.getValue()); 
		} 
		return temp; 
	} 

	/*
	 * https://programming.guide/java/formatting-byte-size-to-human-readable-format.html
	 * https://stackoverflow.com/questions/3758606/how-to-convert-byte-size-into-human-readable-format-in-java
	 */
	public static String humanReadableByteCount(String bytes_string) {
		long bytes = Long.parseLong(bytes_string);
		if (bytes < 1024) return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp-1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
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

	public static String getFileName(String str) {
		String file = str.split("\t")[0];
		String[] filee = file.split("/");
		int len = filee.length;
		return filee[len - 1];
	}



}
