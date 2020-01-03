package data_management;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

public class SHAME {

	public static final String ARCHIVE = "archive";
	public static final String IMPORTANT = "important";
	public static final String USER_FIELD_NAME = "user";
	public static final String ASIZE_FIELD_NAME = "asize";

	
	public static void main(String[] args) throws IOException {

//		/String wd = "C:/Users/iwong/Documents/MGH/scripts/";
		//args = new String[] {wd+"blacklist.txt", wd+"personalDirectories.txt", wd+"ncdutmpout.txt.tsv", wd+"shameout.txt"};
		System.out.println("Java version " + System.getProperty("java.version"));
		if(args.length == 4) {
			// shame(String pathToBlackList, String pathToUserDirList, String pathToNcduParseFile, String outputPath) 
			shame(args[0], args[1], args[2], args[3]);
		}
		else {
			System.out.println("args\n1: pathToBlackList\n2: pathToUserDirLst\n3: pathToNcduParseFile\n4: outputPath");
		}
	}
	
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
			String directoryName = filePath.split("/")[3]; // /data/talkowski/iwong/some_file.txt TODO: test
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
			String directoryName = filePath.split("/")[3]; // /data/talkowski/iwong/some_file.txt TODO: test
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





}
