package data_management;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
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

	public static HashMap<String, String> getBlacklist(String pathToBlacklist) throws FileNotFoundException {
		HashMap<String, String> fileToClassification = new HashMap<>();
		FileInputStream inputStream = null;
		Scanner sc = null;

		inputStream = new FileInputStream(pathToBlacklist);
		sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.charAt(0) == '#') {continue;}

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
			if(line.charAt(0) == '#') {continue;}

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
			if(line.charAt(0) == '#') {continue;}

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

			String currentUser = "";
			String directoryName = filePath.split("/")[2]; // /data/talkowski/iwong/some_file.txt TODO: test
			if(directoryToUser.containsKey(directoryName)) {
				if(blacklist.containsKey(filePath) && (blacklist.get(filePath).equals(ARCHIVE) || blacklist.get(filePath).equals(IMPORTANT)) ) {
					continue;
				} 
				currentUser = directoryToUser.get(directoryName);
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
			if(line.charAt(0) == '#') {continue;}

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

			String currentUser = "";
			String directoryName = filePath.split("/")[2]; // /data/talkowski/iwong/some_file.txt TODO: test
			if(directoryToUser.containsKey(directoryName) == false) { // only change from measureUserDirectories
				if(blacklist.containsKey(filePath) && (blacklist.get(filePath).equals(ARCHIVE) || blacklist.get(filePath).equals(IMPORTANT)) ) {
					continue;
				} 
				currentUser = directoryToUser.get(directoryName);
				userToSize.put(currentUser, userToSize.get(currentUser) + asize);
			} 

		}
		sc.close();
	}



	public static void shame(String pathToBlackList, String pathToUserDirList, String pathToNcduParseFile) throws FileNotFoundException {
		HashMap<String, String> blacklist = getBlacklist(pathToBlackList);
		HashMap<String, String> userDirectories = getUserDirectories(pathToUserDirList);
		HashMap<String, Long> userToSize = createUserToSize(pathToNcduParseFile);

		// i could do this in one pass instead of two for efficiency, however, I'd prefer this to be more modular, and I/O speed in reading files is slow but not that slow
		measureUserDirectories(userToSize, userDirectories, blacklist, pathToNcduParseFile);
		measureProjectDirectories(userToSize, userDirectories, blacklist, pathToNcduParseFile);

		// sort users based on size

		// print sorted users


		// https://www.geeksforgeeks.org/sorting-a-hashmap-according-to-values/
		Map<String, Long> hm1 = sortByValue(userToSize); 

		// print the sorted hashmap 
		for (Map.Entry<String, Long> en : hm1.entrySet()) { 
			System.out.println("Key = " + en.getKey() +  
					", Value = " + en.getValue()); 
		} 



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





}
