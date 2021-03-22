package gCNV;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class gCNVHelperTool {

	
	public String[] args;
	
	public gCNVHelperTool(String[] args) {
		this.args = args;
	}
	
	
	
	public void run() throws IOException, InterruptedException {
		
	}
	
	/**
	 * returns the int median of a double array, for reasons
	 * @param i
	 * @return
	 */
	public String median(double[] i) {
		Arrays.sort(i);
		//		return Integer.toString((int)(i[i.length/2]));
		return Integer.toString((int)(
				i.length%2!=0 ? 
						i[i.length/2] : 
							(i[(i.length-1)/2] + i[i.length/2])/2.0 + 0.5
				));
	}
	
	/**
	 * reads in a GTF format file,
	 * 
	 * @param INPUT_GTF
	 * @throws FileNotFoundException
	 */
	public static ArrayList<Gene> parseGTFFile(String INPUT_GTF) throws FileNotFoundException {
		print("reading annotation file");

		HashMap<String, ArrayList<Integer>> geneToStart = new HashMap<>();
		HashMap<String, ArrayList<Integer>> geneToEnd = new HashMap<>();
		HashMap<String, String> geneToChr = new HashMap<>();

		int gtfChr = 0;
		int gtfFeatureType = 2;
		int gtfStart = 3;
		int gtfEnd = 4;

		Pattern geneNamePattern = Pattern.compile("(?<=gene_name \").+?(?=\")");
		Pattern geneTypePattern = Pattern.compile("(?<=gene_type \").+?(?=\")");

		FileInputStream inputStream = new FileInputStream(INPUT_GTF);
		Scanner gtf = new Scanner(inputStream, "UTF-8");
		while (gtf.hasNext()) {
			String line = gtf.nextLine();
			if (line.startsWith("#"))
				continue;
			String[] linee = line.split("\\t");

			Matcher geneTypeMatcher = geneTypePattern.matcher(line);
			String geneType = geneTypeMatcher.find() ? geneTypeMatcher.group() : "-1";
			String featureType = linee[gtfFeatureType];

			if (geneType.equals("protein_coding") && featureType.equals("exon")) {

				Matcher geneNameMatcher = geneNamePattern.matcher(line);
				String geneName = geneNameMatcher.find() ? geneNameMatcher.group() : "-1";

				if (geneToStart.containsKey(geneName) == false) {
					geneToStart.put(geneName, new ArrayList<>());
					geneToEnd.put(geneName, new ArrayList<>());
					geneToChr.put(geneName, linee[gtfChr]);
				}
				geneToStart.get(geneName).add(Integer.parseInt(linee[gtfStart]));
				geneToEnd.get(geneName).add(Integer.parseInt(linee[gtfEnd]));
			}
		}
		gtf.close();
		gtf = null;

		print("analyzing annotations");

		ArrayList<Gene> genes = new ArrayList<>();
		for (String gene : geneToChr.keySet()) {

			HashMap<Integer, Integer> startIdToIndex = new HashMap<>();
			var currentStarts = geneToStart.get(gene);
			for (int i = 0; i < currentStarts.size(); i++) {
				if (!startIdToIndex.containsKey(currentStarts.get(i))) {
					startIdToIndex.put(currentStarts.get(i), i);
				} else {
					int oldEnd = geneToEnd.get(gene).get(startIdToIndex.get(currentStarts.get(i)));
					int newEnd = geneToEnd.get(gene).get(i);
					if (newEnd > oldEnd) {
						startIdToIndex.put(currentStarts.get(i), i);
					} // else do nothing
				}
			}

			Collections.sort(currentStarts);
			ArrayList<Integer> sortedEnds = new ArrayList<>();
			var currentEnds = geneToEnd.get(gene);
			for (int i = 0; i < currentStarts.size(); i++) {
				sortedEnds.add(currentEnds.get(startIdToIndex.get(currentStarts.get(i))));

			}

			genes.add(new Gene(gene, geneToChr.get(gene), currentStarts, sortedEnds));
		}

		Collections.sort(genes);
		print("read " + genes.size() + " genes");
		return genes;
	}

	
	public static int BinarySearchLowerBound(ArrayList<Integer> list, int value, int start, int end) {
		int low = Math.max(0, start);
		int high = end==-1 ? list.size()-1 : Math.min(list.size()-1, end);

		int ol = low;
		int oh = high;

		while (low <= high) {
			int mid = (low + high) >>> 1;
		int midVal = list.get(mid);
		int cmp = Integer.compare(midVal, value);
		if (cmp < 0)
			low = mid + 1;
		else if (cmp > 0)
			high = mid - 1;
		else {
			while(mid>0 && list.get(mid) >= value) {
				mid--;
			}
			return mid;
		}
		}
		int mid = Math.min(Math.max(ol, low), oh);
		if(value >= list.get(oh)) {
			value = list.get(oh);
		}
		while(mid>0 && list.get(mid) >= value) {
			mid--;
		}
		return mid;
	}
	
	public static int BinarySearchUpperBound(ArrayList<Integer> list, int value, int start, int end) {
		int low = Math.max(0, start);
		int high = end==-1 ? list.size()-1 : Math.min(list.size()-1, end);

		int ol = low;
		int oh = high;

		while (low <= high) {
			int mid = (low + high) >>> 1;
		int midVal = list.get(mid);
		int cmp = Integer.compare(midVal, value);
		if (cmp < 0)
			low = mid + 1;
		else if (cmp > 0)
			high = mid - 1;
		else {
			while(mid<oh && list.get(mid) <= value) {
				mid++;
			}
			return mid;
		}
		}
		int mid = Math.min(Math.max(ol, low), oh);
		while(mid<oh && list.get(mid) <= value) {
			mid++;
		}
		return mid;
	}
	

	/**
	 * @return returns the current JVM heap size in human readable format
	 */
	public static String getHeapSize() {
		return humanReadableByteCount(Runtime.getRuntime().totalMemory());
	}

	/**
	 * @param bytes
	 * @return human readable version of bytes
	 */
	public static String humanReadableByteCount(long bytes) {
		if (bytes < 1024)
			return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp - 1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}

	// returns the mean of an arraylist, using integers, so will trucnate
	public static int mean(ArrayList<Integer> al) {
		int sum = 0;
		for (Integer i : al) {
			sum += i;
		}
		return (sum / al.size());
	}

	/**
	 * returns the median of an integer arraylist
	 * 
	 * @param al
	 * @return
	 */
	public static int median(ArrayList<Integer> al) {
		Collections.sort(al);
		return (al.get(al.size() / 2));
	}

	/**
	 * returns to string arraylist with values separated by sep
	 * 
	 * @param arrayList
	 * @param sep
	 * @return
	 */
	public String arrayToStringWithSep(ArrayList<String> arrayList, String sep) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < arrayList.size() - 1; i++) {
			sb.append(arrayList.get(i));
			sb.append(sep);
		}
		sb.append(arrayList.get(arrayList.size() - 1));
		return sb.toString();
	}

	public static void print(String str) {
		System.out.println(str);
	}

	public static void print(String[] str) {
		for (String s : str) {
			System.out.print(s + "\t");
		}
		System.out.println();
	}

	public static String getTimeHumanReadable(long time) {
		long seconds = time % 60;
		long minutes = (time / 60) % 60;
		long hours = (time / 3600);

		return "Hours: " + hours + "  Min: " + minutes + "  Sec: " + seconds;
	}
	

	/**
	 * takes in a file containing a list of files to download. downloads those
	 * files. this has only been tested on windows, it is VERY DELICATE
	 * 
	 * @param filesFile
	 * @throws InterruptedException
	 */
	public static void downloadFiles(String filesFile, String pathToDownloadToUnix) throws InterruptedException {
		String[] command;
		if (System.getProperty("os.name").contains("indows")) {
			String[] dosCommand = { "bash", "-c", "'cat", filesFile, "|", "gsutil", "-m", "cp", "-I",
					pathToDownloadToUnix + "/", "'" };
			command = dosCommand;

			System.out.println(String.join(" ", command));

			try {
				new ProcessBuilder(command).inheritIO().start().waitFor();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
			System.out.println(
					"WARNING: this download functionality has not been tested on macOS/linux, please report any issues");
			String[] unixCommand = {"/bin/sh", "-c", "cat " + filesFile + " | gsutil -m cp -I " +  pathToDownloadToUnix + "/" };
			command = unixCommand;

			System.out.println(String.join(" ", command));

			try {
				new ProcessBuilder(command).inheritIO().start().waitFor();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

	}
	
	/**
	 * determine if a string is an integer using Java's Integer parse method while
	 * this is not the most efficient way to check for integer-ness, I like that it
	 * handles some exceptions for me, as I will be calling parseInt() eventually
	 * 
	 * @param input
	 * @return
	 */
	public static boolean isInteger(String input) {
		try {
			Integer.parseInt(input);
			return true;
		} catch (NumberFormatException nfe) {
			return false;
		}
	}


	/**
	 * https://stackoverflow.com/questions/852665/command-line-progress-bar-in-java
	 * 
	 * @param remain
	 * @param total
	 */
	public static void progressPercentage(int remain, int total, String message) {
		if (remain > total) {
			throw new IllegalArgumentException();
		}
		int maxBareSize = 10; // 10unit for 100%
		int remainProcent = ((100 * remain) / total) / maxBareSize;
		char defaultChar = '-';
		String icon = "*";
		String bare = new String(new char[maxBareSize]).replace('\0', defaultChar) + "]";
		StringBuilder bareDone = new StringBuilder();
		bareDone.append("[");
		for (int i = 0; i < remainProcent; i++) {
			bareDone.append(icon);
		}
		String bareRemain = bare.substring(remainProcent, bare.length());
		System.out.print("\r" + bareDone + bareRemain + " " + remainProcent * 10 + "%  " + remain + " / " + total + "  "
				+ message);
		if (remain == total) {
			System.out.print("\n");
		}
	}
	

	/**
	 * returns the mean and median of the integer columns for each unique ID of
	 * given column
	 * 
	 * @param input
	 * @param sampleColumn
	 * @throws IOException
	 */
	public ArrayList<DataFrame> getPerSampleMetrics(String INPUT_PATH, String fieldColumn, String OUTPUT_PATH,
			boolean writeFile) throws IOException {
		DataFrame gcnv = new DataFrame(INPUT_PATH, true, "\\t", "#");
		System.out.println(gcnv.nrow() + " rows");
		System.out.println(gcnv.ncol() + " columns");
		System.out.println(gcnv.columnMapping.toString());

		// key: sample name
		// value: hashmap of column names and values
		HashMap<String, HashMap<String, ArrayList<Integer>>> perFieldMap = new HashMap<>();
		HashMap<String, HashMap<String, Integer>> perFieldCalculations = new HashMap<>();

		// data needs at least one entry or will die
		String EXAMPLE_VALUE = gcnv.get(fieldColumn, 0);

		// get a list of all the integer columns
		// note: this assumes a populated field
		ArrayList<String> columnsToMeasure = new ArrayList<>();
		for (String field : gcnv.fieldNames) {
			if (field.equalsIgnoreCase("chr") || field.equalsIgnoreCase("start") || field.equalsIgnoreCase("end")) {
				continue;
			}
			if (isInteger(gcnv.get(field, 0))) {
				columnsToMeasure.add(field);
			}
		}

		// populate hashmap
		for (int i = 0; i < gcnv.nrow(); i++) {
			String currField = gcnv.get(fieldColumn, i); // current field looked at
			if (!perFieldMap.containsKey(currField)) { // if the field value has not been observed yet
				perFieldMap.put(currField, new HashMap<>()); // add the field observation and a new hashmap to the main
				// hashmap
				for (String field : columnsToMeasure) { // populate that value's hashmap with empty arraylists
					perFieldMap.get(currField).put(field, new ArrayList<>());
				}
			}

			// add the value to the appropriate hashmap value
			for (String field : columnsToMeasure) {
				try {
					perFieldMap.get(currField).get(field).add(Integer.parseInt(gcnv.get(field, i)));
				} catch (NumberFormatException nfe) {
					try {
						perFieldMap.get(currField).get(field).add((int) Double.parseDouble(gcnv.get(field, i)));
					} catch (NumberFormatException nfe2) {
						perFieldMap.get(currField).get(field).add(-1);
					}
				}
			}
		}

		// do the calculations
		for (String uniqueValue : perFieldMap.keySet()) {
			perFieldCalculations.put(uniqueValue, new HashMap<>());
			for (String field : columnsToMeasure) {
				Integer mean = mean(perFieldMap.get(uniqueValue).get(field));
				Integer median = median(perFieldMap.get(uniqueValue).get(field));

				perFieldCalculations.get(uniqueValue).put(field + "_MEAN", mean);
				perFieldCalculations.get(uniqueValue).put(field + "_MEDIAN", median);
			}
		}

		// populate data frame values
		ArrayList<ArrayList<String>> values = new ArrayList<>();
		ArrayList<String> columnNames = new ArrayList<>();
		columnNames.addAll(perFieldCalculations.get(EXAMPLE_VALUE).keySet());
		for (int i = 0; i < columnNames.size(); i++) {
			values.add(new ArrayList<>());
		}
		ArrayList<String> nameColumn = new ArrayList<>();
		for (String uniqueValue : perFieldMap.keySet()) {
			nameColumn.add(uniqueValue);
			for (int i = 0; i < columnNames.size(); i++) {
				values.get(i).add(String.valueOf(perFieldCalculations.get(uniqueValue).get(columnNames.get(i))));
			}
		}

		values.add(nameColumn);
		columnNames.add(fieldColumn);
		DataFrame toWrite = new DataFrame();
		toWrite.addColumns(columnNames, values);

		if (writeFile == true) {
			toWrite.writeFile(OUTPUT_PATH, true);
		}

		ArrayList<DataFrame> rList = new ArrayList<>();
		rList.add(gcnv);
		rList.add(toWrite);
		return rList;
	}
	

}
