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
		System.out.println(System.getProperty("java.version"));
		if(args[0].equals("-h") || args[0].equals("--help")) {
			System.out.println("java NCDU_Parser [ncdu_file] [output_file]");
		} else {
			ncduParser(args[0], args[1]);
		}
	}

	public static void ncduParser(String INPUT_PATH, String OUTPUT_BASE) throws IOException {

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
		//filePaths.push("");

		output.write("file\tasize\tdsize\n");

		HashMap<Integer, String> lineMap = new HashMap<Integer, String>();
		HashMap<Long, ArrayList<Integer>> sortMap = new HashMap<Long, ArrayList<Integer>>();
		int lineNumber = 1;

		Pattern ap = Pattern.compile("(?<=asize\":)\\d+");
		Pattern dp = Pattern.compile("(?<=dsize\":)\\d+");
		Pattern np = Pattern.compile("(?<=name\":\").+?(?=\")");
		int count = 0;

		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}

			Matcher am = ap.matcher(line); 
			Matcher dm = dp.matcher(line); 
			Matcher nm = np.matcher(line);  nm.find();

			if(line.contains("[")) {
				filePaths.push(nm.group());
				count = line.length() - line.replaceAll("]", "").length();
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
			sb.append("\n");
			String currPath = sb.toString();
			output.write(currPath);
			

			count = line.length() - line.replaceAll("]", "").length();
			for(int i = 0; i < count; i++) {
				if(filePaths.size() > 0)
					filePaths.pop();
			}

			lineMap.put(lineNumber, currPath);
			Long dSize = Long.parseLong(dsize);
			if(sortMap.containsKey(dSize)) {
				sortMap.get(dSize).add(lineNumber);
			} else {
				sortMap.put(dSize, new ArrayList<Integer>());
			}
			lineNumber++;
		}
		inputStream.close();
		sc.close();
		output.close();

		System.out.println("now sorting");

		output = null;
		file = new File(OUTPUT_BASE + "_sorted.tsv");
		output = new BufferedWriter(new FileWriter(file));

		output.write("file\tasize\tdsize\n");
		ArrayList<Long> sortIndex = new ArrayList<>(sortMap.keySet());
		Collections.sort(sortIndex, Collections.reverseOrder());

		System.out.println("done sorting, now writing");

		for(Long si : sortIndex) {
			for(Integer ln : sortMap.get(si)) {
				output.write(lineMap.get(ln));
			}
		}
		output.close();

	}

}

