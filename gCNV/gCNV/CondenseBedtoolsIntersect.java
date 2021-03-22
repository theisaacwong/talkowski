package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

public class CondenseBedtoolsIntersect extends gCNVHelperTool {

	public CondenseBedtoolsIntersect(String[] args) {
		super(args);
	}

	@Override
	public void run() throws IOException {
		String INPUT_PATH = args[1];
		String OUTPUT_PATH = args[2];
		String columnsToHashOnString = args[3];
		String columnsToMergeString = args[4];
		String columnsToKeepString = args[5];
		this.condenseBedtoolsIntersect(INPUT_PATH, OUTPUT_PATH, columnsToHashOnString, columnsToMergeString, columnsToKeepString);
	}
	

	/**
	 * is string builder creation slower or faster than multiple string.equals()
	 * calls? I'm currently thinking that string.equals() calls might be faster?
	 * I'll probably do both just to see
	 * 
	 * @param INPUT_PATH
	 * @param OUTPUT_PATH
	 * @param columnsToHashOnString
	 * @param columnsToMergeString
	 * @throws IOException
	 */
	public void condenseBedtoolsIntersect(String INPUT_PATH, String OUTPUT_PATH, String columnsToHashOnString,
			String columnsToMergeString, String columnsToKeepString) throws IOException {
		int[] columnsToHash = Arrays.asList(columnsToHashOnString.split(",")).stream().mapToInt(Integer::parseInt)
				.toArray();
		int[] columnsToMerge = Arrays.asList(columnsToMergeString.split(",")).stream().mapToInt(Integer::parseInt)
				.toArray();
		int[] columnsToKeep = Arrays.asList(columnsToKeepString.split(",")).stream().mapToInt(Integer::parseInt)
				.toArray();
		String tab = "\t";
		String newline = "\n";
		int lastColumnToMerge = columnsToMerge[columnsToMerge.length - 1];

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = new FileInputStream(INPUT_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		Map<Integer, Set<String>> columnToMergedStrings = new HashMap<>();
		for (int i : columnsToMerge) {
			columnToMergedStrings.put(i, new HashSet<>());
		}

		String[] previousLinee = sc.nextLine().split("\t"); // previous hash?
		int n_merged = 1;

		while (sc.hasNextLine()) {
			String currentLine = "";
			try {
				currentLine = sc.nextLine();
			} catch (Exception e) {
				System.out.println(e);
			}

			boolean isSame = true;
			String[] currentLinee = currentLine.split("\t");

			// not sure if this is better or worse than comparing concatenated strings to
			// each other
			for (int i : columnsToHash) {
				if (!previousLinee[i].equals(currentLinee[i])) {
					isSame = false;
					break;
				}
			}

			if (isSame) {
				for (int i : columnsToMerge) {
					columnToMergedStrings.get(i).add(currentLinee[i]);
				}
				n_merged++;
			} else {
				StringBuilder lineToWrite = new StringBuilder();
				for (int i : columnsToKeep) {
					lineToWrite.append(previousLinee[i]);
					lineToWrite.append(tab);
				}
				for (int i : columnsToMerge) {
					List<String> values = new ArrayList<>();
					values.addAll(columnToMergedStrings.get(i));
					Collections.sort(values);
					if (values.size() == 0) {
						values.add("None");
					}
					lineToWrite.append(String.join(",", values));
					if (i != lastColumnToMerge) {
						lineToWrite.append(tab);
					}
				}

				lineToWrite.append(tab);
				lineToWrite.append(n_merged);
				lineToWrite.append(newline);

				n_merged = 1;

				output.write(lineToWrite.toString());

				previousLinee = currentLinee;
				for (int i : columnsToMerge) {
					columnToMergedStrings.put(i, new HashSet<>());
				}
			}
		}

		output.close();
		sc.close();
	}

	
}
