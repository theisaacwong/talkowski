package gCNV;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class Defragment extends gCNVHelperTool {

	public Defragment(String[] args) {
		super(args);
	}

	@Override
	public void run() throws IOException {
		String match_output = args[1];
		String defragmentation_output = args[2];
		if(args.length==3) {
			this.defragment(match_output, defragmentation_output);	
		} else {
			String filteredIntervals = args[3];
			this.defragment3(match_output, defragmentation_output, filteredIntervals);
		}
	}
	

	public void defragment3(String input, String output, String filtered_intervals) throws IOException {
		double EXTENSION_VALUE = 0.4; // by what percent of cnv width to extend look behind/ahead distance for overlapping cnvs
		int MINIMUM_QS = 20;
		int CHR = 0;
		int START = 1;
		int END = 2;
		DataFrame gcnv = new DataFrame(input, true, "\\t", "@");
		gcnv.fieldNames[0] = "chr";
		gcnv.columnMapping.put("chr", 0);
		gcnv.columnMapping.remove("#chrom");
		gcnv.sort();

		print("reading filtered intervals");
		//Read filtered intervals, column 1 is name of cluster, column 2 is file
		DataFrame filteredIntervalsManifest = new DataFrame(filtered_intervals, true, "\\t", "#");
		HashMap<String, DataFrame> batchToFilteredIntervalDF = new HashMap<>();
		for(int i = 0; i < filteredIntervalsManifest.nrow(); i++) {
			String batchName = filteredIntervalsManifest.get(i, 0);
			String fiPath = filteredIntervalsManifest.get(i, 1);
			//check that file exists
			if(new File(fiPath).isFile()) {
				batchToFilteredIntervalDF.put(batchName, new DataFrame(fiPath, false, "\\t", "@"));
			}
		}
		print("done");

		// get batch indexes
		HashMap<String, ArrayList<Integer>> batchToIndexes = new HashMap<>(); 
		for(int i = 0; i < gcnv.nrow(); i++) {
			boolean passQS = (Integer.parseInt(gcnv.get("QS", i)) >= MINIMUM_QS);
			if(passQS == false) {
				continue;
			}

			if(!batchToIndexes.containsKey(gcnv.get("batch", i))) {
				batchToIndexes.put(gcnv.get("batch", i), new ArrayList<>());
			}
			batchToIndexes.get(gcnv.get("batch", i)).add(i);
		}

		// check to make sure each 
		for(String batch : batchToIndexes.keySet()) {
			if(!batchToFilteredIntervalDF.containsKey(batch)) {
				print("Warning! " + batch + " does not have a filtered interval");
			}
		}

		ArrayList<Integer> allIndexesDefragged = new ArrayList<>();
		DataFrame allDefragmentedCalls = new DataFrame(gcnv.fieldNames);

		//bin == intervals, 
		for(String batch : batchToIndexes.keySet()) {
			print(batch);
			DataFrame currBin = batchToFilteredIntervalDF.get(batch);

			//			ArrayList<Integer> extendedBinStarts = new ArrayList<>();
			//			ArrayList<Integer> extendedBinEnds = new ArrayList<>();

			HashMap<Integer, Integer> gcnvIndexToExtendedBinStartCoords = new HashMap<>();
			HashMap<Integer, Integer> gcnvIndexToExtendedBinEndCoords = new HashMap<>();

			//convert to intervals coordinates
			HashMap<String, Integer> gCNVtoStartBin = new HashMap<>();
			HashMap<String, Integer> gCNVtoEndBin = new HashMap<>();
			HashMap<String, ArrayList<Integer>> binChrToStarts = new HashMap<>();
			HashMap<String, ArrayList<Integer>> binChrToEnds = new HashMap<>();

			// read the interval/bin file, storing chr-start and chr-end pairings, sort ascending
			for(int i = 0; i < currBin.nrow(); i++) {
				String currChr = currBin.get(i, CHR);
				String currStart = currBin.get(i, START);
				String currEnd = currBin.get(i, END);
				if(!binChrToStarts.containsKey(currChr)) {
					binChrToStarts.put(currChr, new ArrayList<>());
					binChrToEnds.put(currChr, new ArrayList<>());
				}
				binChrToStarts.get(currChr).add(Integer.parseInt(currStart));
				binChrToEnds.get(currChr).add(Integer.parseInt(currEnd));
			}

			for(String chr : binChrToStarts.keySet()) {
				Collections.sort(binChrToStarts.get(chr));
				Collections.sort(binChrToEnds.get(chr));
			}

			// convert each gcnv start/end to bin start/end by finding closest value less than/greater than, then extend
			for(Integer i : batchToIndexes.get(batch)) {
				String currChr = gcnv.get("chr", i);
				Integer currEnd = Integer.parseInt(gcnv.get("end", i));
				Integer currStart = Integer.parseInt(gcnv.get("start", i));

				String currKeyEnd = currChr + "_" + currEnd;
				if(!gCNVtoEndBin.containsKey(currKeyEnd)) {
					int rawIndex = Collections.binarySearch(binChrToEnds.get(currChr), currEnd);
					int index = rawIndex >= 0 ? rawIndex : Math.min(binChrToEnds.get(currChr).size()-1, -1 - rawIndex );
					Integer closestBinEnd = binChrToEnds.get(currChr).get(index);
					gCNVtoEndBin.put(currKeyEnd, closestBinEnd);
				} 

				String currKeyStart = currChr + "_" + currStart;
				if(!gCNVtoStartBin.containsKey(currKeyStart)) {
					int rawIndex = Collections.binarySearch(binChrToStarts.get(currChr), currStart);
					int index = rawIndex >= 0 ? rawIndex : Math.max(0, -2 - rawIndex);
					Integer closestBinStart = binChrToStarts.get(currChr).get(index);
					gCNVtoStartBin.put(currKeyStart, closestBinStart);
				} 

				//width and extension calculation
				double width = gCNVtoEndBin.get(currKeyEnd) - gCNVtoStartBin.get(currKeyStart) + 1;
				//				extendedBinStarts.add((int)(gCNVtoStartBin.get(currKeyStart) - (width * EXTENSION_VALUE)));
				//				extendedBinEnds.add((int) (gCNVtoEndBin.get(currKeyEnd) + (width * EXTENSION_VALUE)));

				gcnvIndexToExtendedBinStartCoords.put(i, (int)(gCNVtoStartBin.get(currKeyStart) - (width * EXTENSION_VALUE)));
				gcnvIndexToExtendedBinEndCoords.put(i, (int) (gCNVtoEndBin.get(currKeyEnd) + (width * EXTENSION_VALUE)));

			}

			// for each unique sample, get map of sample to indexes
			HashMap<String, ArrayList<Integer>> sampleToIndexes = new HashMap<>();
			for(Integer i : batchToIndexes.get(batch)) {
				String currSample = gcnv.get("sample", i);
				if (!sampleToIndexes.containsKey(currSample)) {
					sampleToIndexes.put(currSample, new ArrayList<>());
				}
				sampleToIndexes.get(currSample).add(i);
			}

			DataFrame batchDefragmentedCalls = new DataFrame(gcnv.fieldNames);
			HashSet<Integer> batchIndexesDefragmented = new HashSet<>();


			// Light help you, the below code of Shai'tan's own
			// main defragmentation code - currently not super worried about runtime
			// sampleLoop:
			for (String currentSample : sampleToIndexes.keySet()) {

				// for each CNV
				// cnvLoop:
				for (int i = 0; i < sampleToIndexes.get(currentSample).size();) {
					int x = sampleToIndexes.get(currentSample).get(i); // convert to gcnv row index

					ArrayList<Integer> currentFragment = new ArrayList<>();
					currentFragment.add(x);

					// keep track of largest/current fragment end
					//					int currentFragmentEnd = extendedBinEnds.get(x);
					int currentFragmentEnd = gcnvIndexToExtendedBinEndCoords.get(x);
					int k = ++i; // lookahead amount for current cnv fragment
					while (k < sampleToIndexes.get(currentSample).size()) {
						int y = sampleToIndexes.get(currentSample).get(k);
						// check for overlap
						if (gcnv.get("chr", x).equals(gcnv.get("chr", y)) && gcnv.get("CN", x).equals(gcnv.get("CN", y))
								//								&& currentFragmentEnd >= extendedBinStarts.get(y)
								&& currentFragmentEnd >= gcnvIndexToExtendedBinStartCoords.get(y)
								&& Integer.parseInt(gcnv.get("start", x)) <= Integer.parseInt(gcnv.get("end", y))) {

							currentFragment.add(y);
							//							currentFragmentEnd = Math.max(currentFragmentEnd, extendedBinEnds.get(y));
							currentFragmentEnd = Math.max(currentFragmentEnd, gcnvIndexToExtendedBinEndCoords.get(y));
							batchIndexesDefragmented.add(y);
							k++;
						} else {
							//							System.out.println("nope!");
							break;
						}
					}
					i = k;

					if (currentFragment.size() != 1) {
						batchIndexesDefragmented.add(x);
						int newStart = Integer.parseInt(gcnv.get("start", currentFragment.get(0)));
						int newEnd = 0;
						int newQA = 0;
						int newQS = 0;
						String newQSS = gcnv.get("QSS", currentFragment.get(0));
						String newQSE = gcnv.get("QSE", currentFragment.get(currentFragment.size() - 1));
						for (int n = 0; n < currentFragment.size(); n++) {
							newStart = Math.min(newStart, Integer.parseInt(gcnv.get("start", currentFragment.get(n))));
							newEnd = Math.max(newEnd, Integer.parseInt(gcnv.get("end", currentFragment.get(n))));
							newQA += Integer.parseInt(gcnv.get("QA", currentFragment.get(n)));
							newQS = Math.max(newQS, Integer.parseInt(gcnv.get("QS", currentFragment.get(n))));
						}
						newQA /= currentFragment.size();
						String[] newRow = new String[gcnv.fieldNames.length];
						for (int s = 0; s < newRow.length; s++) {
							newRow[s] = gcnv.get(x)[s];
						}
						newRow[gcnv.columnMapping.get("start")] = Integer.toString(newStart);
						newRow[gcnv.columnMapping.get("end")] = Integer.toString(newEnd);
						newRow[gcnv.columnMapping.get("QA")] = Integer.toString(newQA);
						newRow[gcnv.columnMapping.get("QS")] = Integer.toString(newQS);
						newRow[gcnv.columnMapping.get("QSS")] = newQSS;
						newRow[gcnv.columnMapping.get("QSE")] = newQSE;
						batchDefragmentedCalls.add(newRow);
					} 

				}
			}

			allIndexesDefragged.addAll(batchIndexesDefragmented);
			if(allDefragmentedCalls.df.size() == 0) {
				allDefragmentedCalls = batchDefragmentedCalls;
			} else {
				allDefragmentedCalls.rbind(batchDefragmentedCalls);
			}

		}


		Collections.sort(allIndexesDefragged, Collections.reverseOrder());
		for (int i = 0; i < allIndexesDefragged.size(); i++) {
			gcnv.df.remove((int) (allIndexesDefragged.get(i)));
		}

		ArrayList<String> boolDefragmented_1 = new ArrayList<>();
		for (int i = 0; i < gcnv.nrow(); i++) {
			boolDefragmented_1.add("FALSE");
		}
		gcnv.addColumn("defragmented", boolDefragmented_1);

		ArrayList<String> boolDefragmented_2 = new ArrayList<>();
		for (int i = 0; i < allDefragmentedCalls.nrow(); i++) {
			boolDefragmented_2.add("TRUE");
		}
		allDefragmentedCalls.addColumn("defragmented", boolDefragmented_2);

		gcnv.rbind(allDefragmentedCalls);
		gcnv.sort();
		gcnv.writeFile(output, true);

	}

	/**
	 * 
	 * @param match_output
	 * @param defragmentation_output
	 * @throws IOException
	 */
	public void defragment(String match_output, String defragmentation_output) throws IOException {
		double EXTENSION_VALUE = 0.4;
		int MINIMUM_QS = 20;
		DataFrame gcnv = new DataFrame(match_output, true, "\\t", "@");
		gcnv.fieldNames[0] = "chr";
		gcnv.columnMapping.put("chr", 0);
		gcnv.columnMapping.remove("#chrom");
		gcnv.sort();



		// extend coordinates
		ArrayList<Integer> extendedStarts = new ArrayList<>();
		ArrayList<Integer> extendedEnds = new ArrayList<>();
		// qs pass threshold?
		ArrayList<Boolean> passQS = new ArrayList<>();
		// other metrics
		HashMap<String, Integer> chrToFirstIndex = new HashMap<>();
		HashMap<String, Integer> chrToLastIndex = new HashMap<>();
		for (int i = 0; i < gcnv.nrow(); i++) {
			// extend coordinates
			int currentStart = Integer.parseInt(gcnv.get("start", i));
			int currentEnd = Integer.parseInt(gcnv.get("end", i));
			double width = currentEnd - currentStart + 1;
			int extensionStart = (int) (currentStart - (width * EXTENSION_VALUE));
			int extensionEnd = (int) (currentEnd + (width * EXTENSION_VALUE));
			extendedStarts.add(extensionStart);
			extendedEnds.add(extensionEnd);

			// qs pass threshold
			passQS.add(Integer.parseInt(gcnv.get("QS", i)) >= MINIMUM_QS);
			if (!chrToFirstIndex.containsKey(gcnv.get("chr", i))) {
				chrToFirstIndex.put(gcnv.get("chr", i), i);
			}
			chrToLastIndex.put(gcnv.get("chr", i), i);
		}

		// for each unique sample, get map of sample to indexes
		HashMap<String, ArrayList<Integer>> sampleToIndexes = new HashMap<>();
		for (int i = 0; i < gcnv.nrow(); i++) {

			if (!passQS.get(i)) {
				continue;
			}

			String currSample = gcnv.get("sample", i);
			if (!sampleToIndexes.containsKey(currSample)) {
				sampleToIndexes.put(currSample, new ArrayList<>());
			}
			sampleToIndexes.get(currSample).add(i);
		}

		DataFrame newDefragmentedCalls = new DataFrame(gcnv.fieldNames);
		HashSet<Integer> indexesDefragmented = new HashSet<>();

		// Light help you, the below code of Shai'tan's own
		// main defragmentation code - currently not super worried about runtime
		// sampleLoop:
		for (String currentSample : sampleToIndexes.keySet()) {

			// for each CNV
			// cnvLoop:
			for (int i = 0; i < sampleToIndexes.get(currentSample).size();) {
				int x = sampleToIndexes.get(currentSample).get(i); // convert to gcnv row index

				ArrayList<Integer> currentFragment = new ArrayList<>();
				currentFragment.add(x);

				int currentFragmentEnd = extendedEnds.get(x);
				int k = ++i;
				while (k < sampleToIndexes.get(currentSample).size()) {
					int y = sampleToIndexes.get(currentSample).get(k);
					//					System.out.println("comparing " + x + " vs " + y + " " + currentFragmentEnd + " > " + extendedStarts.get(y) + " = " + (currentFragmentEnd>extendedStarts.get(y)));
					if (gcnv.get("chr", x).equals(gcnv.get("chr", y)) && gcnv.get("CN", x).equals(gcnv.get("CN", y))
							&& currentFragmentEnd >= extendedStarts.get(y)
							&& Integer.parseInt(gcnv.get("start", x)) <= Integer.parseInt(gcnv.get("end", y))) {

						//						System.out.println("overlap!");
						currentFragment.add(y);
						currentFragmentEnd = Math.max(currentFragmentEnd, extendedEnds.get(y));
						indexesDefragmented.add(y);
						k++;
					} else {
						//						System.out.println("nope!");
						break;
					}
				}
				i = k;

				if (currentFragment.size() != 1) {
					indexesDefragmented.add(x);
					int newStart = Integer.parseInt(gcnv.get("start", currentFragment.get(0)));
					int newEnd = 0;
					int newQA = 0;
					int newQS = 0;
					String newQSS = gcnv.get("QSS", currentFragment.get(0));
					String newQSE = gcnv.get("QSE", currentFragment.get(currentFragment.size() - 1));
					for (int n = 0; n < currentFragment.size(); n++) {
						newStart = Math.min(newStart, Integer.parseInt(gcnv.get("start", currentFragment.get(n))));
						newEnd = Math.max(newEnd, Integer.parseInt(gcnv.get("end", currentFragment.get(n))));
						newQA += Integer.parseInt(gcnv.get("QA", currentFragment.get(n)));
						newQS = Math.max(newQS, Integer.parseInt(gcnv.get("QS", currentFragment.get(n))));
					}
					newQA /= currentFragment.size();
					String[] newRow = new String[gcnv.fieldNames.length];
					for (int s = 0; s < newRow.length; s++) {
						newRow[s] = gcnv.get(x)[s];
					}
					newRow[gcnv.columnMapping.get("start")] = Integer.toString(newStart);
					newRow[gcnv.columnMapping.get("end")] = Integer.toString(newEnd);
					newRow[gcnv.columnMapping.get("QA")] = Integer.toString(newQA);
					newRow[gcnv.columnMapping.get("QS")] = Integer.toString(newQS);
					newRow[gcnv.columnMapping.get("QSS")] = newQSS;
					newRow[gcnv.columnMapping.get("QSE")] = newQSE;
					newDefragmentedCalls.add(newRow);
				} else {
					// newDefragmentedCalls.add(gcnv.get(x));
				}

			}
		}

		ArrayList<Integer> indexesDefragged = new ArrayList<>();
		indexesDefragged.addAll(indexesDefragmented);
		Collections.sort(indexesDefragged, Collections.reverseOrder());
		for (int i = 0; i < indexesDefragged.size(); i++) {
			gcnv.df.remove((int) (indexesDefragged.get(i)));
		}

		ArrayList<String> boolDefragmented_1 = new ArrayList<>();
		for (int i = 0; i < gcnv.nrow(); i++) {
			boolDefragmented_1.add("FALSE");
		}
		gcnv.addColumn("defragmented", boolDefragmented_1);

		ArrayList<String> boolDefragmented_2 = new ArrayList<>();
		for (int i = 0; i < newDefragmentedCalls.nrow(); i++) {
			boolDefragmented_2.add("TRUE");
		}
		newDefragmentedCalls.addColumn("defragmented", boolDefragmented_2);

		gcnv.rbind(newDefragmentedCalls);
		gcnv.sort();
		gcnv.writeFile(defragmentation_output, true);
	}
}
