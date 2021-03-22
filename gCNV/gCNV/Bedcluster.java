package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.TreeSet;

public class Bedcluster extends gCNVHelperTool{

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH};
	
	public Bedcluster(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException, InterruptedException {
		String INPUT = args.get(INPUT_PATH);
		String OUTPUT = args.get(OUTPUT_PATH);
		
		if(args.contains("extra") == false) {
			this.bedcluster(INPUT, OUTPUT, "suffix_", false, "");
		} else {
			String prefix = args.get("extra").split(",")[0];
			String META_PATH = args.get("extra").split(",")[1];
			this.bedcluster(INPUT_PATH, OUTPUT_PATH, prefix, true, META_PATH);
		}
		
//		if(args.size() == 3) {
//			this.bedcluster(INPUT_PATH, OUTPUT_PATH, "suffix_", false, "");
//		} else if(args.size() == 4) {
//			String prefix = args[3];
//			this.bedcluster(INPUT_PATH, OUTPUT_PATH, prefix, false, "");
//		} else if(args.length == 5) {	
//			String prefix = args[3];
//			String META_PATH = args[4];
//			this.bedcluster(INPUT_PATH, OUTPUT_PATH, prefix, true, META_PATH);
//		}
	}
	

	/**
	 * 
	 * @param C2I - modified
	 * @param CC
	 * @param I2I
	 * @param CI
	 * @param indexesClustered - modified
	 * @param indexToCluster - modified
	 * @param clusterToClusterLink - modified
	 */
	public void addIndexToCluster(HashMap<String, HashSet<Integer>> C2I, String CC, HashMap<Integer, TreeSet<Integer>> I2I, int CI, HashSet<Integer> indexesClustered, HashMap<Integer, String> indexToCluster, HashMap<String, String> clusterToClusterLink, HashMap<String, TreeSet<Integer>> observedOverlaps, HashMap<Integer, ArrayList<String>> indexToObservedOverlapsForI2I ) {

		PriorityQueue<Integer> toClusterQueue = new PriorityQueue<Integer>();
		HashSet<Integer> indexesQueued = new HashSet<>();
		HashSet<Integer> currentIndexesClustered = new HashSet<>();

		toClusterQueue.add(CI);

		while(toClusterQueue.size() != 0) {
			int currentIndex = toClusterQueue.poll();


			TreeSet<Integer> I2IChildren = new TreeSet<>();
			for(Integer child : I2I.get(currentIndex)) {
				I2IChildren.add(child);
			}


			if(indexToObservedOverlapsForI2I.containsKey(currentIndex)) {
				for(String observedOverlap : indexToObservedOverlapsForI2I.get(currentIndex) ) {
					I2IChildren.addAll(observedOverlaps.get(observedOverlap));
				}	
			}


			for(Integer child : I2IChildren) {

				// check if any child indexes have already been assigned a different cluster
				if(indexToCluster.containsKey(child) && (indexToCluster.get(child).equals(CC) == false) && clusterToClusterLink.containsKey(CC)==false) {
					String correctCC = indexToCluster.get(child);
					//					print("correction!: " + CC + " -> " + correctCC);
					clusterToClusterLink.put(CC, correctCC);
				}

				if(indexesClustered.contains(child)==false && indexesQueued.contains(child)==false && currentIndexesClustered.contains(child)==false) {
					toClusterQueue.add(child);
					indexesQueued.add(child);
				}
			}



			indexToCluster.put(currentIndex, CC);
			C2I.get(CC).add(currentIndex);
			currentIndexesClustered.add(currentIndex);

		}

		indexesClustered.addAll(currentIndexesClustered);

	}

	
	public void bedcluster(String input, String output_file, String prefixString, boolean writeIndexFile, String indexFilePath) throws IOException, InterruptedException {

		//read in bed file
		DataFrame df = new DataFrame(input, true, "\\t", "@");
		df.columnMapping.put("chr", 0);
		df.fieldNames[0] = "chr";

		//sort bed file
		df.sort();

		// set overlap percent threshold
		double PERCENT_OVERLAP = 0.8;
		String variantName = "variant_name";
		boolean writeIndexMap = writeIndexFile;

		// create a map of CHR to line number of first observation
		HashMap<String, Integer> chrToFirstIndex = new HashMap<>();

		// create a map of CHR to line number of last observation
		HashMap<String, Integer> chrToLastIndex = new HashMap<>();

		// arraylist of integer starts and ends by line number (index)
		ArrayList<Integer> starts = new ArrayList<>();
		ArrayList<Integer> ends = new ArrayList<>();

		// populate starts, ends, chrToFirstIndex, chrToLastIndex
		for(int i = 0; i < df.nrow(); i++) {
			chrToLastIndex.put(df.get("chr", i), i);
			if (!chrToFirstIndex.containsKey(df.get("chr", i))) {
				chrToFirstIndex.put(df.get("chr", i), i);
			}

			starts.add(Integer.parseInt(df.get("start", i)));
			ends.add(Integer.parseInt(df.get("end", i)));
		}

		// map of end coordinate to last instance, could do the same for starts?
		HashMap<String, ArrayList<Integer>> chrToSortedEnds = new HashMap<>();

		// map of chr to map of end to first index, when the DF is sorted, some END coordinates appear multiple times, this maps the first line number that each end appears at
		HashMap<String, HashMap<Integer, Integer>> chrToEndToFirstIndex = new HashMap<>();

		for(String chr : chrToFirstIndex.keySet()) {
			HashMap<Integer, Integer> endToFirstIndex = new HashMap<>();
			int chrStart = chrToFirstIndex.get(chr);
			int chrEnd = chrToLastIndex.get(chr);
			ArrayList<Integer> rawEnds = new ArrayList<>();
			for(int i = chrEnd; i >= chrStart; i--) {
				int currEnd = ends.get(i);
				rawEnds.add(currEnd);
				endToFirstIndex.put(currEnd, i);
			}
			Collections.sort(rawEnds);
			chrToSortedEnds.put(chr, rawEnds);
			chrToEndToFirstIndex.put(chr, endToFirstIndex);
		}

		// keep track of which indexes have been checked
		HashSet<Integer> checkedIndexes = new HashSet<>();

		// map of which indexes link to which other indexes
		HashMap<Integer, TreeSet<Integer>> indexToIndexes = new HashMap<>();
		HashMap<Integer, ArrayList<String>> indexToObservedOverlapsForI2I = new HashMap<>();

		// map of which START/END coordinates have been observed so that searching only needs to be done once
		HashMap<String, Integer> visitedStarts = new HashMap<>();
		HashMap<String, Integer> visitedEnds = new HashMap<>();

		// map of which intervals are already known to intersect with a given interval
		HashMap<String, TreeSet<Integer>> observedOverlaps = new HashMap<>();

		// keep track of which CHR is currently being iterated over
		HashSet<String> newChrReset = new HashSet<>();

		// keep track of which indexes have already been clustered
		HashSet<Integer> indexesClustered = new HashSet<>();

		// map of cluster name to indexes
		HashMap<String, HashSet<Integer>> clusterToIndexes = new HashMap<>();

		String prefix = prefixString;
		int counter = 0;
		String currentCluster = prefix + counter;

		BufferedWriter output = null;
		if(writeIndexMap) {
			File file = new File(indexFilePath);
			output = new BufferedWriter(new FileWriter(file));			
		}

		// link clusters together if they have an overlapping interval
		HashMap<String, String> clusterToClusterLink = new HashMap<>();
		HashMap<Integer, String> indexToCluster = new HashMap<>();

		for(int i = 0; i < df.nrow(); i++) {
			int currStart = starts.get(i);
			int currEnd = ends.get(i);
			String currChr = df.get("chr", i);
			String currType = df.get("svtype", i);

			newChrReset.add(currChr);
			if(newChrReset.size() > 1) {
				System.out.print(currChr);

				newChrReset = new HashSet<>();
				visitedStarts.clear();
				visitedEnds.clear();


				ArrayList<Integer> indexKeys = new ArrayList<>();
				indexKeys.addAll(indexToIndexes.keySet());
				Collections.sort(indexKeys);

				for(Integer k : indexKeys){
					if(!indexesClustered.contains(k)) {
						clusterToIndexes.put(currentCluster, new HashSet<>());
						addIndexToCluster(clusterToIndexes, currentCluster, indexToIndexes, k, indexesClustered, indexToCluster, clusterToClusterLink, observedOverlaps, indexToObservedOverlapsForI2I);
						counter++;
						currentCluster = prefix + counter;
					}
				}

				if(writeIndexMap) {
					ArrayList<Integer> temp = new ArrayList<>();
					temp.addAll(indexToIndexes.keySet());
					Collections.sort(temp);
					for(int x : temp) {
						output.write(x + "\t" + indexToIndexes.get(x).toString() +  "\n");
					}					
				}

				indexToObservedOverlapsForI2I.clear();
				observedOverlaps.clear();
				indexToIndexes.clear();
				indexToCluster.clear();
			}

			indexToIndexes.put(i, new TreeSet<>());
			checkedIndexes.add(i);

			String startKey = currChr + "_" + currStart;
			int indexS;
			if(visitedStarts.containsKey(startKey)) {
				indexS = visitedStarts.get(startKey);
			} else {
				int indexOfSortedTargetEnd = BinarySearchLowerBound(chrToSortedEnds.get(currChr), currStart, -1, -1);
				int targetEnd = chrToSortedEnds.get(currChr).get(indexOfSortedTargetEnd);
				indexS = chrToEndToFirstIndex.get(currChr).get(targetEnd);
				visitedStarts.put(startKey, indexS);
			}

			String endKey = currChr + "_" + currEnd;
			int indexE;
			if(visitedEnds.containsKey(endKey)) {
				indexE = visitedEnds.get(endKey);
			} else {
				int chrStart = chrToFirstIndex.get(df.get("chr", i));
				int chrEnd = chrToLastIndex.get(df.get("chr", i));
				indexE = BinarySearchUpperBound(starts, currEnd, chrStart, chrEnd);
				visitedEnds.put(endKey, indexE);
			}

			String intervalKey = currChr +  "_" + currStart + "_" + currEnd + "_" + currType;
			if(!observedOverlaps.containsKey(intervalKey)) {
				for(int n = indexS; n <= indexE; n++) {
					if(!checkedIndexes.contains(n) && currEnd >= starts.get(n) && currStart <= ends.get(n) && currType.equals(df.get("svtype", n)) && reciprocalOverlap(currStart, currEnd, starts.get(n), ends.get(n), PERCENT_OVERLAP)) {
						indexToIndexes.get(i).add(n);
					}
				}
				observedOverlaps.put(intervalKey, new TreeSet<>());
				observedOverlaps.get(intervalKey).addAll(indexToIndexes.get(i));

			} else {
				observedOverlaps.get(intervalKey).remove(Integer.valueOf(i));
				//				indexToIndexes.get(i).addAll(observedOverlaps.get(intervalKey));
				if(indexToObservedOverlapsForI2I.containsKey(i)) {
					indexToObservedOverlapsForI2I.get(i).add(intervalKey);
				} else {
					indexToObservedOverlapsForI2I.put(i, new ArrayList<>());
					indexToObservedOverlapsForI2I.get(i).add(intervalKey);
				}
			}

		}

		visitedStarts.clear();
		visitedEnds.clear();

		ArrayList<Integer> indexKeys = new ArrayList<>();
		indexKeys.addAll(indexToIndexes.keySet());
		Collections.sort(indexKeys);
		for(Integer k : indexKeys){
			if(!indexesClustered.contains(k)) {
				clusterToIndexes.put(currentCluster, new HashSet<>());
				addIndexToCluster(clusterToIndexes, currentCluster, indexToIndexes, k, indexesClustered, indexToCluster, clusterToClusterLink, observedOverlaps, indexToObservedOverlapsForI2I);
				counter++;
				currentCluster = prefix + counter;

			}

		}
		observedOverlaps.clear();
		indexToObservedOverlapsForI2I.clear();

		if(writeIndexMap) {
			ArrayList<Integer> temp = new ArrayList<>();
			temp.addAll(indexToIndexes.keySet());
			Collections.sort(temp);
			for(int x : temp) {
				output.write(x + "\t" + indexToIndexes.get(x).toString() +  "\n");
			}
			output.close();			
		}
		indexToIndexes.clear();


		//post-process

		ArrayList<String> cluster_labels = new ArrayList<>();
		ArrayList<String> newStarts = new ArrayList<>();
		ArrayList<String> newEnds = new ArrayList<>();

		ArrayList<String> ID = new ArrayList<>();
		for(int i = 0; i < df.nrow(); i++) {
			cluster_labels.add("Error");
			newStarts.add("Error");
			newEnds.add("Error");
			ID.add(Integer.toString(i));
		}
		for(String clusterName : clusterToIndexes.keySet()) {
			String mainClusterName = clusterName;
			while(clusterToClusterLink.containsKey(mainClusterName)) {
				mainClusterName = clusterToClusterLink.get(mainClusterName);
			}

			for(int i : clusterToIndexes.get(clusterName)) {
				cluster_labels.set(i, mainClusterName);
			}
		}

		ArrayList<String> columnNamesToAdd = new ArrayList<>();
		ArrayList<ArrayList<String>> columnValuesToAdd = new ArrayList<>();

		columnNamesToAdd.add(variantName);
		columnNamesToAdd.add("ID");

		columnValuesToAdd.add(cluster_labels);
		columnValuesToAdd.add(ID);
		df.addColumns(columnNamesToAdd, columnValuesToAdd);

		// post process
		// merge the same variant within sample
		HashMap<String, HashMap<String, ArrayList<Integer>>> sampleToVariantToIndexes = new HashMap<>();
		for(int i = 0; i < df.size(); i++) {

			String sample = df.get("sample", i);
			if(!sampleToVariantToIndexes.containsKey(sample)) {
				sampleToVariantToIndexes.put(sample, new HashMap<>());
			}

			String variant = df.get(variantName, i);
			if(!sampleToVariantToIndexes.get(sample).containsKey(variant)){
				sampleToVariantToIndexes.get(sample).put(variant, new ArrayList<>());
			}

			sampleToVariantToIndexes.get(sample).get(variant).add(i);
		}


		//it looks like a lot of nested loops, however there will only be i iterations total
		HashSet<Integer> rowsToDelete = new HashSet<>();
		ArrayList<String[]> rowsToAdd = new ArrayList<>();
		for(String sample : sampleToVariantToIndexes.keySet()) {
			for(String variant : sampleToVariantToIndexes.get(sample).keySet()) {

				if(sampleToVariantToIndexes.get(sample).get(variant).size() != 1) {
					int minStart = Integer.MAX_VALUE;
					int maxEnd = Integer.MIN_VALUE;
					for(int i : sampleToVariantToIndexes.get(sample).get(variant)) {
						minStart = Math.min(minStart, starts.get(i));
						maxEnd = Math.max(maxEnd, ends.get(i));
					}
					String[] mergedRow = df.get(sampleToVariantToIndexes.get(sample).get(variant).get(0)).clone();
					mergedRow[df.columnMapping.get("start")] = Integer.toString(minStart);
					mergedRow[df.columnMapping.get("end")] = Integer.toString(maxEnd);

					rowsToAdd.add(mergedRow);
					rowsToDelete.addAll( sampleToVariantToIndexes.get(sample).get(variant));
				}
			}
		}

		ArrayList<Integer> rowsToDeleteList = new ArrayList<>();
		rowsToDeleteList.addAll(rowsToDelete);
		Collections.sort(rowsToDeleteList, Collections.reverseOrder());
		for(int i = 0; i < rowsToDeleteList.size(); i++) {
			df.df.remove((int)rowsToDeleteList.get(i));
		}
		for(String[] newRow : rowsToAdd) {
			df.df.add(newRow);
		}


		//		//calculate rmsstd, convert to median coordinates for start and end
		//		// should probably make this a method
		//		DecimalFormat format = new DecimalFormat("#.####");
		//		HashMap<String, String> variantToRmsstd = new HashMap<>();
		//		HashMap<String, ArrayList<Integer>> variantToIndexes = new HashMap<>();
		//		HashMap<String, String> variantToMedianStartCoordinate = new HashMap<>();
		//		HashMap<String, String> variantToMedianEndCoordinate = new HashMap<>();
		//		for(int i = 0; i < df.size(); i++) {
		//			String variant = df.get(variantName, i);
		//			if(!variantToIndexes.containsKey(variant)) {
		//				variantToIndexes.put(variant, new ArrayList<>());
		//			}
		//			variantToIndexes.get(variant).add(i);
		//		}
		//		for(String variant : variantToIndexes.keySet()) {
		//			int n = variantToIndexes.get(variant).size();
		//			double[] vstarts = new double[n];
		//			double[] vends = new double[n];
		//			for(int i = 0; i < n; i++) {
		//				vstarts[i] = Double.parseDouble(df.get("start", variantToIndexes.get(variant).get(i)));
		//				vends[i] = Double.parseDouble(df.get("end", variantToIndexes.get(variant).get(i)));
		//			}
		//			variantToMedianStartCoordinate.put(variant, median(vstarts));
		//			variantToMedianEndCoordinate.put(variant, median(vends));
		//			variantToRmsstd.put(variant, format.format((rmsstd(vstarts, vends))));
		//		}
		//		ArrayList<String> rmsstd = new ArrayList<>();
		//		ArrayList<String> medStarts = new ArrayList<>();
		//		ArrayList<String> medEnds = new ArrayList<>();
		//		for(int i = 0; i < df.size(); i++) {
		//			rmsstd.add(variantToRmsstd.get(df.get(variantName, i)));
		//			medStarts.add(variantToMedianStartCoordinate.get(df.get(variantName, i)));
		//			medEnds.add(variantToMedianEndCoordinate.get(df.get(variantName, i)));
		//		}
		//
		//		columnNamesToAdd = new ArrayList<>();
		//		columnValuesToAdd = new ArrayList<>();
		//		
		//		columnNamesToAdd.add("rmsstd");
		//		columnNamesToAdd.add("medStart");
		//		columnNamesToAdd.add("medEnd");
		//		
		//		columnValuesToAdd.add(rmsstd);
		//		columnValuesToAdd.add(medStarts);
		//		columnValuesToAdd.add(medEnds);
		//		
		//		df.addColumns(columnNamesToAdd, columnValuesToAdd);


		//calculate rmsstd
		DecimalFormat format = new DecimalFormat("#.####");
		HashMap<String, String> variantToRmsstd = new HashMap<>();
		HashMap<String, ArrayList<Integer>> variantToIndexes = new HashMap<>();
		for(int i = 0; i < df.size(); i++) {
			String variant = df.get(variantName, i);
			if(!variantToIndexes.containsKey(variant)) {
				variantToIndexes.put(variant, new ArrayList<>());
			}
			variantToIndexes.get(variant).add(i);
		}
		for(String variant : variantToIndexes.keySet()) {
			int n = variantToIndexes.get(variant).size();
			double[] vstarts = new double[n];
			double[] vends = new double[n];
			for(int i = 0; i < n; i++) {
				vstarts[i] = Double.parseDouble(df.get("start", variantToIndexes.get(variant).get(i)));
				vends[i] = Double.parseDouble(df.get("end", variantToIndexes.get(variant).get(i)));
			}
			variantToRmsstd.put(variant, format.format((rmsstd(vstarts, vends))));
		}
		ArrayList<String> rmsstd = new ArrayList<>();
		for(int i = 0; i < df.size(); i++) {
			rmsstd.add(variantToRmsstd.get(df.get(variantName, i)));
		}

		df.addColumn("rmsstd", rmsstd);


		df.writeFile(output_file, true);

	}
	
	public double rmsstd(double[] starts, double[] ends) {
		double SS = meanSS(starts) + meanSS(ends);
		return Math.sqrt(SS);
	}
	
	public double mean(double[] x) {
		double sum = 0;
		for(double d : x) {
			sum += d;
		}
		return sum/x.length;
	}

	public double meanSS(double[] x) {
		double mu = mean(x);

		double sumSS = 0;
		for(double d : x) {
			sumSS += Math.pow((d - mu), 2);
		}

		return sumSS/x.length;
	}
	
	public static boolean reciprocalOverlap(int s1, int e1, int s2, int e2, double percent) {
		double nOverlapx100 = Math.max(0, Math.min(e1, e2) - Math.max(s1, s2) + 1);
		double width1 = e1 - s1;
		double width2 = e2 - s2;
		return ((nOverlapx100 / width1 >= percent) && (nOverlapx100 / width2 >= percent));
	}

}
