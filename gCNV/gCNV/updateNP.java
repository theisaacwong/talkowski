package gCNV;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class updateNP extends gCNVHelperTool {

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH, MANIFEST_PATH};
	
	public updateNP(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException {
		String input = args.get(INPUT_PATH);
		String output = args.get(OUTPUT_PATH);
		String filteredIntervals = args.get(MANIFEST_PATH);
		this.updateNPfield(input, output, filteredIntervals);
	}


	/**
	 * update NP from filtered interval file
	 * @param input
	 * @param output
	 * @throws IOException 
	 */
	public void updateNPfield(String input, String output, String filteredIntervals) throws IOException {
		DataFrame gcnv = new DataFrame(input, true, "\\t", "@");	

		print("reading filtered intervals");
		HashMap<String, Integer> colMapsIF = new HashMap<>();
		colMapsIF.put("chr", 0);
		colMapsIF.put("start", 1);
		colMapsIF.put("end", 2);
		//Read filtered intervals, column 1 is name of cluster, column 2 is file
		DataFrame filteredIntervalsManifest = new DataFrame(filteredIntervals, true, "\\t", "#");
		HashMap<String, DataFrame> batchToFilteredIntervalDF = new HashMap<>();
		for(int i = 0; i < filteredIntervalsManifest.nrow(); i++) {
			String batchName = filteredIntervalsManifest.get(i, 0);
			String fiPath = filteredIntervalsManifest.get(i, 1);
			//check that file exists
			if(new File(fiPath).isFile()) {
				batchToFilteredIntervalDF.put(batchName, new DataFrame(fiPath, false, "\\t", "@"));
				batchToFilteredIntervalDF.get(batchName).columnMapping = colMapsIF;
			} else {
				System.out.println("WARNING: File not found:  " + fiPath);
			}
		}
		
		HashMap<String, HashMap<String, ArrayList<Integer>>> batchToChrToIndexes = new HashMap<>();
		HashMap<String, HashMap<String, ArrayList<Integer>>> batchToChrToStarts = new HashMap<>();
		for(String batch : batchToFilteredIntervalDF.keySet()) {
			batchToChrToIndexes.put(batch, new HashMap<>());
			batchToChrToStarts.put(batch, new HashMap<>());
			for(int i = 0; i < batchToFilteredIntervalDF.get(batch).nrow(); i++) {
				String currChr = batchToFilteredIntervalDF.get(batch).get("chr", i);
				if(!batchToChrToIndexes.get(batch).containsKey(currChr)) {
					batchToChrToIndexes.get(batch).put(currChr, new ArrayList<>());
					batchToChrToStarts.get(batch).put(currChr, new ArrayList<>());
				}
				batchToChrToIndexes.get(batch).get(currChr).add(i);
				batchToChrToStarts.get(batch).get(currChr).add(Integer.parseInt(batchToFilteredIntervalDF.get(batch).get("start", i)));
			}
		}
		
		print("done");
		
//		ArrayList<String> newNP = new ArrayList<>();
		for(int i = 0; i < gcnv.nrow(); i++) {
			String currBatch = gcnv.get("batch", i);
			DataFrame currIF = batchToFilteredIntervalDF.get(currBatch);
			int count = 0;
			String currChr = gcnv.get("chr", i);
			int currStart = Integer.parseInt(gcnv.get("start", i));
			int currEnd = Integer.parseInt(gcnv.get("end", i));
			
			ArrayList<Integer> currStarts = batchToChrToStarts.get(currBatch).get(currChr);
			ArrayList<Integer> currIndexes = batchToChrToIndexes.get(currBatch).get(currChr);
			int index = Math.max(Collections.binarySearch(currStarts, currStart) - 10, 0);
			
			for(int j = index; j < currIndexes.size(); j++) {
				int k = currIndexes.get(j);
				boolean boolStart = Integer.parseInt(currIF.get("start", k)) <= currEnd;
				boolean boolEnd = Integer.parseInt(currIF.get("end", k)) >= currStart;
				if(boolStart && boolEnd) {
					count++;
				}
			}
			
//			newNP.add(Integer.toString(count));
			gcnv.df.get(i)[gcnv.columnMapping.get("NP")] = Integer.toString(count);
		}
		
//		gcnv.addColumn("newNP", newNP);
		gcnv.writeFile(output, true);
		
	}
	
}
