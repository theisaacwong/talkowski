package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Filter extends gCNVHelperTool {

	public Filter(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException {
		String input = args[1];
		String output = args[2];
		this.filter(input, output);
	}


	/**
	 * apply quality filters, does not count X/Y chroms
	 * @param input
	 * @param output
	 * @throws IOException 
	 */
	public void filter(String input, String output) throws IOException {
		DataFrame gcnv = new DataFrame(input, true, "\\t", "@");	

		boolean[] passQS = new boolean[gcnv.nrow()];
		boolean[] passFREQ = new boolean[gcnv.nrow()];
		boolean[] lt100rawCalls = new boolean[gcnv.nrow()];
		boolean[] lt10highQSRareCalls = new boolean[gcnv.nrow()];

		HashMap<String, Integer> sampleToNRawCalls = new HashMap<>();
		for(int i = 0; i < gcnv.nrow(); i++) {
			String currSample = gcnv.get("sample", i);
			if(sampleToNRawCalls.containsKey(currSample) == false) {
				sampleToNRawCalls.put(currSample, 0);
			}
			if(!(gcnv.get("chr", i).equals("chrX") || gcnv.get("chr", i).equals("chrY")) ) {
				sampleToNRawCalls.put(currSample,
						1 + sampleToNRawCalls.get(currSample));	
			}

		}

		for(int i = 0; i < gcnv.nrow(); i++) {
			passQS[i] = (gcnv.get("svtype", i).equals("DUP") && Integer.parseInt(gcnv.get("QS", i)) >= 50) | 
					(gcnv.get("svtype", i).equals("DEL") && Integer.parseInt(gcnv.get("QS", i)) >= 100) | 
					(gcnv.get("CN", i).equals("0") && Integer.parseInt(gcnv.get("QS", i)) >= 400);
			passFREQ[i] = Double.parseDouble(gcnv.get("vaf", i)) <= 0.01;
			lt100rawCalls[i] = sampleToNRawCalls.get(gcnv.get("sample", i)) <= 100;
		}

		HashMap<String, Integer> sampleToNPseudoHQCalls = new HashMap<>();
		for(int i = 0; i < gcnv.nrow(); i++) {
			String currSample = gcnv.get("sample", i);
			if(sampleToNPseudoHQCalls.containsKey(currSample) == false) {
				sampleToNPseudoHQCalls.put(currSample, 0);
			}
			if(passQS[i] && passFREQ[i] && !(gcnv.get("chr", i).equals("chrX") || gcnv.get("chr", i).equals("chrY"))) {
				sampleToNPseudoHQCalls.put( currSample,
						1 + sampleToNPseudoHQCalls.get(currSample));	
			}

		}
		for(int i = 0; i < gcnv.nrow(); i++) {
			lt10highQSRareCalls[i] = sampleToNPseudoHQCalls.get(gcnv.get("sample", i)) <= 10;
		}

		ArrayList<String> PASS_HQ = new ArrayList<>();
		ArrayList<String> PASS_FREQ = new ArrayList<>();
		ArrayList<String> LT100_RAW_CALLS = new ArrayList<>();
		ArrayList<String> LT10_HIGH_QS_RARE_CALLS = new ArrayList<>();
		ArrayList<String> PASS_SAMPLE = new ArrayList<>();
		ArrayList<String> HIGH_QUALITY = new ArrayList<>();

		for(int i = 0; i < gcnv.nrow(); i++) {
			PASS_HQ.add(passQS[i] ? "TRUE" : "FALSE");
			PASS_FREQ.add(passFREQ[i] ? "TRUE" : "FALSE");
			LT100_RAW_CALLS.add(lt100rawCalls[i] ? "TRUE" : "FALSE");
			LT10_HIGH_QS_RARE_CALLS.add(lt10highQSRareCalls[i] ? "TRUE" : "FALSE");
			PASS_SAMPLE.add(lt100rawCalls[i] && lt10highQSRareCalls[i]? "TRUE" : "FALSE");
			HIGH_QUALITY.add(passQS[i] && passFREQ[i] && lt100rawCalls[i] && lt10highQSRareCalls[i]? "TRUE" : "FALSE");
		}

		ArrayList<String> columnNamesToAdd = new ArrayList<>();
		ArrayList<ArrayList<String>> columnValuesToAdd = new ArrayList<>();

		columnNamesToAdd.add("lt100_raw_calls");
		columnNamesToAdd.add("lt10_highQS_rare_calls");
		columnNamesToAdd.add("PASS_SAMPLE");
		columnNamesToAdd.add("PASS_FREQ");
		columnNamesToAdd.add("PASS_QS");
		columnNamesToAdd.add("HIGH_QUALITY");

		columnValuesToAdd.add(LT100_RAW_CALLS);
		columnValuesToAdd.add(LT10_HIGH_QS_RARE_CALLS);
		columnValuesToAdd.add(PASS_SAMPLE);
		columnValuesToAdd.add(PASS_FREQ);
		columnValuesToAdd.add(PASS_HQ);
		columnValuesToAdd.add(HIGH_QUALITY);

		gcnv.addColumns(columnNamesToAdd, columnValuesToAdd);
		gcnv.writeFile(output, true);

	}
	
}
