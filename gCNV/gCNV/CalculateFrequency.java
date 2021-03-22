package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class CalculateFrequency extends gCNVHelperTool{

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH, COLUMN_NAME};
	
	public CalculateFrequency(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException {
		String INPUT = args.get(INPUT_PATH);
		String variantColumn = args.get(COLUMN_NAME);
		String OUTPUT = args.get(OUTPUT_PATH);
		this.calculateFrequency(INPUT, variantColumn, OUTPUT);
	}
	
	
	public void calculateFrequency(String input, String variantColumn, String output) throws IOException {

		DataFrame df = new DataFrame(input, true, "\\t", "@");
		//		df.columnMapping.put("chr", 0);
		//		df.fieldNames[0] = "chr";
		df.sort();

		HashMap<String, Double> variantToCount = new HashMap<>();
		HashSet<String> samples = new HashSet<>();

		for(int i = 0; i < df.size(); i++) {
			String variant = df.get(variantColumn, i);
			if(!variantToCount.containsKey(variant)) {
				variantToCount.put(variant, 0.0);
			}
			variantToCount.put(variant, variantToCount.get(variant)+1);
			samples.add(df.get("sample", i));
		}

		double nSamples = samples.size();
		//DecimalFormat format = new DecimalFormat("#.####");

		HashMap<String, String> variantToVAF = new HashMap<>();
		for(String variant : variantToCount.keySet()) {
			//variantToVAF.put(variant, format.format(variantToCount.get(variant)/nSamples));
			variantToVAF.put(variant, Double.toString(variantToCount.get(variant)/nSamples));
		}

		ArrayList<String> vaf = new ArrayList<>();
		ArrayList<String> vac = new ArrayList<>();
		for(int i = 0; i < df.size(); i++) {
			String variant = df.get(variantColumn, i);
			vaf.add(variantToVAF.get(variant));
			vac.add(Double.toString(variantToCount.get(variant)));
		}


		ArrayList<String> columnNamesToAdd = new ArrayList<>();
		ArrayList<ArrayList<String>> columnValuesToAdd = new ArrayList<>();

		columnNamesToAdd.add("vaf");
		columnNamesToAdd.add("vac");

		columnValuesToAdd.add(vaf);
		columnValuesToAdd.add(vac);
		df.addColumns(columnNamesToAdd, columnValuesToAdd);
		df.writeFile(output, true);
	}

	
}
