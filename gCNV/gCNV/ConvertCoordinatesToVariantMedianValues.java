package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class ConvertCoordinatesToVariantMedianValues extends gCNVHelperTool  {

	public ConvertCoordinatesToVariantMedianValues(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException {
		String variantName = args[1];
		String INPUT_PATH = args[2];
		String OUTPUT_PATH = args[3];
		this.convertCoordinatesToVariantMedianValues(variantName, INPUT_PATH, OUTPUT_PATH);
	}
	
	public void convertCoordinatesToVariantMedianValues(String variantName, String input, String output) throws IOException {
		DataFrame df = new DataFrame(input, true, "\\t", "@");

		HashMap<String, ArrayList<Integer>> variantToIndexes = new HashMap<>();
		HashMap<String, String> variantToMedianStartCoordinate = new HashMap<>();
		HashMap<String, String> variantToMedianEndCoordinate = new HashMap<>();
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
			variantToMedianStartCoordinate.put(variant, median(vstarts));
			variantToMedianEndCoordinate.put(variant, median(vends));
		}
		ArrayList<String> medStarts = new ArrayList<>();
		ArrayList<String> medEnds = new ArrayList<>();
		for(int i = 0; i < df.size(); i++) {
			medStarts.add(variantToMedianStartCoordinate.get(df.get(variantName, i)));
			medEnds.add(variantToMedianEndCoordinate.get(df.get(variantName, i)));
		}

		int s = df.columnMapping.get("start");
		int e = df.columnMapping.get("end");

		for(int i = 0; i < df.nrow(); i++) {
			df.df.get(i)[s] = medStarts.get(i);
			df.df.get(i)[e] = medEnds.get(i);
		}

		df.writeFile(output, true);
	}
	
	
	
}
