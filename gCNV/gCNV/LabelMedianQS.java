package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class LabelMedianQS extends gCNVHelperTool {

	public LabelMedianQS(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException, InterruptedException {
		String INPUT_PATH = args[1];
		String variantColumn = args[2];
		String qsColumn = args[3];
		String OUTPUT_PATH = args[4];
		this.labelMedianQS(INPUT_PATH, variantColumn, qsColumn, OUTPUT_PATH);
	}
	

	/**
	 * 
	 * @param INPUT_PATH
	 * @param variantColumn
	 * @param qsColumn
	 * @param OUTPUT_PATH
	 * @throws IOException
	 */
	public void labelMedianQS(String INPUT_PATH, String variantColumn, String qsColumn, String OUTPUT_PATH)
			throws IOException {
		System.out.println("generating metadata");
		ArrayList<DataFrame> info = getPerSampleMetrics(INPUT_PATH, variantColumn, "null", false);
		DataFrame gcnv = info.get(0);
		DataFrame metaData = info.get(1);
		System.out.println("labeling metadata");
		labelMedianQS(gcnv, variantColumn, OUTPUT_PATH, metaData, variantColumn, qsColumn + "_MEDIAN");
	}

	/**
	 * 
	 * @param INPUT_PATH
	 * @param variantColumn
	 * @param qsColumn
	 * @param OUTPUT_PATH
	 * @param metaData
	 * @throws IOException
	 */
	public void labelMedianQS(DataFrame gcnv, String variantColumn, String OUTPUT_PATH, DataFrame metaData, String metaDataVariantColmnm, String meteDataQSMedColumn) throws IOException {
		ArrayList<String> med_QS_values = new ArrayList<>();
		HashMap<String, String> nameToValue = new HashMap<>();
		for (int i = 0; i < metaData.nrow(); i++) {
			nameToValue.put(metaData.get(metaDataVariantColmnm, i), metaData.get(meteDataQSMedColumn, i));
		}
		for (int i = 0; i < gcnv.nrow(); i++) {
			med_QS_values.add(nameToValue.get(gcnv.get(variantColumn, i)));
		}

		gcnv.addColumn(meteDataQSMedColumn, med_QS_values);

		System.out.println("writing to file");
		gcnv.writeFile(OUTPUT_PATH, true);
	}
	
	

}
