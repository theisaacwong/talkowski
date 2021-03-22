package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class Simplify extends gCNVHelperTool {
	
	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH, COLUMN_TO_AGGREGATE_BY, COLUMNS_TO_AGGREGATE, COLUMN_TO_SPLIT_BY};

	public Simplify(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException, InterruptedException {
		String input = args.get(INPUT_PATH);
		String OUTPUT = args.get(OUTPUT_PATH);
		String columnToAggregateBy = args.get(COLUMN_TO_AGGREGATE_BY);
		String columnsToAggregate = args.get(COLUMNS_TO_AGGREGATE);
		String columnToSplitBy = args.get(COLUMN_TO_SPLIT_BY);
		this.simplify(input, OUTPUT, columnToAggregateBy, columnsToAggregate, columnToSplitBy);
	}
	

	public void simplify(String input, String OUTPUT_PATH, String columnToAggregateBy, String columnsToAggregate, String columnToSplitBy) throws IOException {
		DataFrame gcnv = new DataFrame(input, true, "\\t", "#");
		simplify(gcnv, OUTPUT_PATH, columnToAggregateBy, columnsToAggregate, columnToSplitBy);
	}

	/**
	 * simplify so that one per sample / aggregate by sample user can specify which
	 * columns to aggregate, eg : gene_1 user can specify conditional columns to
	 * aggregate on, eg : DUP_gene_1, DEL_gene_1 User needs to specify which columns
	 * to aggregate on, can be with
	 * 
	 * @param gcnvInput           - callset path
	 * @param columnToAggregateBy - default: "sample" , aggregate other columns by
	 *                            this column
	 * @param columnsToAggregate  - comma separated list of column numbers to
	 *                            aggregate
	 * @param columnToSplitBy     - default "svtype", split aggregation into
	 *                            categories
	 * @throws IOException
	 */
	public void simplify(DataFrame gcnv, String OUTPUT_PATH, String columnToAggregateBy, String columnsToAggregate, String columnToSplitBy) throws IOException {
		HashSet<String> splitTypes = new HashSet<>();
		for (int i = 0; i < gcnv.nrow(); i++) {
			splitTypes.add(gcnv.get(columnToSplitBy, i));
		}

		// key : sample
		// value: sample values
		// key: column name
		// value: set of value, list of genes for this column
		HashMap<String, HashMap<String, HashSet<String>>> sampleToField = new HashMap<>();

		String[] columnsToAggregateArray = columnsToAggregate.split(",");
		ArrayList<String> columnNamesToAggregate = new ArrayList<>();
		for (String s : columnsToAggregateArray) {
			columnNamesToAggregate.add(gcnv.fieldNames[Integer.parseInt(s)]);
		}

		ArrayList<String> columnsToGenerate = new ArrayList<>();
		for (String columnName : columnNamesToAggregate) {
			for (String type : splitTypes) {
				columnsToGenerate.add(columnName + "_" + type);
			}
		}

		for (int i = 0; i < gcnv.nrow(); i++) {
			String currentSample = gcnv.get(columnToAggregateBy, i);

			if (sampleToField.containsKey(currentSample) == false) {
				sampleToField.put(currentSample, new HashMap<>());
				for (String columnToGenerate : columnsToGenerate) {
					sampleToField.get(currentSample).put(columnToGenerate, new HashSet<String>());
				}
			}

			var currentFields = sampleToField.get(currentSample);

			for (String columnName : columnNamesToAggregate) {
				String columnToAddTo = columnName + "_" + gcnv.get(columnToSplitBy, i);
				String currentValue = gcnv.get(columnName, i);
				currentFields.get(columnToAddTo).add(currentValue);
			}
		}

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		String header = columnToAggregateBy + "\t" + String.join("\t", columnsToGenerate) + "\n";
		output.write(header);
		for (String sample : sampleToField.keySet()) {
			StringBuilder line = new StringBuilder();
			line.append(sample);
			line.append("\t");

			ArrayList<String> values = new ArrayList<>();
			for (String column : columnsToGenerate) {
				ArrayList<String> currentGenes = new ArrayList<>();
				currentGenes.addAll(sampleToField.get(sample).get(column));
				currentGenes.remove("None");
				Collections.sort(currentGenes);
				if (currentGenes.size() == 0) {
					currentGenes.add("None");
				}
				values.add(String.join(",", currentGenes));
			}
			line.append(String.join("\t", values));
			line.append("\n");
			output.write(line.toString());
		}
		output.close();
	}
	
	
}
