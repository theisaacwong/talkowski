package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class JointCallVCF extends gCNVHelperTool {

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH, CLASS_COLUMN, OBSERVATION_COLUMN, COLUMN_NAME};

	public JointCallVCF(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException, InterruptedException {
		String INPUT = args.get(INPUT_PATH);
		String OUTPUT = args.get(OUTPUT_PATH);
		String classColumn = args.get(CLASS_COLUMN);
		String observationColumn = args.get(OBSERVATION_COLUMN);
		String sep = ",";
		String svtypeColumn = args.get(COLUMN_NAME);
		this.jointCallVCF(INPUT, OUTPUT, classColumn, observationColumn, sep, svtypeColumn);
	}


	/**
	 * write gcnv to joint vcf format
	 * 
	 * @param INPUT_PATH
	 * @param OUTPUT_PATH
	 * @param classColumn
	 * @param observationColumn
	 * @param sep
	 * @throws IOException
	 */
	public void jointCallVCF(String INPUT_PATH, String OUTPUT_PATH, String classColumn, String observationColumn,
			String sep, String svtypeColumn) throws IOException {
		DataFrame gcnv = new DataFrame(INPUT_PATH, true, "\t", "#");

		int svtypeColumnInt = 0;
		for (int i = 0; i < gcnv.ncol(); i++) {
			if (gcnv.fieldNames[i].equals(svtypeColumn)) {
				svtypeColumnInt = i;
			}
		}

		HashMap<String, Integer> classToFirstObservationLine = new HashMap<>();
		HashMap<String, ArrayList<String>> classToObersvations = new HashMap<>();
		HashMap<String, Boolean> classObserved = new HashMap<>();

		for (int i = 0; i < gcnv.nrow(); i++) {
			String currentClass = gcnv.get(classColumn, i);
			if (!classObserved.containsKey(currentClass)) {
				classObserved.put(currentClass, true);
				classToFirstObservationLine.put(currentClass, i);
				classToObersvations.put(currentClass, new ArrayList<>());
			}

			String observedObservation = gcnv.get(observationColumn, i);
			classToObersvations.get(currentClass).add(observedObservation);
		}

		ArrayList<String> chr = new ArrayList<>();
		ArrayList<String> start = new ArrayList<>();
		ArrayList<String> end = new ArrayList<>();
		ArrayList<String> svtype = new ArrayList<>();

		ArrayList<String> columnToAdd = new ArrayList<>();
		ArrayList<String> classes = new ArrayList<>();
		classes.addAll(classObserved.keySet());
		for (int i = 0; i < classes.size(); i++) {
			String currentClass = classes.get(i);
			columnToAdd.add(arrayToStringWithSep(classToObersvations.get(currentClass), sep));

			chr.add(gcnv.get(classToFirstObservationLine.get(currentClass), 0));
			start.add(gcnv.get(classToFirstObservationLine.get(currentClass), 1));
			end.add(gcnv.get(classToFirstObservationLine.get(currentClass), 2));
			svtype.add(gcnv.get(classToFirstObservationLine.get(currentClass), svtypeColumnInt));
		}

		ArrayList<ArrayList<String>> columnsToAdd = new ArrayList<>();
		ArrayList<String> columnNames = new ArrayList<>();

		columnsToAdd.add(chr);
		columnNames.add(gcnv.fieldNames[0]);
		columnsToAdd.add(start);
		columnNames.add(gcnv.fieldNames[1]);
		columnsToAdd.add(end);
		columnNames.add(gcnv.fieldNames[2]);
		columnsToAdd.add(classes);
		columnNames.add(classColumn);
		columnsToAdd.add(svtype);
		columnNames.add(gcnv.fieldNames[svtypeColumnInt]);
		columnsToAdd.add(columnToAdd);
		columnNames.add(observationColumn);

		DataFrame toWrite = new DataFrame();
		toWrite.addColumns(columnNames, columnsToAdd);
		toWrite.writeFile(OUTPUT_PATH, true);
	}

	
}
