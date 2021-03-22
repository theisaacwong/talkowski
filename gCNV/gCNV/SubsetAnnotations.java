package gCNV;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class SubsetAnnotations extends gCNVHelperTool {

	public SubsetAnnotations(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException, InterruptedException {
		String gcnvInput = args[1];
		String output_1 = args[2];
		String output_2 = args[3];
		String sourceColumnName = args[4];
		ArrayList<String> annotationSubsets = new ArrayList<>();
		for (int i = 5; i < args.length; i++) {
			annotationSubsets.add(args[i]);
		}
		this.subsetAnnotations(gcnvInput, output_1, output_2, sourceColumnName, annotationSubsets);
	}

	/**
	 * given a list of gene lists, add columns for annotations found needs columns :
	 * genes_any_overalp, genes_strict_overlap
	 * 
	 * @param gcnvInput
	 * @param output
	 * @param sourceColumnName
	 * @param annotationSubsets
	 * @throws IOException
	 */
	public void subsetAnnotations(String gcnvInput, String output_1, String output_2, String sourceColumnName, ArrayList<String> annotationSubsets) throws IOException {
		print("reading callset");
		DataFrame gcnv = new DataFrame(gcnvInput, true, "\\t", "#");

		// populate hashmap/set with list of all genes
		print("reading " + annotationSubsets.size() + " gene lists");
		HashMap<String, HashSet<String>> listNameToGenes = new HashMap<>();
		for (String geneList : annotationSubsets) {

			listNameToGenes.put(geneList, new HashSet<>());

			FileInputStream inputStream = new FileInputStream(geneList);
			Scanner sc = new Scanner(inputStream, "UTF-8");
			String line = "";
			while (sc.hasNextLine()) {
				line = sc.nextLine();
				listNameToGenes.get(geneList).add(line);
			}
			sc.close();
		}
		print("finished reading " + listNameToGenes.size() + " lists");
		print("marking annotations");

		HashMap<String, ArrayList<String>> columnNameToAnnotatedGenes = new HashMap<>();
		for (String geneListName : listNameToGenes.keySet()) {
			columnNameToAnnotatedGenes.put(geneListName, new ArrayList<>());
		}

		// for each row in gcnv callset
		for (int i = 0; i < gcnv.nrow(); i++) {

			// get array of genes found
			String[] currentGenes = gcnv.get(sourceColumnName, i).split(",");

			// for each gene list
			for (String geneListName : listNameToGenes.keySet()) {

				ArrayList<String> foundGenes = new ArrayList<>();
				// check to see if each current row gene is contained in the list
				for (String gene : currentGenes) {
					if (listNameToGenes.get(geneListName).contains(gene)) {
						foundGenes.add(gene);
					}
				}

				if (foundGenes.size() == 0) {
					foundGenes.add("None");
				}
				columnNameToAnnotatedGenes.get(geneListName).add(String.join(",", foundGenes));
			}
		}
		print("writing annotations");
		ArrayList<ArrayList<String>> columnsToAdd = new ArrayList<>();
		ArrayList<String> columnNames = new ArrayList<>();

		for (String annotationSubsetName : annotationSubsets) {
			columnsToAdd.add(columnNameToAnnotatedGenes.get(annotationSubsetName));
			String[] tempNameArray = annotationSubsetName.split("/|\\\\");
			String baseName = tempNameArray[tempNameArray.length - 1];
			columnNames.add(baseName);
		}

		gcnv.addColumns(columnNames, columnsToAdd);
		gcnv.writeFile(output_1, true);

		String OUTPUT_PATH = output_2;
		String columnToAggregateBy = "sample";
		ArrayList<String> cols = new ArrayList<>();
		cols.add(Integer.toString(gcnv.columnMapping.get("genes_any_overlap")));
		cols.add(Integer.toString(gcnv.columnMapping.get("genes_strict_overlap")));
		for (String col : columnNames) {
			cols.add(Integer.toString(gcnv.columnMapping.get(col)));
		}
		String columnsToAggregate = String.join(",", cols);
		String columnToSplitBy = "svtype";
		
		Simplify simplify = new Simplify(args);
		simplify.simplify(gcnv, OUTPUT_PATH, columnToAggregateBy, columnsToAggregate, columnToSplitBy);
	}
	
	
	

}
