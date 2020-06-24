package gCNV;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

/**
 * 
 * @author Isaac Wong
 * @since December 2, 2019 Massachusetts General Hospital Center for Genomic
 *        Medicine Talkowski Lab
 * 
 *        This tool is designed as a bridge between the various steps of gCNV
 * 
 *        Contact Information website: https://talkowski.mgh.harvard.edu/ email:
 *        iwong@broadinstitute.org github: https://github.com/theisaacwong/
 *
 */
public class gCNV_helper {

	public String[] initializationArgs;
	public String date;

	public static long startTime;
	public static long endTime;

	public final String CHR = "CHR";
	public final String START = "START";
	public final String END = "END";
	public static final String VERSION = "2.18";

	public gCNV_helper(String[] args) {
		initializationArgs = args;
		date = Calendar.getInstance().getTime().toString();
	}

	public static void main(String[] args) throws IOException, InterruptedException {

		gCNV_helper g = new gCNV_helper(args);
		System.out.println(g.toString());
		System.out.println("version " + VERSION);
		try {
			g.checkVersion();
		} catch (Exception e) {
		}
		System.out.println("Java version: " + System.getProperty("java.version"));
		System.out.println("Heap Size: " + getHeapSize());
		start();
		g.runBetter(args);
		stop();
	}

	public void convertToEnsemble(String GCNV_INPUT, String OUTPUT, String geneColumnName, String gencodeGTF) throws IOException {
		print("reading annotation file");

		HashMap<String, String> nameToEnsembleID = new HashMap<>();

		Pattern geneNamePattern = Pattern.compile("(?<=gene_name \").+?(?=\")");
		Pattern ensemblePattern = Pattern.compile("(?<=gene_id \").+?(?=\")");

		FileInputStream inputStream = new FileInputStream(gencodeGTF);
		Scanner gtf = new Scanner(inputStream, "UTF-8");
		while (gtf.hasNext()) {
			String line = gtf.nextLine();
			if (line.startsWith("#"))
				continue;

			Matcher geneNameMatcher = geneNamePattern.matcher(line);
			String geneName = geneNameMatcher.find() ? geneNameMatcher.group() : "-1";

			Matcher ensembleMatcher = ensemblePattern.matcher(line);
			String ensemble = ensembleMatcher.find() ? ensembleMatcher.group() : "-1";

			nameToEnsembleID.put(geneName, ensemble);

		}
		gtf.close();

		nameToEnsembleID.put("None", "None");

		DataFrame gcnv = new DataFrame(GCNV_INPUT, true, "\\t", "#");
		ArrayList<String> ensembleIDcolummn = new ArrayList<>();
		for (int i = 0; i < gcnv.nrow(); i++) {
			ArrayList<String> ensembleIDs = new ArrayList<>();
			String[] geneNames = gcnv.get(geneColumnName, i).split(",");
			for (String geneName : geneNames) {
				ensembleIDs.add(nameToEnsembleID.get(geneName));
			}
			ensembleIDcolummn.add(String.join(",", ensembleIDs));
		}

		gcnv.addColumn(geneColumnName + "_Ensemble_ID", ensembleIDcolummn);
		gcnv.writeFile(OUTPUT, true);
	}

	public void countExons(String GCNV_INPUT, String genesColumnName, String gencodeGTF, String output) throws IOException {
		ArrayList<Gene> gtfGenes = parseGTFFile(gencodeGTF);
		HashMap<String, Gene> geneNameToGene = new HashMap<>();
		for (int i = 0; i < gtfGenes.size(); i++) {
			geneNameToGene.put(gtfGenes.get(i).name, gtfGenes.get(i));
		}

		DataFrame gcnv = new DataFrame(GCNV_INPUT, true, "\\t", "#");
		ArrayList<String> nExonsColumn = new ArrayList<>();
		ArrayList<String> totalExons = new ArrayList<>();
		for (int i = 0; i < gcnv.nrow(); i++) {
			String[] genes = gcnv.get(genesColumnName, i).split(",");
			if (genes[0].equalsIgnoreCase("None")) {
				nExonsColumn.add("NA");
				totalExons.add("NA");
				continue;
			}
			ArrayList<String> nExons = new ArrayList<>();
			int gStart = Integer.parseInt(gcnv.get("start", i));
			int gEnd = Integer.parseInt(gcnv.get("end", i));
			int N = 0;
			for (String gene : genes) {

				int exons = geneNameToGene.containsKey(gene) ? geneNameToGene.get(gene).getNExons(gStart, gEnd) : 0;
				nExons.add(Integer.toString(exons));
				N += exons;

			}
			totalExons.add(Integer.toString(N));
			nExonsColumn.add(String.join(",", nExons));
		}

		ArrayList<String> columnNamesToAdd = new ArrayList<>();
		columnNamesToAdd.add(genesColumnName + "_exonsPerGene");
		columnNamesToAdd.add(genesColumnName + "_totalExons");

		ArrayList<ArrayList<String>> columnValuesToAddArrayList = new ArrayList<>();
		columnValuesToAddArrayList.add(nExonsColumn);
		columnValuesToAddArrayList.add(totalExons);

		gcnv.addColumns(columnNamesToAdd, columnValuesToAddArrayList);

		gcnv.writeFile(output, true);
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

	/**
	 * given a list of gene lists, see if the source annotation gtf contains them,
	 * for checking purposes
	 * 
	 * @param gtfFile
	 * @param annotationSubsets
	 * @throws FileNotFoundException
	 */
	public void validateSubsetAnnotations(String gtfFile, ArrayList<String> annotationSubsets) throws FileNotFoundException {
		ArrayList<Gene> gtfGenes = parseGTFFile(gtfFile);
		HashSet<String> geneNames = new HashSet<>();
		for (Gene gene : gtfGenes) {
			geneNames.add(gene.name);
		}

		System.out.println();

		for (String s : annotationSubsets) {

			FileInputStream inputStream = new FileInputStream(s);
			Scanner sc = new Scanner(inputStream, "UTF-8");
			String line = "";
			System.out.print(s + ":\t");
			boolean pass = true;
			while (sc.hasNextLine()) {
				line = sc.nextLine();

				if (geneNames.contains(line) == false) {
					if (pass)
						System.out.print("FAIL\n");
					System.out.println(line);
					pass = false;
				}
			}

			if (pass) {
				System.out.print("PASS\n");
			}
			sc.close();
		}
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
		simplify(gcnv, OUTPUT_PATH, columnToAggregateBy, columnsToAggregate, columnToSplitBy);
	}

	/**
	 * TODO: figure out if the filtered intervals are needed or not
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
		// for each batch
		// batchLoop:
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

	/**
	 * reads in a GTF format file,
	 * 
	 * @param INPUT_GTF
	 * @throws FileNotFoundException
	 */
	public static ArrayList<Gene> parseGTFFile(String INPUT_GTF) throws FileNotFoundException {
		print("reading annotation file");

		HashMap<String, ArrayList<Integer>> geneToStart = new HashMap<>();
		HashMap<String, ArrayList<Integer>> geneToEnd = new HashMap<>();
		HashMap<String, String> geneToChr = new HashMap<>();

		int gtfChr = 0;
		int gtfFeatureType = 2;
		int gtfStart = 3;
		int gtfEnd = 4;

		Pattern geneNamePattern = Pattern.compile("(?<=gene_name \").+?(?=\")");
		Pattern geneTypePattern = Pattern.compile("(?<=gene_type \").+?(?=\")");

		FileInputStream inputStream = new FileInputStream(INPUT_GTF);
		Scanner gtf = new Scanner(inputStream, "UTF-8");
		while (gtf.hasNext()) {
			String line = gtf.nextLine();
			if (line.startsWith("#"))
				continue;
			String[] linee = line.split("\\t");

			Matcher geneTypeMatcher = geneTypePattern.matcher(line);
			String geneType = geneTypeMatcher.find() ? geneTypeMatcher.group() : "-1";
			String featureType = linee[gtfFeatureType];

			if (geneType.equals("protein_coding") && featureType.equals("exon")) {

				Matcher geneNameMatcher = geneNamePattern.matcher(line);
				String geneName = geneNameMatcher.find() ? geneNameMatcher.group() : "-1";

				if (geneToStart.containsKey(geneName) == false) {
					geneToStart.put(geneName, new ArrayList<>());
					geneToEnd.put(geneName, new ArrayList<>());
					geneToChr.put(geneName, linee[gtfChr]);
				}
				geneToStart.get(geneName).add(Integer.parseInt(linee[gtfStart]));
				geneToEnd.get(geneName).add(Integer.parseInt(linee[gtfEnd]));
			}
		}
		gtf.close();
		gtf = null;

		print("analyzing annotations");

		ArrayList<Gene> genes = new ArrayList<>();
		for (String gene : geneToChr.keySet()) {

			HashMap<Integer, Integer> startIdToIndex = new HashMap<>();
			var currentStarts = geneToStart.get(gene);
			for (int i = 0; i < currentStarts.size(); i++) {
				if (!startIdToIndex.containsKey(currentStarts.get(i))) {
					startIdToIndex.put(currentStarts.get(i), i);
				} else {
					int oldEnd = geneToEnd.get(gene).get(startIdToIndex.get(currentStarts.get(i)));
					int newEnd = geneToEnd.get(gene).get(i);
					if (newEnd > oldEnd) {
						startIdToIndex.put(currentStarts.get(i), i);
					} // else do nothing
				}
			}

			Collections.sort(currentStarts);
			ArrayList<Integer> sortedEnds = new ArrayList<>();
			var currentEnds = geneToEnd.get(gene);
			for (int i = 0; i < currentStarts.size(); i++) {
				sortedEnds.add(currentEnds.get(startIdToIndex.get(currentStarts.get(i))));

			}

			genes.add(new Gene(gene, geneToChr.get(gene), currentStarts, sortedEnds));
		}

		Collections.sort(genes);
		print("read " + genes.size() + " genes");
		return genes;
	}

	public void annotateGenesPercentBased(ArrayList<Gene> genes, String GCNV_INPUT, String OUTPUT_PATH, String typeColumnName) throws IOException {
		ArrayList<Integer> annoStarts = new ArrayList<>();
		ArrayList<Integer> annoEnds = new ArrayList<>();
		HashMap<String, Integer> chrToIndex = new HashMap<>();
		HashMap<String, String> variantToGene = new HashMap<>();
		HashMap<String, Integer> chrToLastIndex = new HashMap<>();

		for (int i = 0; i < genes.size(); i++) {
			annoStarts.add(genes.get(i).minStart);
			annoEnds.add(genes.get(i).maxEnd);
			chrToLastIndex.put(genes.get(i).chr, i);
			if (!chrToIndex.containsKey(genes.get(i).chr)) {
				chrToIndex.put(genes.get(i).chr, i);
			}
		}

		System.out.println("writing intersections");
		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = new FileInputStream(GCNV_INPUT);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		//String[] gColumnMapping = gColumnString.split(",");
		int gcnvChr = 0; //Integer.parseInt(gColumnMapping[0]) - 1;
		int gcnvStart = 1; //Integer.parseInt(gColumnMapping[1]) - 1;
		int gcnvEnd = 2; //Integer.parseInt(gColumnMapping[2]) - 1;
		int gcnvType = 4; //Integer.parseInt(gColumnMapping[3]) - 1;

		line = sc.nextLine();
		if (line.split("\\t")[gcnvChr].equals("chr") || line.split("\\t")[gcnvChr].equals("CHROM")) {
			output.write(line + "\tgenes_any_overlap" + "\n");
			String[] linee = line.split("\\t");
			for(int i = 0; i < linee.length; i++) {
				if(linee[i].equals(typeColumnName)) {
					gcnvType = i; break;
				}
			}
		} else {
			sc.close();
			sc = null;
			sc = new Scanner(inputStream, "UTF-8");
		}
		
		while (sc.hasNextLine()) {
			line = sc.nextLine();
			String[] linee = line.split("\\t");

			String varName = linee[gcnvChr] + "_" + linee[gcnvStart] + "_" + linee[gcnvEnd] + "_" + linee[gcnvType];
			if (variantToGene.containsKey(varName)) {
				String overlappingGenes = variantToGene.get(varName);
				output.write(line + "\t" + overlappingGenes + "\n");
			} else {

				String gChr = linee[gcnvChr];
				int gStart = Integer.parseInt(linee[gcnvStart]);
				int gEnd = Integer.parseInt(linee[gcnvEnd]);

				if (!chrToIndex.containsKey(gChr)) {
					continue;
				}

				HashSet<String> annos = new HashSet<>();
				int firstAnnoEndGTgStart = chrToIndex.get(gChr);
				int lastAnnoStartGTgEnd = chrToLastIndex.get(gChr);
				while ((firstAnnoEndGTgStart < genes.size()) && gStart > annoEnds.get(firstAnnoEndGTgStart)
						&& gChr.equals(genes.get(firstAnnoEndGTgStart).chr)) {
					firstAnnoEndGTgStart++;
				}
				while ((lastAnnoStartGTgEnd > 0) && gEnd < annoStarts.get(lastAnnoStartGTgEnd)
						&& gChr.equals(genes.get(lastAnnoStartGTgEnd).chr)) {
					lastAnnoStartGTgEnd--;
				}

				int buffer = 3;
				firstAnnoEndGTgStart = Math.max(0, firstAnnoEndGTgStart - buffer);
				lastAnnoStartGTgEnd = Math.min(genes.size() - 1, lastAnnoStartGTgEnd + buffer);
				// System.out.println("n overlaps: " + (lastAnnoStartGTgEnd -
				// firstAnnoEndGTgStart));
				for (int i = firstAnnoEndGTgStart; i <= lastAnnoStartGTgEnd; i++) {

					double percentOverlap = genes.get(i).calculateOverlapPercent(gStart, gEnd);
					// System.out.println("%: " + percentOverlap);
					if (linee[gcnvType].equals("DEL") && percentOverlap >= 0.1 && gChr.equals(genes.get(i).chr)) {
						annos.add(genes.get(i).name);
					} else if (linee[gcnvType].equals("DUP") && percentOverlap >= 0.75 && gChr.equals(genes.get(i).chr)) {
						annos.add(genes.get(i).name);
					}

				}

				ArrayList<String> annosal = new ArrayList<>();
				annosal.addAll(annos);
				Collections.sort(annosal);

				String overlappingGenes = String.join(",", annosal);
				if (overlappingGenes.equals("")) {
					overlappingGenes = "None";
				}
				variantToGene.put(varName, overlappingGenes);
				output.write(line + "\t" + overlappingGenes + "\n");

			}

		}
		output.close();
		sc.close();
	}

	/**
	 * is string builder creation slower or faster than multiple string.equals()
	 * calls? I'm currently thinking that string.equals() calls might be faster?
	 * I'll probably do both just to see
	 * 
	 * @param INPUT_PATH
	 * @param OUTPUT_PATH
	 * @param columnsToHashOnString
	 * @param columnsToMergeString
	 * @throws IOException
	 */
	public void condenseBedtoolsIntersect(String INPUT_PATH, String OUTPUT_PATH, String columnsToHashOnString,
			String columnsToMergeString, String columnsToKeepString) throws IOException {
		int[] columnsToHash = Arrays.asList(columnsToHashOnString.split(",")).stream().mapToInt(Integer::parseInt)
				.toArray();
		int[] columnsToMerge = Arrays.asList(columnsToMergeString.split(",")).stream().mapToInt(Integer::parseInt)
				.toArray();
		int[] columnsToKeep = Arrays.asList(columnsToKeepString.split(",")).stream().mapToInt(Integer::parseInt)
				.toArray();
		String tab = "\t";
		String newline = "\n";
		int lastColumnToMerge = columnsToMerge[columnsToMerge.length - 1];

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = new FileInputStream(INPUT_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		Map<Integer, Set<String>> columnToMergedStrings = new HashMap<>();
		for (int i : columnsToMerge) {
			columnToMergedStrings.put(i, new HashSet<>());
		}

		String[] previousLinee = sc.nextLine().split("\t"); // previous hash?
		int n_merged = 1;

		while (sc.hasNextLine()) {
			String currentLine = "";
			try {
				currentLine = sc.nextLine();
			} catch (Exception e) {
				System.out.println(e);
			}

			boolean isSame = true;
			String[] currentLinee = currentLine.split("\t");

			// not sure if this is better or worse than comparing concatenated strings to
			// each other
			for (int i : columnsToHash) {
				if (!previousLinee[i].equals(currentLinee[i])) {
					isSame = false;
					break;
				}
			}

			if (isSame) {
				for (int i : columnsToMerge) {
					columnToMergedStrings.get(i).add(currentLinee[i]);
				}
				n_merged++;
			} else {
				StringBuilder lineToWrite = new StringBuilder();
				for (int i : columnsToKeep) {
					lineToWrite.append(previousLinee[i]);
					lineToWrite.append(tab);
				}
				for (int i : columnsToMerge) {
					List<String> values = new ArrayList<>();
					values.addAll(columnToMergedStrings.get(i));
					Collections.sort(values);
					if (values.size() == 0) {
						values.add("None");
					}
					lineToWrite.append(String.join(",", values));
					if (i != lastColumnToMerge) {
						lineToWrite.append(tab);
					}
				}

				lineToWrite.append(tab);
				lineToWrite.append(n_merged);
				lineToWrite.append(newline);

				n_merged = 1;

				output.write(lineToWrite.toString());

				previousLinee = currentLinee;
				for (int i : columnsToMerge) {
					columnToMergedStrings.put(i, new HashSet<>());
				}
			}
		}

		output.close();
		sc.close();
	}

	public void annotateGenesAnyOverlap(ArrayList<Gene> genes, String GCNV_INPUT, String OUTPUT_PATH, String typeColumnName) throws IOException {
		ArrayList<Integer> annoStarts = new ArrayList<>();
		ArrayList<Integer> annoEnds = new ArrayList<>();
		HashMap<String, Integer> chrToIndex = new HashMap<>();
		HashMap<String, String> variantToGene = new HashMap<>();
		HashMap<String, Integer> chrToLastIndex = new HashMap<>();

		for (int i = 0; i < genes.size(); i++) {
			annoStarts.add(genes.get(i).minStart);
			annoEnds.add(genes.get(i).maxEnd);
			chrToLastIndex.put(genes.get(i).chr, i);
			if (!chrToIndex.containsKey(genes.get(i).chr)) {
				chrToIndex.put(genes.get(i).chr, i);
			}
		}

		System.out.println("writing intersections");
		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = new FileInputStream(GCNV_INPUT);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		//String[] gColumnMapping = gColumnString.split(",");
		int gcnvChr = 0; //Integer.parseInt(gColumnMapping[0]) - 1;
		int gcnvStart = 1; //Integer.parseInt(gColumnMapping[1]) - 1;
		int gcnvEnd = 2; //Integer.parseInt(gColumnMapping[2]) - 1;
		int gcnvType = 4; //Integer.parseInt(gColumnMapping[3]) - 1;

		line = sc.nextLine();
		if (line.split("\\t")[gcnvChr].equals("chr") || line.split("\\t")[gcnvChr].equals("CHROM")) {
			output.write(line + "\tgenes_any_overlap" + "\n");
			String[] linee = line.split("\\t");
			for(int i = 0; i < linee.length; i++) {
				if(linee[i].equals(typeColumnName)) {
					gcnvType = i; break;
				}
			}
		} else {
			sc.close();
			sc = null;
			sc = new Scanner(inputStream, "UTF-8");
		}

		while (sc.hasNextLine()) {
			line = sc.nextLine();
			String[] linee = line.split("\\t");

			String varName = linee[gcnvChr] + "_" + linee[gcnvStart] + "_" + linee[gcnvEnd] + "_" + linee[gcnvType];
			if (variantToGene.containsKey(varName)) {
				String overlappingGenes = variantToGene.get(varName);
				output.write(line + "\t" + overlappingGenes + "\n");
			} else {

				String gChr = linee[gcnvChr];
				int gStart = Integer.parseInt(linee[gcnvStart]);
				int gEnd = Integer.parseInt(linee[gcnvEnd]);

				if (!chrToIndex.containsKey(gChr)) {
					continue;
				}

				HashSet<String> annos = new HashSet<>();
				int firstAnnoEndGTgStart = chrToIndex.get(gChr);
				int lastAnnoStartGTgEnd = chrToLastIndex.get(gChr);
				while ((firstAnnoEndGTgStart < genes.size()) && gStart > annoEnds.get(firstAnnoEndGTgStart)
						&& gChr.equals(genes.get(firstAnnoEndGTgStart).chr)) {
					firstAnnoEndGTgStart++;
				}
				while ((lastAnnoStartGTgEnd > 0) && gEnd < annoStarts.get(lastAnnoStartGTgEnd)
						&& gChr.equals(genes.get(lastAnnoStartGTgEnd).chr)) {
					lastAnnoStartGTgEnd--;
				}

				int buffer = 3;
				firstAnnoEndGTgStart = Math.max(0, firstAnnoEndGTgStart - buffer);
				lastAnnoStartGTgEnd = Math.min(genes.size() - 1, lastAnnoStartGTgEnd + buffer);
				for (int i = firstAnnoEndGTgStart; i <= lastAnnoStartGTgEnd; i++) {

					if (genes.get(i).calculateOverlapPercent(gStart, gEnd) > 0 && gChr.equals(genes.get(i).chr)) {
						annos.add(genes.get(i).name);
					}

				}

				ArrayList<String> annosal = new ArrayList<>();
				annosal.addAll(annos);
				Collections.sort(annosal);

				String overlappingGenes = String.join(",", annosal);
				if (overlappingGenes.equals("")) {
					overlappingGenes = "None";
				}
				variantToGene.put(varName, overlappingGenes);
				output.write(line + "\t" + overlappingGenes + "\n");

			}

		}
		output.close();
		sc.close();
	}


	public void jointCallVCF(String INPUT_PATH, String OUTPUT_PATH, String classColumn, String observationColumn)
			throws IOException {
		String defaultSep = ",";
		String defaultSvtypeCol = "svtype";
		this.jointCallVCF(INPUT_PATH, OUTPUT_PATH, classColumn, observationColumn, defaultSep, defaultSvtypeCol);
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

	/**
	 * returns the mean and median of the integer columns for each unique ID of
	 * given column
	 * 
	 * @param input
	 * @param sampleColumn
	 * @throws IOException
	 */
	public ArrayList<DataFrame> getPerSampleMetrics(String INPUT_PATH, String fieldColumn, String OUTPUT_PATH,
			boolean writeFile) throws IOException {
		DataFrame gcnv = new DataFrame(INPUT_PATH, true, "\\t", "#");
		System.out.println(gcnv.nrow() + " rows");
		System.out.println(gcnv.ncol() + " columns");
		System.out.println(gcnv.columnMapping.toString());

		// key: sample name
		// value: hashmap of column names and values
		HashMap<String, HashMap<String, ArrayList<Integer>>> perFieldMap = new HashMap<>();
		HashMap<String, HashMap<String, Integer>> perFieldCalculations = new HashMap<>();

		// data needs at least one entry or will die
		String EXAMPLE_VALUE = gcnv.get(fieldColumn, 0);

		// get a list of all the integer columns
		// note: this assumes a populated field
		ArrayList<String> columnsToMeasure = new ArrayList<>();
		for (String field : gcnv.fieldNames) {
			if (field.equalsIgnoreCase("chr") || field.equalsIgnoreCase("start") || field.equalsIgnoreCase("end")) {
				continue;
			}
			if (isInteger(gcnv.get(field, 0))) {
				columnsToMeasure.add(field);
			}
		}

		// populate hashmap
		for (int i = 0; i < gcnv.nrow(); i++) {
			String currField = gcnv.get(fieldColumn, i); // current field looked at
			if (!perFieldMap.containsKey(currField)) { // if the field value has not been observed yet
				perFieldMap.put(currField, new HashMap<>()); // add the field observation and a new hashmap to the main
																// hashmap
				for (String field : columnsToMeasure) { // populate that value's hashmap with empty arraylists
					perFieldMap.get(currField).put(field, new ArrayList<>());
				}
			}

			// add the value to the appropriate hashmap value
			for (String field : columnsToMeasure) {
				try {
					perFieldMap.get(currField).get(field).add(Integer.parseInt(gcnv.get(field, i)));
				} catch (NumberFormatException nfe) {
					try {
						perFieldMap.get(currField).get(field).add((int) Double.parseDouble(gcnv.get(field, i)));
					} catch (NumberFormatException nfe2) {
						perFieldMap.get(currField).get(field).add(-1);
					}
				}
			}
		}

		// do the calculations
		for (String uniqueValue : perFieldMap.keySet()) {
			perFieldCalculations.put(uniqueValue, new HashMap<>());
			for (String field : columnsToMeasure) {
				Integer mean = mean(perFieldMap.get(uniqueValue).get(field));
				Integer median = median(perFieldMap.get(uniqueValue).get(field));

				perFieldCalculations.get(uniqueValue).put(field + "_MEAN", mean);
				perFieldCalculations.get(uniqueValue).put(field + "_MEDIAN", median);
			}
		}

		// populate data frame values
		ArrayList<ArrayList<String>> values = new ArrayList<>();
		ArrayList<String> columnNames = new ArrayList<>();
		columnNames.addAll(perFieldCalculations.get(EXAMPLE_VALUE).keySet());
		for (int i = 0; i < columnNames.size(); i++) {
			values.add(new ArrayList<>());
		}
		ArrayList<String> nameColumn = new ArrayList<>();
		for (String uniqueValue : perFieldMap.keySet()) {
			nameColumn.add(uniqueValue);
			for (int i = 0; i < columnNames.size(); i++) {
				values.get(i).add(String.valueOf(perFieldCalculations.get(uniqueValue).get(columnNames.get(i))));
			}
		}

		values.add(nameColumn);
		columnNames.add(fieldColumn);
		DataFrame toWrite = new DataFrame();
		toWrite.addColumns(columnNames, values);

		if (writeFile == true) {
			toWrite.writeFile(OUTPUT_PATH, true);
		}

		ArrayList<DataFrame> rList = new ArrayList<>();
		rList.add(gcnv);
		rList.add(toWrite);
		return rList;
	}

	/**
	 * determine if a string is an integer using Java's Integer parse method while
	 * this is not the most efficient way to check for integer-ness, I like that it
	 * handles some exceptions for me, as I will be calling parseInt() eventually
	 * 
	 * @param input
	 * @return
	 */
	public static boolean isInteger(String input) {
		try {
			Integer.parseInt(input);
			return true;
		} catch (NumberFormatException nfe) {
			return false;
		}
	}

	public String toString() {
		return String.join(" ", initializationArgs) + "\n" + date;
	}

	// TODO: MAKE GENERIC
	/**
	 * currently, I am hard-coding things, but in the future I can make the maps
	 * dynamically
	 * 
	 * @param svtk_input
	 * @param svtk_output
	 * @param match_output
	 * @throws IOException
	 */
	public void svtkMatch(String svtk_input, String svtk_output, String match_output) throws IOException {
		// read in the data frames for the svtk input and output
		DataFrame svtkInput = new DataFrame(svtk_input, true, "\\t", "@");
		DataFrame svtkOutput = new DataFrame(svtk_output, true, "\\t", "@");

		// create hashmaps to map 'name' to CN/QS/etc values
		HashMap<String, String> CN_map = new HashMap<>();
		HashMap<String, String> GT_map = new HashMap<>();
		HashMap<String, String> NP_map = new HashMap<>();
		HashMap<String, String> QA_map = new HashMap<>();
		HashMap<String, String> QS_map = new HashMap<>();
		HashMap<String, String> QSE_map = new HashMap<>();
		HashMap<String, String> QSS_map = new HashMap<>();
		HashMap<String, String> ploidy_map = new HashMap<>();
		HashMap<String, String> strand_map = new HashMap<>();
		for (int i = 0; i < svtkInput.nrow(); i++) {
			CN_map.put(svtkInput.get("name", i), svtkInput.get("CN", i));
			GT_map.put(svtkInput.get("name", i), svtkInput.get("GT", i));
			NP_map.put(svtkInput.get("name", i), svtkInput.get("NP", i));
			QA_map.put(svtkInput.get("name", i), svtkInput.get("QA", i));
			QS_map.put(svtkInput.get("name", i), svtkInput.get("QS", i));
			QSE_map.put(svtkInput.get("name", i), svtkInput.get("QSE", i));
			QSS_map.put(svtkInput.get("name", i), svtkInput.get("QSS", i));
			ploidy_map.put(svtkInput.get("name", i), svtkInput.get("ploidy", i));
			strand_map.put(svtkInput.get("name", i), svtkInput.get("strand", i));
		}

		int oldnrow = svtkOutput.nrow();
		int counter = 0;
		// svtk will sometimes have two callnames in on index, so duplicate the line for
		// each element in the bin
		for (int i = 0; i < svtkOutput.nrow(); i++) {
			if (svtkOutput.get("call_name", i).contains(",")) {
				String[] call_names = svtkOutput.get("call_name", i).split(",");
				svtkOutput.get(i)[svtkOutput.columnMapping.get("call_name")] = call_names[0];
				for (int k = 1; k < call_names.length; k++) {
					String[] newRow = new String[svtkOutput.get(i).length];
					for (int j = 0; j < newRow.length; j++) {
						newRow[j] = svtkOutput.get(i, j);
					}
					newRow[svtkOutput.columnMapping.get("call_name")] = call_names[k];
					svtkOutput.df.add(newRow); // shouldn't need to i--
					counter++;
				}
			}
		}
		System.out.println("oldnrow: " + oldnrow + "\tnewnrow: " + svtkOutput.nrow() + "\trowsadded: " + counter);

		// create arraylists for all the new fields to be added
		ArrayList<String> CN_newField = new ArrayList<>();
		ArrayList<String> GT_newField = new ArrayList<>();
		ArrayList<String> NP_newField = new ArrayList<>();
		ArrayList<String> QA_newField = new ArrayList<>();
		ArrayList<String> QS_newField = new ArrayList<>();
		ArrayList<String> QSE_newField = new ArrayList<>();
		ArrayList<String> QSS_newField = new ArrayList<>();
		ArrayList<String> ploidy_newField = new ArrayList<>();
		ArrayList<String> strand_newField = new ArrayList<>();

		for (int i = 0; i < svtkOutput.nrow(); i++) {
			CN_newField.add(CN_map.get(svtkOutput.get("call_name", i)));
			GT_newField.add(GT_map.get(svtkOutput.get("call_name", i)));
			NP_newField.add(NP_map.get(svtkOutput.get("call_name", i)));
			QA_newField.add(QA_map.get(svtkOutput.get("call_name", i)));
			QS_newField.add(QS_map.get(svtkOutput.get("call_name", i)));
			QSE_newField.add(QSE_map.get(svtkOutput.get("call_name", i)));
			QSS_newField.add(QSS_map.get(svtkOutput.get("call_name", i)));
			ploidy_newField.add(ploidy_map.get(svtkOutput.get("call_name", i)));
			strand_newField.add(strand_map.get(svtkOutput.get("call_name", i)));
		}

		ArrayList<String> columnNames = new ArrayList<>();
		columnNames.add("CN");
		columnNames.add("GT");
		columnNames.add("NP");
		columnNames.add("QA");
		columnNames.add("QS");
		columnNames.add("QSE");
		columnNames.add("QSS");
		columnNames.add("ploidy");
		columnNames.add("strand");

		ArrayList<ArrayList<String>> columnValues = new ArrayList<>();
		columnValues.add(CN_newField);
		columnValues.add(GT_newField);
		columnValues.add(NP_newField);
		columnValues.add(QA_newField);
		columnValues.add(QS_newField);
		columnValues.add(QSE_newField);
		columnValues.add(QSS_newField);
		columnValues.add(ploidy_newField);
		columnValues.add(strand_newField);

		svtkOutput.addColumns(columnNames, columnValues);
		svtkOutput.writeFile(match_output, true);
	}

	/**
	 * work in progress, you probably should not use this all this does is call svtk
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void svtk(String svtk_input, String svtk_output) throws IOException, InterruptedException {
		Runtime r = Runtime.getRuntime();
		Process p = r.exec("svtk bedcluster " + svtk_input + " " + svtk_output); // linux
		p.waitFor();
	}

	/**
	 * @deprecated A wrapper to call the full convertVCFsToBEDFormat() using generic
	 *             arguements
	 * @param wd
	 * @param output
	 * @throws IOException
	 */
	public void convertVCFsToBEDFormat(String wd, String output) throws IOException {
		convertVCFsToBEDFormat(wd, output, "genotyped-segments-", ".vcf");
	}

	/**
	 * @deprecated
	 * @param wd          - the working directory where VCFs are stored. VCFs can be
	 *                    in sub-directories.
	 * @param output      - The output path for the final consolidated BED file
	 * @param prefixRegex - prefix to trim from file name, eg "genotyped-segments-"
	 * @param suffixRegex - suffix used to identify VCF files, used also to trim
	 *                    from file name. eg ".vcf"
	 * @throws IOException
	 */
	public void convertVCFsToBEDFormat_old(String wd, String output, String prefixRegex, String suffixRegex)
			throws IOException {
		System.out.println("WARNING: this method is deprecated");
		File[] directories = new File(wd).listFiles(File::isDirectory);
		System.out.println(Arrays.toString(directories));
		System.out.println("n directories: " + directories.length);

		ArrayList<String> all_bed_paths = new ArrayList<>();

		// look through all directories in current working directory
		for (int i = 0; i < directories.length; i++) {
			int cnvNameCounter = 1; // might need to be long
			String currentCluster = directories[i].getAbsolutePath();

			ArrayList<Path> currentClusterVCFsPATH = new ArrayList<>();
			ArrayList<String> currentClusterVCFsNAME = new ArrayList<>();
			ArrayList<String> sampleNames = new ArrayList<>();

			// get a path to all the vcf files, similar to 'ls wd/*vcf'
			Path p = Paths.get(currentCluster);
			final int maxDepth = 1;
			Stream<Path> matches = Files.find(p, maxDepth,
					(path, basicFileAttributes) -> String.valueOf(path).endsWith(suffixRegex));
			matches.filter(s -> s.getFileName().toString().endsWith(suffixRegex)).forEach(currentClusterVCFsPATH::add);
			matches.close();
			for (Path fp : currentClusterVCFsPATH) {
				currentClusterVCFsNAME.add(fp.toAbsolutePath().toString());
				sampleNames.add(fp.getFileName().toString().replaceAll(suffixRegex, ""));
			}

			// if there are no vcfs in this folder, move on to the next folder
			if (currentClusterVCFsPATH.size() < 2) {
				continue;
			}

			// an arraylist containing vcf dataframes
			ArrayList<DataFrame> VCFsArrayList = new ArrayList<>();
			for (String vcfFile : currentClusterVCFsNAME) {
				VCFsArrayList.add(new DataFrame(vcfFile, true, "\\t", "##"));
			}

			// create an array list of to store the new column names, eg "GT, CN, NP, QA"
			// etc
			String[] newColumns = VCFsArrayList.get(0).get("FORMAT", 0).split(":"); // an array list of the new field
																					// names, eg 'GT', 'CN', etc
			ArrayList<String> newColumnNames = new ArrayList<>();
			for (int j = 0; j < newColumns.length; j++) {
				newColumnNames.add(newColumns[j]);
			}

			// to store the full merged bed file
			DataFrame fullClusterBed = new DataFrame();

			for (int k = 0; k < VCFsArrayList.size(); k++) {
				DataFrame currentVCF = new DataFrame();
				ArrayList<String> chr = VCFsArrayList.get(k).getColumn("#CHROM");
				ArrayList<String> start = VCFsArrayList.get(k).getColumn("POS");
				ArrayList<String> end = VCFsArrayList.get(k).getColumn("ID");
				for (int j = 0; j < end.size(); j++) { // Lord help me
					end.set(j, end.get(j).split("_")[end.get(j).split("_").length - 1]);
				}
				ArrayList<String> name = new ArrayList<>();
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					name.add(directories[i].getName() + "_cnv_" + cnvNameCounter);
					cnvNameCounter++;
				}
				ArrayList<String> sample = new ArrayList<>();
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					sample.add(sampleNames.get(k).replace(prefixRegex, "").replace(suffixRegex, ""));
				}

				// create a list of lists to store the new column values
				ArrayList<ArrayList<String>> newColumnValues = new ArrayList<>();
				// this is the raw string from the VCF, the last field name is the raw strings
				// of interest, so length-1 is the right index
				ArrayList<String> columnValues = VCFsArrayList.get(k)
						.getColumn(VCFsArrayList.get(k).fieldNames[VCFsArrayList.get(k).fieldNames.length - 1]);

				// sanity check
				if (columnValues.size() != VCFsArrayList.get(k).nrow())
					System.out.println("ERROR 2482469776");

				// add a list to store every new field we will be parsing/adding
				for (int j = 0; j < newColumnNames.size(); j++) {
					newColumnValues.add(new ArrayList<>());
				}

				// for every row in the current VCF
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					// for every new column
					for (int l = 0; l < newColumnNames.size(); l++) {
						// to the target column add the value in the split[] index
						newColumnValues.get(l).add(columnValues.get(j).split(":")[l]);
					}
				}

				ArrayList<String> columnNamesToAdd = new ArrayList<>();
				ArrayList<ArrayList<String>> columnValuesToAdd = new ArrayList<>();

				columnNamesToAdd.add("chr");
				columnValuesToAdd.add(chr);
				columnNamesToAdd.add("start");
				columnValuesToAdd.add(start);
				columnNamesToAdd.add("end");
				columnValuesToAdd.add(end);
				columnNamesToAdd.add("name");
				columnValuesToAdd.add(name);
				columnNamesToAdd.add("sample");
				columnValuesToAdd.add(sample);
				columnNamesToAdd.addAll(newColumnNames);
				columnValuesToAdd.addAll(newColumnValues);
				boolean temp_1 = currentVCF.addColumns(columnNamesToAdd, columnValuesToAdd);
				if (!temp_1) {
					System.out.println("eror 42721191233: " + currentClusterVCFsNAME.get(k));
				}
				// System.out.println(temp_1);
				// label dups and dels
				ArrayList<String> svtype = new ArrayList<>();
				// VCFsArrayList.get(k).nrow() should now be equal to currentVCF?
				// System.out.println(currentVCF.columnMapping.toString());
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					// System.out.println(currentVCF.get("chr").size());
					String currChr = currentVCF.getColumn("chr").get(j);
					String svt = "NA_j";
					if (currChr.contains("X") || currChr.contains("x")) {
						svt = Integer.parseInt(currentVCF.getColumn("CN").get(j)) >= 2 ? "DUP" : "DEL";
					} else if (currChr.contains("Y") || currChr.contains("y")) {
						svt = Integer.parseInt(currentVCF.getColumn("CN").get(j)) >= 1 ? "DUP" : "DEL";
					} else {
						svt = Integer.parseInt(currentVCF.getColumn("CN").get(j)) >= 2 ? "DUP" : "DEL";
					}
					svtype.add(svt);
				}
				currentVCF.addColumn("svtype", svtype);

				// remove rows where GT=0, this is not the most efficient way to remove indexes
				// from arraylist
				ArrayList<Integer> gtIsZero = new ArrayList<>();
				for (int j = currentVCF.nrow() - 1; j >= 0; j--) {
					if (currentVCF.get("GT", j).equals("0")) {
						gtIsZero.add(j);
					}
				}
				Collections.sort(gtIsZero, Collections.reverseOrder());
				for (int j : gtIsZero) {
					currentVCF.df.remove(j);
				}

				if (fullClusterBed.df.size() == 0) {
					fullClusterBed = currentVCF;
				} else {
					boolean ctrl = fullClusterBed.rbind(currentVCF);
					if (ctrl == false) {
						System.out.println("error 23089438  " + k + " " + i);
						System.out.println();
					}
				}
			}

			System.out.println(currentCluster + "\\" + directories[i].getName() + ".bed");
			fullClusterBed.writeBed(currentCluster + "\\" + directories[i].getName() + ".java.bed");

			all_bed_paths.add(currentCluster + "\\" + directories[i].getName() + ".java.bed");
		}

		DataFrame fullMergedBed = new DataFrame(all_bed_paths.get(0), true, "\\t", "@");
		for (int i = 1; i < all_bed_paths.size(); i++) {
			fullMergedBed.rbind(new DataFrame(all_bed_paths.get(i), true, "\\t", "@"));
		}
		System.out.println(wd + output);
		fullMergedBed.writeFile(wd + output, true);

	}

	/**
	 * 
	 * @param wd          - the working directory where VCFs are stored. VCFs can be
	 *                    in sub-directories.
	 * @param output      - The output path for the final consolidated BED file
	 * @param prefixRegex - prefix to trim from file name, eg "genotyped-segments-"
	 * @param suffixRegex - suffix used to identify VCF files, used also to trim
	 *                    from file name. eg ".vcf"
	 * @throws IOException
	 */
	public void convertVCFsToBEDFormat(String wd, String output, String prefixRegex, String suffixRegex)
			throws IOException {
		// withploidy
		File[] directories = new File(wd).listFiles(File::isDirectory);
		System.out.println(Arrays.toString(directories));
		System.out.println("n directories: " + directories.length);

		ArrayList<String> all_bed_paths = new ArrayList<>();

		// look through all directories in current working directory
		for (int i = 0; i < directories.length; i++) {
			int cnvNameCounter = 1; // might need to be long
			String currentCluster = directories[i].getAbsolutePath();

			ArrayList<Path> currentClusterVCFsPATH = new ArrayList<>();
			ArrayList<String> currentClusterVCFsNAME = new ArrayList<>();
			ArrayList<String> sampleNames = new ArrayList<>();

			// get a path to all the vcf files, similar to 'ls wd/*vcf'
			Path p = Paths.get(currentCluster);
			final int maxDepth = 1;
			Stream<Path> matches = Files.find(p, maxDepth,
					(path, basicFileAttributes) -> String.valueOf(path).endsWith(suffixRegex));
			matches.filter(s -> s.getFileName().toString().endsWith(suffixRegex)).forEach(currentClusterVCFsPATH::add);
			matches.close();
			for (Path fp : currentClusterVCFsPATH) {
				currentClusterVCFsNAME.add(fp.toAbsolutePath().toString());
				sampleNames.add(fp.getFileName().toString().replaceAll(suffixRegex, ""));
			}

			// if there are no vcfs in this folder, move on to the next folder
			if (currentClusterVCFsPATH.size() < 2) {
				continue;
			}

			// an arraylist containing vcf dataframes
			ArrayList<DataFrame> VCFsArrayList = new ArrayList<>();
			for (String vcfFile : currentClusterVCFsNAME) {
				VCFsArrayList.add(new DataFrame(vcfFile, true, "\\t", "##"));
			}

			// create an array list of to store the new column names, eg "GT, CN, NP, QA"
			// etc
			String[] newColumns = VCFsArrayList.get(0).get("FORMAT", 0).split(":"); // an array list of the new field
																					// names, eg 'GT', 'CN', etc
			ArrayList<String> newColumnNames = new ArrayList<>();
			for (int j = 0; j < newColumns.length; j++) {
				newColumnNames.add(newColumns[j]);
			}

			// to store the full merged bed file
			DataFrame fullClusterBed = new DataFrame();

			for (int k = 0; k < VCFsArrayList.size(); k++) {
				DataFrame currentVCF = new DataFrame();
				ArrayList<String> chr = VCFsArrayList.get(k).getColumn("#CHROM");
				ArrayList<String> start = VCFsArrayList.get(k).getColumn("POS");
				ArrayList<String> end = VCFsArrayList.get(k).getColumn("ID");
				for (int j = 0; j < end.size(); j++) { // Lord help me
					end.set(j, end.get(j).split("_")[end.get(j).split("_").length - 1]);
				}
				ArrayList<String> name = new ArrayList<>();
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					name.add(directories[i].getName() + "_cnv_" + cnvNameCounter);
					cnvNameCounter++;
				}
				ArrayList<String> sample = new ArrayList<>();
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					sample.add(sampleNames.get(k).replace(prefixRegex, "").replace(suffixRegex, ""));
				}

				// create a list of lists to store the new column values
				ArrayList<ArrayList<String>> newColumnValues = new ArrayList<>();
				// this is the raw string from the VCF, the last field name is the raw strings
				// of interest, so length-1 is the right index
				ArrayList<String> columnValues = VCFsArrayList.get(k)
						.getColumn(VCFsArrayList.get(k).fieldNames[VCFsArrayList.get(k).fieldNames.length - 1]);

				// sanity check
				if (columnValues.size() != VCFsArrayList.get(k).nrow())
					System.out.println("ERROR 2482469776");

				// add a list to store every new field we will be parsing/adding
				for (int j = 0; j < newColumnNames.size(); j++) {
					newColumnValues.add(new ArrayList<>());
				}

				// for every row in the current VCF
				for (int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
					// for every new column
					for (int l = 0; l < newColumnNames.size(); l++) {
						// to the target column add the value in the split[] index
						newColumnValues.get(l).add(columnValues.get(j).split(":")[l]);
					}
				}

				ArrayList<String> columnNamesToAdd = new ArrayList<>();
				ArrayList<ArrayList<String>> columnValuesToAdd = new ArrayList<>();

				columnNamesToAdd.add("chr");
				columnValuesToAdd.add(chr);
				columnNamesToAdd.add("start");
				columnValuesToAdd.add(start);
				columnNamesToAdd.add("end");
				columnValuesToAdd.add(end);
				columnNamesToAdd.add("name");
				columnValuesToAdd.add(name);
				columnNamesToAdd.add("sample");
				columnValuesToAdd.add(sample);
				columnNamesToAdd.addAll(newColumnNames);
				columnValuesToAdd.addAll(newColumnValues);
				boolean temp_1 = currentVCF.addColumns(columnNamesToAdd, columnValuesToAdd);
				if (!temp_1) {
					System.out.println("eror 42721191233: " + currentClusterVCFsNAME.get(k));
				}
				// System.out.println(temp_1);
				// label dups and dels
				// ArrayList<String> svtype = new ArrayList<>();
				// //VCFsArrayList.get(k).nrow() should now be equal to currentVCF?
				// //System.out.println(currentVCF.columnMapping.toString());
				// for(int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
				// //System.out.println(currentVCF.get("chr").size());
				// String currChr = currentVCF.get("chr").get(j);
				// String svt = "NA_j";
				// if(currChr.contains("X") || currChr.contains("x")) {
				// svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 2 ? "DUP" : "DEL";
				// } else if(currChr.contains("Y") || currChr.contains("y")) {
				// svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 1 ? "DUP" : "DEL";
				// } else {
				// svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 2 ? "DUP" : "DEL";
				// }
				// svtype.add(svt);
				// }
				// currentVCF.addColumn("svtype", svtype);

				// calculate x and y ploidy
				int sumWidthX = 0;
				int sumProdX = 0;
				int sumWidthY = 0;
				int sumProdY = 0;
				for (int j = 0; j < currentVCF.nrow(); j++) {
					if (currentVCF.get("chr", j).contains("x") || currentVCF.get("chr", j).contains("X")) {
						int currWidth = (Integer.parseInt(currentVCF.get("end", j))
								- Integer.parseInt(currentVCF.get("start", j)));
						int currCN = Integer.parseInt(currentVCF.get("CN", j));
						sumWidthX += currWidth;
						sumProdX += currWidth * currCN;
					} else if (currentVCF.get("chr", j).contains("y") || currentVCF.get("chr", j).contains("Y")) {
						int currWidth = (Integer.parseInt(currentVCF.get("end", j))
								- Integer.parseInt(currentVCF.get("start", j)));
						int currCN = Integer.parseInt(currentVCF.get("CN", j));
						sumWidthY += currWidth;
						sumProdY += currWidth * currCN;
					}
				}
				int ploidy_x = (int) ((double) (sumProdX) / sumWidthX + 0.5);
				int ploidy_y = (int) ((double) (sumProdY) / sumWidthY + 0.5);
				if (ploidy_x + ploidy_y != 2) {
					if (ploidy_y > 0)
						ploidy_x = 1;
					if (ploidy_y == 0)
						ploidy_x = 2;
					if (ploidy_y > 1)
						ploidy_y = 1;
				}
				String ploidyX = "" + ploidy_x;
				String ploidyY = "" + ploidy_y;

				ArrayList<String> ploidy = new ArrayList<>();
				ArrayList<String> strand = new ArrayList<>();
				ArrayList<String> call = new ArrayList<>();
				for (int j = 0; j < currentVCF.nrow(); j++) {
					if (currentVCF.get("chr", j).contains("x") || currentVCF.get("chr", j).contains("X")) {
						ploidy.add(ploidyX);
					} else if (currentVCF.get("chr", j).contains("y") || currentVCF.get("chr", j).contains("Y")) {
						ploidy.add(ploidyY);
					} else {
						ploidy.add("2");
					}
				}
				currentVCF.addColumn("ploidy", ploidy);

				// remove rows where GT=0 and where CN==ploidy, this is not the most efficient
				// way to remove indexes from arraylist
				ArrayList<Integer> gtIsZeroOrCNisPloidy = new ArrayList<>();
				for (int j = currentVCF.nrow() - 1; j >= 0; j--) {
					if (currentVCF.get("GT", j).equals("0")
							|| currentVCF.get("CN", j).equals(currentVCF.get("ploidy", j))) {
						gtIsZeroOrCNisPloidy.add(j);
					}
				}
				Collections.sort(gtIsZeroOrCNisPloidy, Collections.reverseOrder());
				for (int j : gtIsZeroOrCNisPloidy) {
					currentVCF.df.remove(j);
				}

				// currentVCF.get("ploidy", j)
				for (int j = 0; j < currentVCF.nrow(); j++) {
					if (Integer.parseInt(currentVCF.get("CN", j)) > Integer.parseInt(currentVCF.get("ploidy", j))) {
						call.add("DUP");
						strand.add("+");
					} else {
						call.add("DEL");
						strand.add("-");
					}
				}
				currentVCF.addColumn("svtype", call);
				currentVCF.addColumn("strand", strand);

				if (fullClusterBed.df.size() == 0) {
					fullClusterBed = currentVCF;
				} else {
					boolean ctrl = fullClusterBed.rbind(currentVCF);
					if (ctrl == false) {
						System.out.println("error 23089438  " + k + " " + i);
						System.out.println();
					}
				}
			}

			System.out.println(currentCluster + "\\" + directories[i].getName() + ".bed");
			fullClusterBed.writeBed(currentCluster + "\\" + directories[i].getName() + ".ploidy.java.bed");

			all_bed_paths.add(currentCluster + "\\" + directories[i].getName() + ".ploidy.java.bed");
		}

		DataFrame fullMergedBed = new DataFrame(all_bed_paths.get(0), true, "\\t", "@");
		for (int i = 1; i < all_bed_paths.size(); i++) {
			fullMergedBed.rbind(new DataFrame(all_bed_paths.get(i), true, "\\t", "@"));
		}
		System.out.println(wd + output);
		fullMergedBed.writeFile(wd + output, true);

	}

	/**
	 * 
	 * @param entityPath
	 * @param wd
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void downloadSegmentsVCFs(String entityPath, String wd) throws IOException, InterruptedException {
		downloadSegmentsVCFs(entityPath, wd, "segments_vcfs");
	}

	/**
	 * @param entityPath               - full path to the entity file, eg
	 *                                 "sample_set_entity.tsv"
	 * @param wd                       - working directory, where the files should
	 *                                 be downloaded to
	 * @param segments_vcfs_columnName - the name of the column in the entity file
	 *                                 containing the segment_vcfs, eg
	 *                                 "segments_vcfs"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void downloadSegmentsVCFs(String entityPath, String wd, String segments_vcfs_columnName)
			throws IOException, InterruptedException {
		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
		for (int i = 0; i < sampleSetEntity.nrow(); i++) {
			String individualToDownload = sampleSetEntity.get("entity:sample_set_id", i);
			String pathInGoogleBucket = sampleSetEntity.get(segments_vcfs_columnName, i).replaceAll("\\[|\\]|\"", "")
					.replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");

			// if(pathInGoogleBucket.trim().equals("") || (files.length==1 &&
			// files[0].trim().equals("")) || files.length==0) {continue;}
			if (files.length == 1 && (files[0].trim().equals("") || files[0].trim().equals("0"))) {
				continue;
			}

			String temp_suffix = "_files_to_download.txt";

			File file = new File(wd + individualToDownload + temp_suffix);
			BufferedWriter output = new BufferedWriter(new FileWriter(file));
			for (String filePath : files) {
				output.write(filePath + "\n");
			}
			output.close();
			String pathToDownloadTo = wd + individualToDownload + "/";
			File dirToDownloadTo = new File(pathToDownloadTo);
			if (!dirToDownloadTo.exists()) {
				dirToDownloadTo.mkdir();
			}

			String wdUnix = wd.replace("C:", "/mnt/c");
			String pathToDownloadToUnix = wdUnix + individualToDownload;

			downloadFiles(wdUnix + individualToDownload + temp_suffix, pathToDownloadToUnix);

			// //this is only tested on windows, and it works so far. ITS VERY DELICATE.
			String[] command;

			// unzip vcfs, working as of 12/3/19
			if (System.getProperty("os.name").contains("Windows")) {
				String[] dosCommand = { "bash", "-c", "'gunzip", "-k", pathToDownloadToUnix + "/*'" };
				command = dosCommand;
			} else {
				String[] unixCommand = { "gunzip", "-k" + pathToDownloadToUnix + "/*'" };
				command = unixCommand;
			}
			System.out.println(String.join(" ", command));
			try {
				new ProcessBuilder(command).inheritIO().start().waitFor();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * this is just a wrapper for the buffered version, which sets the buffer size
	 * to 20
	 * 
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @param countsRegex
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrixBuffered(String sourceFolder, String OUTPUT_PATH, String countsRegex)
			throws IOException, InterruptedException {
		getCountsMatrixBuffered(sourceFolder, OUTPUT_PATH, countsRegex, 20);
	}

	/**
	 * write buffer after every n files is read per thread, uses less memory as less
	 * count data is held in memory at any given time however, uses more I/O, which
	 * will be slower overall
	 * 
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @param countsRegex
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrixBuffered(String sourceFolder, String OUTPUT_PATH, String countsRegex, int BUFFER_SIZE)
			throws IOException, InterruptedException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();

		System.out.print("finding files.\t");
		Path p = Paths.get(sourceFolder);
		final int maxDepth = 10;
		Stream<Path> matches = Files.find(p, maxDepth,
				(path, basicFileAttributes) -> String.valueOf(path).endsWith(countsRegex));
		matches.filter(s -> s.getFileName().toString().contains(countsRegex)).forEach(barcodeCountsPaths::add);
		matches.close();
		for (Path fp : barcodeCountsPaths) {
			barcodeCountsFiles.add(fp.toAbsolutePath().toString());
			sampleNames.add(fp.getFileName().toString().replaceAll(countsRegex, ""));
		}
		System.out.println("found " + sampleNames.size() + " files");

		System.out.println("reading files. \t");

		// toRead functions as a synchronized queue so each thread knows which file to
		// read next
		List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
		for (int i = 0; i < sampleNames.size(); i++) {
			toRead.add(i);
		}
		int N_THREADS = Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		int totalNFiles = toRead.size();
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {
					File file = new File(OUTPUT_PATH + "_" + Thread.currentThread().getId());
					BufferedWriter output = null;
					try {
						output = new BufferedWriter(new FileWriter(file));
					} catch (IOException e1) {
						e1.printStackTrace();
					}

					HashMap<String, String> doneReading = new HashMap<>(); // map of sample name to line string

					int bufferCounter = 0;

					ArrayList<String> labels = new ArrayList<>();
					DataFrame labelsDF = null;
					try {
						labelsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
					} catch (IOException e1) {
						e1.printStackTrace();
					}
					for (int i = 0; i < labelsDF.nrow(); i++) {
						labels.add(labelsDF.get("CONTIG", i) + "_" + labelsDF.get("START", i) + "_"
								+ labelsDF.get("END", i));
					}
					labelsDF = null;

					for (int i = 0; i < labels.size(); i++) {
						try {
							output.write(labels.get(i));
						} catch (IOException e) {
							e.printStackTrace();
						}
						if (i != labels.size() - 1) {
							try {
								output.write("\t");
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
					}
					try {
						output.write("\n");
					} catch (IOException e1) {
						e1.printStackTrace();
					}

					while (toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading.put(sampleNames.get(currentFile),
									String.join("\t", countsDF.getColumn("COUNT")));
							progressPercentage(totalNFiles - toRead.size(), totalNFiles,
									barcodeCountsFiles.get(currentFile));
							countsDF = null; // java gc is a fickle mistress

							if (bufferCounter >= BUFFER_SIZE) {
								StringBuilder line = new StringBuilder();
								for (String snKey : doneReading.keySet()) {
									line.append(snKey); // sample name key
									line.append("\t");
									line.append(doneReading.get(snKey));
									line.append("\n");
								}
								output.write(line.toString());
								doneReading = new HashMap<>();
								bufferCounter = 0;
							} else {
								bufferCounter++;
							}

						} catch (IOException e) {
							e.printStackTrace();
						}
					}

					// flush remaining buffer
					StringBuilder line = new StringBuilder();
					for (String snKey : doneReading.keySet()) {
						line.append(snKey); // sample name key
						line.append("\t");
						line.append(doneReading.get(snKey));
						line.append("\n");
					}
					try {
						output.write(line.toString());
					} catch (IOException e2) {
						e2.printStackTrace();
					}
					doneReading = new HashMap<>();
					bufferCounter = 0;

					try {
						output.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

	}

	/**
	 * this is just a wrapper for getCountsMatrix() which sets the default regex to
	 * ".barcode.counts.tsv"
	 * 
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrix(String sourceFolder, String OUTPUT_PATH) throws IOException, InterruptedException {
		getCountsMatrix(sourceFolder, OUTPUT_PATH, ".barcode.counts.tsv");
	}

	/**
	 * This is an improvement upon the previous version as it now more accurately
	 * forces GC to remove the DF after its been read and only keep the important
	 * field Reads in barcode counts files in parallel, and writes a counts matrix
	 * 
	 * @param sourceFolder - directory where files are located in, files can be in
	 *                     sub-directories.
	 * @param OUTPUT_PATH  - the full output path where the matrix file will be
	 *                     written to
	 * @param countsRegex  - the regex suffix to identify counts files
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrix(String sourceFolder, String OUTPUT_PATH, String countsRegex)
			throws IOException, InterruptedException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();

		System.out.print("finding files.\t");
		Path p = Paths.get(sourceFolder);
		final int maxDepth = 10;
		Stream<Path> matches = Files.find(p, maxDepth,
				(path, basicFileAttributes) -> String.valueOf(path).endsWith(countsRegex));
		matches.filter(s -> s.getFileName().toString().contains(countsRegex)).forEach(barcodeCountsPaths::add);
		matches.close();
		for (Path fp : barcodeCountsPaths) {
			barcodeCountsFiles.add(fp.toAbsolutePath().toString());
			sampleNames.add(fp.getFileName().toString().replaceAll(countsRegex, ""));
		}
		System.out.println("found " + sampleNames.size() + " files");

		System.out.println("reading files. \t");

		// toRead functions as a synchronized queue so each thread knows which file to
		// read next
		// each dataframe is then mapped to a unique index number so that the arraylist
		// of dataframes
		// can be assembled again in order even though each one is read out of order
		List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
		Map<Integer, String> doneReading2 = new ConcurrentHashMap<>();
		for (int i = 0; i < sampleNames.size(); i++) {
			toRead.add(i);
		}
		int N_THREADS = Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		int totalNFiles = toRead.size();
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {
					while (toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading2.put(currentFile, String.join("\t", countsDF.getColumn("COUNT")));
							progressPercentage(totalNFiles - toRead.size(), totalNFiles,
									barcodeCountsFiles.get(currentFile));
							countsDF = null; // java gc is a fickle mistress
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

		ArrayList<String> countsArrayList = new ArrayList<>();
		for (int i = 0; i < doneReading2.size(); i++) {
			countsArrayList.add(doneReading2.get(i));
		}

		System.out.println("done reading files");

		// 'labels' is the exon labels, / column names
		System.out.print("generating counts matrix. \t");
		ArrayList<String> labels = new ArrayList<>();
		DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
		for (int i = 0; i < countsDF.nrow(); i++) {
			labels.add(countsDF.get("CONTIG", i) + "_" + countsDF.get("START", i) + "_" + countsDF.get("END", i));
		}
		System.out.println("done generating counts matrix");

		// for sanity checking, number should all appropriately match
		System.out.println("sampleNames.size()\t" + sampleNames.size());
		System.out.println("countsArrayList.size()\t" + countsArrayList.size());
		System.out.println("labels.size()\t" + labels.size());
		System.out.println("writing counts matrix");

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		for (int i = 0; i < labels.size(); i++) {
			output.write(labels.get(i));
			if (i != labels.size() - 1) {
				output.write("\t");
			}
		}
		output.write("\n");

		if (countsArrayList.size() != sampleNames.size()) {
			System.out.println("error 51");
		}

		// written in such a way as to be readable in R, adding rownames for the sample
		// ID
		for (int i = 0; i < sampleNames.size(); i++) {
			StringBuilder line = new StringBuilder();
			line.append(sampleNames.get(i) + "\t");
			line.append(countsArrayList.get(i));
			if (i != sampleNames.size() - 1) {
				line.append("\n");
			}
			output.write(line.toString());
		}
		output.close();
	}

	/**
	 * a wrapper for the full getBarcodeCounts() with generic args
	 * 
	 * @param entityPath
	 * @param wd
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getBarcodeCounts(String entityPath, String wd) throws IOException, InterruptedException {
		getBarcodeCounts(entityPath, wd, "output_counts_barcode");
	}

	/**
	 * Downloads the barcode counts file from terra/firecloud
	 * 
	 * @param entityPath                  - full path to the entity file, eg
	 *                                    "/documents/sample_set_entity.tsv"
	 * @param wd                          - the directory to download the counts
	 *                                    files to
	 * @param output_counts_barcode_regex - the name of the column containing the
	 *                                    counts files, eg "output_counts_barcode"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getBarcodeCounts(String entityPath, String wd, String output_counts_barcode_regex)
			throws IOException, InterruptedException {
		DataFrame entityDF = new DataFrame(entityPath, true, "\\t", "@");
		System.out.println(entityDF.nrow());
		System.out.println(entityDF.columnMapping.toString());
		for (int i = 0; i < entityDF.nrow(); i++) {
			String individualToDownload = entityDF.get("entity:sample_set_id", i);
			String pathInGoogleBucket = entityDF.get(output_counts_barcode_regex, i).replaceAll("\\[|\\]|\"", "")
					.replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");

			if (files.length == 1 && (files[0].trim().equals("") || files[0].trim().equals("0"))) {
				continue;
			}

			String temp_suffix = "_temp_files_to_download.txt";

			File file = new File(wd + individualToDownload + temp_suffix);
			BufferedWriter output = new BufferedWriter(new FileWriter(file));
			for (String filePath : files) {
				output.write(filePath + "\n");
			}
			output.close();

			String pathToDownloadTo = wd + individualToDownload + "/";
			File dirToDownloadTo = new File(pathToDownloadTo);
			if (!dirToDownloadTo.exists()) {
				dirToDownloadTo.mkdir();
			}

			String wdUnix = wd.replace("C:", "/mnt/c");
			String pathToDownloadToUnix = wdUnix + individualToDownload;

			downloadFiles(wdUnix + individualToDownload + temp_suffix, pathToDownloadToUnix);

		}
	}

	/**
	 * takes in a file containing a list of files to download. downloads those
	 * files. this has only been tested on windows, it is VERY DELICATE
	 * 
	 * @param filesFile
	 * @throws InterruptedException
	 */
	public static void downloadFiles(String filesFile, String pathToDownloadToUnix) throws InterruptedException {
		String[] command;
		if (System.getProperty("os.name").contains("Windows")) {
			String[] dosCommand = { "bash", "-c", "'cat", filesFile, "|", "gsutil", "-m", "cp", "-I",
					pathToDownloadToUnix, "'" };
			command = dosCommand;
		} else {
			System.out.println(
					"WARNING: this download functionality has not been tested on macOS, please report any issues");
			String[] unixCommand = { "cat", filesFile, "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix };
			command = unixCommand;
		}
		System.out.println(String.join(" ", command));

		try {
			new ProcessBuilder(command).inheritIO().start().waitFor();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * https://stackoverflow.com/questions/852665/command-line-progress-bar-in-java
	 * 
	 * @param remain
	 * @param total
	 */
	public static void progressPercentage(int remain, int total, String message) {
		if (remain > total) {
			throw new IllegalArgumentException();
		}
		int maxBareSize = 10; // 10unit for 100%
		int remainProcent = ((100 * remain) / total) / maxBareSize;
		char defaultChar = '-';
		String icon = "*";
		String bare = new String(new char[maxBareSize]).replace('\0', defaultChar) + "]";
		StringBuilder bareDone = new StringBuilder();
		bareDone.append("[");
		for (int i = 0; i < remainProcent; i++) {
			bareDone.append(icon);
		}
		String bareRemain = bare.substring(remainProcent, bare.length());
		System.out.print("\r" + bareDone + bareRemain + " " + remainProcent * 10 + "%  " + remain + " / " + total + "  "
				+ message);
		if (remain == total) {
			System.out.print("\n");
		}
	}

	/**
	 * @return returns the current JVM heap size in human readable format
	 */
	public static String getHeapSize() {
		return humanReadableByteCount(Runtime.getRuntime().totalMemory());
	}

	/**
	 * @param bytes
	 * @return human readable version of bytes
	 */
	public static String humanReadableByteCount(long bytes) {
		if (bytes < 1024)
			return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp - 1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}

	// returns the mean of an arraylist, using integers, so will trucnate
	public static int mean(ArrayList<Integer> al) {
		int sum = 0;
		for (Integer i : al) {
			sum += i;
		}
		return (sum / al.size());
	}

	/**
	 * returns the median of an integer arraylist
	 * 
	 * @param al
	 * @return
	 */
	public static int median(ArrayList<Integer> al) {
		Collections.sort(al);
		return (al.get(al.size() / 2));
	}

	/**
	 * returns to string arraylist with values separated by sep
	 * 
	 * @param arrayList
	 * @param sep
	 * @return
	 */
	public String arrayToStringWithSep(ArrayList<String> arrayList, String sep) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < arrayList.size() - 1; i++) {
			sb.append(arrayList.get(i));
			sb.append(sep);
		}
		sb.append(arrayList.get(arrayList.size() - 1));
		return sb.toString();
	}

	public static void print(String str) {
		System.out.println(str);
	}

	public static void print(String[] str) {
		for (String s : str) {
			System.out.print(s + "\t");
		}
		System.out.println();
	}

	public static String getTimeHumanReadable(long time) {
		long seconds = time % 60;
		long minutes = (time / 60) % 60;
		long hours = (time / 3600);

		return "Hours: " + hours + "  Min: " + minutes + "  Sec: " + seconds;
	}

	public static void start() {
		startTime = System.nanoTime();
	}

	public static void stop(String str) {
		System.out.println(str);
		stop();
	}

	public static void stop() {
		endTime = System.nanoTime();
		long secondsElapsed = (long) ((endTime - startTime) / 1000000000);
		System.out.println("Time elapsed: " + secondsElapsed + "s");
		System.out.println(getTimeHumanReadable(secondsElapsed));
	}

	public void checkVersion() throws IOException {
		// https://stackoverflow.com/questions/9489726/get-raw-text-from-html
		URL u = new URL("https://raw.githubusercontent.com/theisaacwong/talkowski/master/gCNV/Readme.md");
		URLConnection conn = u.openConnection();
		conn.setReadTimeout(15000);
		BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
		StringBuffer buffer = new StringBuffer();
		String inputLine;
		while ((inputLine = in.readLine()) != null)
			buffer.append(inputLine);
		in.close();
		String line = buffer.toString();
		Pattern versionPattern = Pattern.compile("(?<=Version: )\\d+.\\d+");
		Matcher versionMatcher = versionPattern.matcher(line);
		String versionString = versionMatcher.find() ? versionMatcher.group() : "-1";
		if (versionString.trim().equals(VERSION)) {
			System.out.println("(Most current version)");
		} else {
			System.out.println("Waring. You are currently running version " + VERSION + ". The most recent version is "
					+ versionString + ". Updating is strongly recommended.");
		}
	}

	/**
	 * this started out as only three options but has quickly spiralled out of
	 * control, I'm so sorry
	 */
	public void printOptions() {
		System.out.println("java -jar gCNV_helper.jar [Command] [required argument(s)] ");
		System.out.println();
		InputStream in = this.getClass().getResourceAsStream("/helpMenu"); 
		Scanner sc = new Scanner(in, "UTF-8");
		String line = "";
		while(sc.hasNextLine()) {
			line = sc.nextLine();
			System.out.println(line);
		}
		sc.close();
	}

	public void runBetter(String[] args) throws IOException, InterruptedException {
		System.out.println();
		if (args.length == 0 || args[0].contains("-help") || args[0].contains("-h")) {
			printOptions();
			return;
		}
		String toolName = args[0];
		switch (toolName) {
		case "annotateWithGenes" -> {
			String mode = args[1];
			var genes = parseGTFFile(args[2]);
			String GCNV_INPUT = args[3];
			String OUTPUT_PATH = args[4];
			String typeColumnName = args[5];
			if(mode.equals("strict")) {
				this.annotateGenesPercentBased(genes, GCNV_INPUT, OUTPUT_PATH, typeColumnName);
			} else if(mode.equals("any")) {
				this.annotateGenesAnyOverlap(genes, GCNV_INPUT, OUTPUT_PATH, typeColumnName);
			} else {
				System.out.println("wrong input");
			}
		}
		case "condenseBedtoolsIntersect" -> {
			String INPUT_PATH = args[1];
			String OUTPUT_PATH = args[2];
			String columnsToHashOnString = args[3];
			String columnsToMergeString = args[4];
			String columnsToKeepString = args[5];
			this.condenseBedtoolsIntersect(INPUT_PATH, OUTPUT_PATH, columnsToHashOnString, columnsToMergeString, columnsToKeepString);
		}
		case "convertToEnsemble" -> {
			String GCNV_INPUT = args[1];
			String OUTPUT = args[2];
			String geneColumnName = args[3];
			String gencodeGTF = args[4];
			this.convertToEnsemble(GCNV_INPUT, OUTPUT, geneColumnName, gencodeGTF);
		}
		case "convertVCFsToBEDFormat" -> {
			String wd = args[1];
			String output = args[2];
			String prefixRegex = args[3];
			String suffixRegex = args[4];
			this.convertVCFsToBEDFormat(wd, output, prefixRegex, suffixRegex);
		}
		case "countExons" -> {
			String GCNV_INPUT = args[1];
			String genesColumnName = args[2];
			String gencodeGTF = args[3];
			String output = args[4];
			this.countExons(GCNV_INPUT, genesColumnName, gencodeGTF, output);
		}
		case "defragment" -> {
			String match_output = args[1];
			String defragmentation_output = args[2];
			this.defragment(match_output, defragmentation_output);
		}
		case "downloadSegmentsVCFs" -> {
			String entityPath = args[1];
			String wd = args [2];
			String segments_vcfs_columnName = args[3];
			this.downloadSegmentsVCFs(entityPath, wd, segments_vcfs_columnName);
		}
		case "getBarcodeCounts" -> {
			String entityPath = args[1];
			String wd = args[2];
			String output_counts_barcode_regex = args[3];
			this.getBarcodeCounts(entityPath, wd, output_counts_barcode_regex);
		}
		case "getCountsMatrix" -> {
			String sourceFolder = args[1];
			String OUTPUT_PATH = args[2];
			String countsRegex = args[3];
			this.getCountsMatrix(sourceFolder, OUTPUT_PATH, countsRegex);
		}
		case "getCountsMatrixBuffered" -> {
			String sourceFolder = args[1];
			String OUTPUT_PATH = args[2];
			String countsRegex = args[3];
			int BUFFER_SIZE = Integer.parseInt(args[4]);
			this.getCountsMatrixBuffered(sourceFolder, OUTPUT_PATH, countsRegex, BUFFER_SIZE);
		}
		case "getPerSampleMetrics" -> {
			String INPUT_PATH = args[1];
			String fieldColumn = args[2];
			String OUTPUT_PATH = args[3];
			boolean writeFile = true;
			this.getPerSampleMetrics(INPUT_PATH, fieldColumn, OUTPUT_PATH, writeFile);
		}
		case "jointCallVCF" -> {
			String INPUT_PATH = args[1];
			String OUTPUT_PATH = args[2];
			String classColumn = args[3];
			String observationColumn = args[4];
			String sep = args[5];
			String svtypeColumn = args[6];
			this.jointCallVCF(INPUT_PATH, OUTPUT_PATH, classColumn, observationColumn, sep, svtypeColumn);
		}
		case "labelMedianQS" -> {
			String INPUT_PATH = args[1];
			String variantColumn = args[2];
			String qsColumn = args[3];
			String OUTPUT_PATH = args[4];
			this.labelMedianQS(INPUT_PATH, variantColumn, qsColumn, OUTPUT_PATH);
		}
		case "simplify" -> {
			String input = args[1];
			String OUTPUT_PATH = args[2];
			String columnToAggregateBy = args[3];
			String columnsToAggregate = args[4];
			String columnToSplitBy = args[5];
			this.simplify(input, OUTPUT_PATH, columnToAggregateBy, columnsToAggregate, columnToSplitBy);
		}
		case "subsetAnnotations" -> {
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
		case "validateSubsetAnnotations" -> {
			String gtfFile = args[1];
			ArrayList<String> annotationSubsets = new ArrayList<>();
			for (int i = 2; i < args.length; i++) {
				annotationSubsets.add(args[i]);
			}
			this.validateSubsetAnnotations(gtfFile, annotationSubsets);
		}
        default -> {
            System.out.println("unknown command: " + toolName);
        }
		
		}
		
		
	}
	
	// /public void annotateWithGenes(String GCNV_INPUT, String ANNO_INPUT, String
	// OUTPUT_PATH)
	public void run(String[] args) throws IOException, InterruptedException {
		System.out.println();
		if (args.length == 0 || args[0].contains("-help") || args[0].contains("-h")) {
			printOptions();
		} else if (args[0].equals("getBarcodeCounts")) {
			if (args.length == 3) {
				getBarcodeCounts(args[1], args[2]);
			} else if (args.length == 4) {
				getBarcodeCounts(args[1], args[2], args[3]);
			}
		} else if (args[0].equals("getCountsMatrix")) {
			if (args.length == 3) {
				getCountsMatrix(args[1], args[2]);
			} else if (args.length == 4) {
				getCountsMatrix(args[1], args[2], args[3]);
			}
		} else if (args[0].equals("downloadSegmentsVCFs")) {
			if (args.length == 3) {
				downloadSegmentsVCFs(args[1], args[2]);
			} else if (args.length == 4) {
				downloadSegmentsVCFs(args[1], args[2], args[3]);
			}
		} else if (args[0].equals("convertVCFsToBEDFormat")) {
			if (args.length == 3) {
				convertVCFsToBEDFormat(args[1], args[2]);
			} else if (args.length == 5) {
				convertVCFsToBEDFormat(args[1], args[2], args[3], args[4]);
			}
		} else if (args[0].equals("svtkMatch")) {
			if (args.length == 4) {
				svtkMatch(args[1], args[2], args[3]);
			}
		} else if (args[0].equals("getCountsMatrixBuffered")) {
			if (args.length == 4) {
				getCountsMatrixBuffered(args[1], args[2], args[3]);
			} else if (args.length == 5) {
				getCountsMatrixBuffered(args[1], args[2], args[3], Integer.parseInt(args[4]));
			}
		} else if (args[0].equals("getPerSampleMetrics")) {
			if (args.length == 4) {
				getPerSampleMetrics(args[1], args[2], args[3], true);
			}
		} else if (args[0].equals("labelMedianQS")) {
			if (args.length == 5) {
				labelMedianQS(args[1], args[2], args[3], args[4]);
			}
		} else if (args[0].equals("jointCallVCF")) {
			if (args.length == 7) {
				jointCallVCF(args[1], args[2], args[3], args[4], args[5], args[6]);
			} else if (args.length == 5) {
				jointCallVCF(args[1], args[2], args[3], args[4]);
			}
		} else if (args[0].equals("annotateWithGenes")) {
			if (args[1].contains("strict")) {
				ArrayList<Gene> genes = parseGTFFile(args[3]);
				if (args.length == 6) {
					annotateGenesPercentBased(genes, args[2], args[4], args[5]);
				} else if (args.length == 5) {
					annotateGenesPercentBased(genes, args[2], args[4], "1,2,3,5");
				}

			} else if (args[1].contains("any")) {
				ArrayList<Gene> genes = parseGTFFile(args[3]);
				if (args.length == 6) {
					annotateGenesAnyOverlap(genes, args[2], args[4], args[5]);
				} else if (args.length == 5) {
					annotateGenesAnyOverlap(genes, args[2], args[4], "1,2,3,5");
				}
			}
		} else if (args[0].equals("condenseBedtoolsIntersect")) {
			condenseBedtoolsIntersect(args[1], args[2], args[3], args[4], args[5]);
		} else if (args[0].equals("defragment")) {
			this.defragment(args[1], args[2]);
		} else if (args[0].equals("subsetAnnotations")) {
			ArrayList<String> annotationSubsets = new ArrayList<>();
			for (int i = 4; i < args.length; i++) {
				annotationSubsets.add(args[i]);
			}
			this.subsetAnnotations(args[1], args[2], args[3], args[4], annotationSubsets);
		} else if (args[0].equals("validateSubsetAnnotations")) {
			ArrayList<String> annotationSubsets = new ArrayList<>();
			for (int i = 2; i < args.length; i++) {
				annotationSubsets.add(args[i]);
			}
			// this.subsetAnnotations(args[1], args[2], args[3], args[4],
			// annotationSubsets);
			this.validateSubsetAnnotations(args[1], annotationSubsets);
		} else if (args[0].equals("convertToEnsemble")) {
			this.convertToEnsemble(args[1], args[2], args[3], args[4]);

		} else if (args[0].equals("countExons")) {
			this.countExons(args[1], args[2], args[3], args[4]);
		} else {
			System.out.println("unknown command");
		}
	}

}
