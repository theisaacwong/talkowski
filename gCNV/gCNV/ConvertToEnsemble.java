package gCNV;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ConvertToEnsemble extends gCNVHelperTool {

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH, COLUMN_NAME, GTF_PATH};

	public ConvertToEnsemble(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException {
		String GCNV_INPUT = args.get(INPUT_PATH);
		String OUTPUT = args.get(OUTPUT_PATH);
		String geneColumnName = args.get(COLUMN_NAME);
		String gencodeGTF = args.get(GTF_PATH);
		this.convertToEnsemble(GCNV_INPUT, OUTPUT, geneColumnName, gencodeGTF);
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

}
