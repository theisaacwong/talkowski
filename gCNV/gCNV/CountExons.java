package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CountExons extends gCNVHelperTool {

	public CountExons(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException {
		String GCNV_INPUT = args[1];
		String genesColumnName = args[2];
		String gencodeGTF = args[3];
		String output = args[4];
		this.countExons(GCNV_INPUT, genesColumnName, gencodeGTF, output);
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
				nExonsColumn.add("0");
				totalExons.add("0");
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

}
