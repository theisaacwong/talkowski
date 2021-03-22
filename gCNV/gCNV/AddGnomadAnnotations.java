package gCNV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class AddGnomadAnnotations extends gCNVHelperTool{

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH, GTF_PATH, COLUMN_NAME};
	
	public AddGnomadAnnotations(ArgParser args) {
		super(args, inputs);
	}


	@Override
	public void run() throws IOException {
		String GCNV_INPUT = args.get(INPUT_PATH);
		String OUTPUT = args.get(OUTPUT_PATH);
		String gencodeGTF = args.get(GTF_PATH);
		String geneColumnName = args.get(COLUMN_NAME);
		this.addGnomadAnnotations(GCNV_INPUT, OUTPUT, gencodeGTF, geneColumnName);
	}
	
	public void addGnomadAnnotations(String GCNV_INPUT, String OUTPUT, String gencodeGTF, String geneColumnName) throws IOException {
		ArrayList<Gene> geneList = parseGTFFile(gencodeGTF);
		HashMap<String, Gene> genes = new HashMap<>();
		for(int i = 0; i < geneList.size(); i++) {
			genes.put(geneList.get(i).name, geneList.get(i));
		}

		print("readint callset");
		DataFrame gcnv = new DataFrame(GCNV_INPUT, true, "\\t", "#");

		print("writing annotations");
		ArrayList<String> annotationColumn = new ArrayList<>();

		for(int i = 0; i < gcnv.nrow(); i++) {
			String[] currentGenes = gcnv.get(geneColumnName, i).split(",");
			ArrayList<String> annotations = new ArrayList<>();

			int gStart = Integer.parseInt(gcnv.get("start", i));
			int gEnd = Integer.parseInt(gcnv.get("end", i));
			String gSvtype = gcnv.get("svtype", i);

			for(String gene : currentGenes) {
				if(gene.equals("None")) {
					annotations.add("None");
				} else {
					annotations.add(genes.get(gene).getGnomadSchemeAnnotation(gStart, gEnd, gSvtype));	
				}
			}

			annotationColumn.add(String.join(",", annotations));
		}


		gcnv.addColumn(geneColumnName + "_gnomAD_annotation", annotationColumn);
		gcnv.writeFile(OUTPUT, true);
	}
	
}
