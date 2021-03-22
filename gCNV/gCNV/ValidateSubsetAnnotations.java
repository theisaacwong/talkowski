package gCNV;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Scanner;

public class ValidateSubsetAnnotations extends gCNVHelperTool {

	public static String[] inputs = new String[] {GENE_SET_FILE, GTF_PATH};
	
	public ValidateSubsetAnnotations(ArgParser args, String[] toolArgs) {
		super(args, toolArgs);
	}
	
	public ValidateSubsetAnnotations(ArgParser args) {
		super(args, inputs);
	}

	@Override
	public void run() throws IOException, InterruptedException {
		String gtfFile = args.get(GTF_PATH);
		
		String annotationSubset = args.get(GENE_SET_FILE);
		ArrayList<String> annotationSubsets = new ArrayList<>();
		annotationSubsets.add(annotationSubset);
		if(args.contains("extra")) {
			annotationSubsets.addAll(Arrays.asList(args.get("extra").split(",")));
		}
		
		this.validateSubsetAnnotations(gtfFile, annotationSubsets);
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

	

}
