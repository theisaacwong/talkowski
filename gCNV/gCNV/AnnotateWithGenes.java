package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class AnnotateWithGenes extends gCNVHelperTool {

	public AnnotateWithGenes(String[] args) {
		super(args);
	}
	
	@Override
	public void run() throws IOException {
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
			output.write(line + "\tgenes_strict_overlap" + "\n");
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

}
