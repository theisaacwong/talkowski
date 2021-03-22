package gCNV;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Stream;

public class ConvertVCFsToBEDFormat extends gCNVHelperTool {

	public ConvertVCFsToBEDFormat(String[] args) {
		super(args);
	}

	@Override
	public void run() throws IOException {
		String wd = args[1];
		String output = args[2];
		String prefixRegex = args[3];
		String suffixRegex = args[4];
		this.convertVCFsToBEDFormat(wd, output, prefixRegex, suffixRegex);
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
	
}
