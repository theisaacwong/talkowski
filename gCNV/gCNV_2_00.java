package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

/**
 * 
 * @author Isaac Wong
 * December 2, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski Lab
 * 
 * This tool is designed as a bridge between the various steps of gCNV
 *
 */
public class gCNV_2_00 {
	
	
	public String[] initializationArgs;
	public String date;
	
	public gCNV_2_00(String[] args) {
		initializationArgs = args;
		date = Calendar.getInstance().getTime().toString();
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		
		gCNV_2_00 g = new gCNV_2_00(args);
		System.out.println(g.toString());
		if(args[0].contains("-help") || args[0].contains("-h")) {
			System.out.println("-getBarcodeCounts [entityPath] [working-directory]");
			System.out.println("-getCountsMatrix [sourceFolder] [working-OUTPUT_PATH]");
			System.out.println("-getCountsMatrixParallel [sourceFolder] [working-OUTPUT_PATH]");
			System.out.println("-convertVCFsToBEDFormat [working-diectory] [output-path]");
			System.out.println("-svtkMatch [svtk_input] [svtk_output] [output_path]");
		} else if(args[0].equals("-getBarcodeCounts")) {
			g.getBarcodeCounts(args[1], args[2]);
		} else if(args[0].equals("-getCountsMatrix")) {
			g.getCountsMatrix(args[1], args[2]);
		} else if(args[0].equals("-getCountsMatrixParallel")) {
			g.getCountsMatrixParallel(args[1], args[2]);
		} else if(args[0].equals("-convertVCFsToBEDFormat")) {
			g.convertVCFsToBEDFormat(args[1], args[2]);
		} else if(args[0].equals("-svtkMatch")) {
			g.svtkMatch(args[1], args[2], args[3]);
		}
		
		//g.getBarcodeCounts("C:/Users/iwong/Documents/MGH/CMG_10_29_19/counting_crams_2019_11_12/sample_set_entity.tsv", "C:/Users/iwong/Documents/temp/");
		//g.getCountsMatrix("C:/Users/iwong/Documents/temp/", "C:/Users/iwong/Documents/temp/temp_count_mat.tsv");
		//g.getCountsMatrixParallel("C:/Users/iwong/Documents/temp/", "C:/Users/iwong/Documents/temp/temp_count_mat.tsv");
//		g.convertVCFsToBEDFormat("C:/Users/iwong/Documents/MGH/CMG_10_29_19/cohort_mode_2019_11_13/sample_set_entity.tsv", "C:/Users/iwong/Documents/MGH/CMG_10_29_19/cohort_mode_2019_11_13/", "svtk_input.bed");
//		g.svtk(svtk_input, svtk_output);
//		g.svtkMatch("C:/Users/iwong/Documents/MGH/CMG_10_29_19/cohort_mode_2019_11_13/svtk_input_java.bed", "C:/Users/iwong/Documents/MGH/CMG_10_29_19/cohort_mode_2019_11_13/svtk_output_java_testing.bed", "C:/Users/iwong/Documents/MGH/CMG_10_29_19/cohort_mode_2019_11_13/svtk_match_testing.tsv");
		
	}
	
	public String toString() {
		return String.join(" ", initializationArgs) + "\n" + date;
	}
	
	
	/**
	 * currently, I am hard-coding things, but in the future I can make the maps dynamically
	 * @param svtk_input
	 * @param svtk_output
	 * @param match_output
	 * @throws IOException
	 */
	public void svtkMatch(String svtk_input, String svtk_output, String match_output) throws IOException {
		DataFrame svtkInput = new DataFrame(svtk_input, true, "\\t", "@");
		DataFrame svtkOutput = new DataFrame(svtk_output, true, "\\t", "@");
		
		HashMap<String, String> CN_map = new HashMap<>();
		HashMap<String, String> GT_map = new HashMap<>();
		HashMap<String, String> NP_map = new HashMap<>();
		HashMap<String, String> QA_map = new HashMap<>();
		HashMap<String, String> QS_map = new HashMap<>();
		HashMap<String, String> QSE_map = new HashMap<>();
		HashMap<String, String> QSS_map = new HashMap<>();
		for(int i = 0; i < svtkInput.nrow(); i++) {
			CN_map.put(svtkInput.get("name", i), svtkInput.get("CN", i));
			GT_map.put(svtkInput.get("name", i), svtkInput.get("GT", i));
			NP_map.put(svtkInput.get("name", i), svtkInput.get("NP", i));
			QA_map.put(svtkInput.get("name", i), svtkInput.get("QA", i));
			QS_map.put(svtkInput.get("name", i), svtkInput.get("QS", i));
			QSE_map.put(svtkInput.get("name", i), svtkInput.get("QSE", i));
			QSS_map.put(svtkInput.get("name", i), svtkInput.get("QSS", i));
		}
		
		int oldnrow = svtkOutput.nrow();
		int counter = 0;
		for(int i = 0; i < svtkOutput.nrow(); i++) {
			if(svtkOutput.get("call_name", i).contains(",")) {
				String[] call_names = svtkOutput.get("call_name", i).split(",");
				svtkOutput.get(i)[svtkOutput.columnMapping.get("call_name")] = call_names[0];
				for(int k = 1; k < call_names.length; k++) {
					String[] newRow = new String[svtkOutput.get(i).length];
					for(int j = 0; j < newRow.length; j++) {
						newRow[j] = svtkOutput.get(i, j);
					}
					newRow[svtkOutput.columnMapping.get("call_name")] = call_names[k];
					svtkOutput.df.add(newRow); // shouldn't need to i--
					counter++;
				}
			}
		}
		System.out.println("oldnrow: " + oldnrow + "\tnewnrow: " + svtkOutput.nrow() + "\trowsadded: " + counter);
		
		ArrayList<String> CN_newField = new ArrayList<>();
		ArrayList<String> GT_newField = new ArrayList<>();
		ArrayList<String> NP_newField = new ArrayList<>();
		ArrayList<String> QA_newField = new ArrayList<>();
		ArrayList<String> QS_newField = new ArrayList<>();
		ArrayList<String> QSE_newField = new ArrayList<>();
		ArrayList<String> QSS_newField = new ArrayList<>();
		
		for(int i = 0; i < svtkOutput.nrow(); i++){
			CN_newField.add(CN_map.get(svtkOutput.get("call_name", i)));
			GT_newField.add(GT_map.get(svtkOutput.get("call_name", i)));
			NP_newField.add(NP_map.get(svtkOutput.get("call_name", i)));
			QA_newField.add(QA_map.get(svtkOutput.get("call_name", i)));
			QS_newField.add(QS_map.get(svtkOutput.get("call_name", i)));
			QSE_newField.add(QSE_map.get(svtkOutput.get("call_name", i)));
			QSS_newField.add(QSS_map.get(svtkOutput.get("call_name", i)));
		}
		
		ArrayList<String> columnNames = new ArrayList<>();
		columnNames.add("CN");
		columnNames.add("GT");
		columnNames.add("NP");
		columnNames.add("QA");
		columnNames.add("QS");
		columnNames.add("QSE");
		columnNames.add("QSS");
		
		ArrayList<ArrayList<String>> columnValues = new ArrayList<>();
		columnValues.add(CN_newField);
		columnValues.add(GT_newField);
		columnValues.add(NP_newField);
		columnValues.add(QA_newField);
		columnValues.add(QS_newField);
		columnValues.add(QSE_newField);
		columnValues.add(QSS_newField);
		
		svtkOutput.addColumns(columnNames, columnValues);
		svtkOutput.writeFile(match_output, true);

	}
	
	/**
	 * work in progress, you probably should not use this
	 * all this does is call svtk
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void svtk(String svtk_input, String svtk_output) throws IOException, InterruptedException {
		Runtime r = Runtime.getRuntime();
		Process p = r.exec("svtk bedcluster " + svtk_input + " " + svtk_output); // linux
//		Process p = r.exec("bash -c 'svtk bedcluster  " + svtk_input + " " + svtk_output + "'"); // windows
		p.waitFor();
	}
	
	/**
	 * TODO: make run parallel??, currenty, it runs in < 1 minute, so not a priority
	 * @param entityPath
	 * @param membershipPath
	 * @param wd
	 * @throws IOException 
	 */
	public void convertVCFsToBEDFormat(String wd, String output) throws IOException {
//		DataFrame sampleSetMembership = new DataFrame(membershipPath, true, "\\t", "@");	
//		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
		
		File[] directories = new File(wd).listFiles(File::isDirectory);
		System.out.println(Arrays.toString(directories));
		System.out.println("n directories: " + directories.length);
		
		ArrayList<String> all_bed_paths = new ArrayList<>();
		
		
		for(int i = 0; i < directories.length; i++) {
			int cnvNameCounter = 1; // might need to be long
			//String currentCluster = sampleSetEntity.get("entity:sample_set_id", i);
			String currentCluster = directories[i].getAbsolutePath();
			
			ArrayList<Path> currentClusterVCFsPATH = new ArrayList<>();
			ArrayList<String> currentClusterVCFsNAME = new ArrayList<>();
			ArrayList<String> sampleNames = new ArrayList<>();
			
	        Path p = Paths.get(currentCluster);
	        final int maxDepth = 1;
	        Stream<Path> matches = Files.find(p, maxDepth, (path, basicFileAttributes) -> String.valueOf(path).endsWith(".vcf"));
	        matches.filter(s->s.getFileName().toString().contains("vcf")).forEach(currentClusterVCFsPATH::add);
	        matches.close();
	        for(Path fp : currentClusterVCFsPATH) {
	        	currentClusterVCFsNAME.add(fp.toAbsolutePath().toString());
	        	sampleNames.add(fp.getFileName().toString().replaceAll(".tsv", ""));
	        }
	        
	        ArrayList<DataFrame> VCFsArrayList = new ArrayList<>();
	        for(String vcfFile : currentClusterVCFsNAME) {
	        	//System.out.println(iter_c + " " + barcodeCountFile); iter_c++;
	        	VCFsArrayList.add(new DataFrame(vcfFile, true, "\\t", "##"));
	        }
	        
	        System.out.println(currentCluster + " " + 
	        		currentClusterVCFsPATH.size() + " " + 
	        		currentClusterVCFsNAME.size() + " " + 
	        		sampleNames.size() + " " + 
	        		VCFsArrayList.size());

	        // create an array list of to store the new column names, eg "GT, CN, NP, QA" etc
	        String[] newColumns = VCFsArrayList.get(0).get("FORMAT", 0).split(":"); // an array list of the new field names, eg 'GT', 'CN', etc
	        ArrayList<String> newColumnNames = new ArrayList<>();
        	for(int j = 0; j < newColumns.length; j++) {
        		newColumnNames.add(newColumns[j]);
        	}
        	
	        DataFrame fullClusterBed = new DataFrame();
	        
	        for(int k = 0; k < VCFsArrayList.size(); k++) {
	        	DataFrame currentVCF = new DataFrame();
	        	ArrayList<String> chr = VCFsArrayList.get(k).get("#CHROM");
	        	ArrayList<String> start = VCFsArrayList.get(k).get("POS");
	        	ArrayList<String> end = VCFsArrayList.get(k).get("ID");
	        	for(int j = 0; j < end.size(); j++) { // Lord help me
	        		end.set(j, end.get(j).split("_")[end.get(j).split("_").length-1]);
	        	}
	        	ArrayList<String > name = new ArrayList<>();
	        	for(int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
	        		name.add(directories[i].getName() + "_cnv_" + cnvNameCounter);
	        		cnvNameCounter++;
	        	}
	        	ArrayList<String> sample = new ArrayList<>();
	        	for(int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
	        		sample.add(sampleNames.get(k).replace("genotyped-segments-", "").replace(".vcf", ""));
	        	}
	        	
	        	// create a list of lists to store the new column values
	        	ArrayList<ArrayList<String>> newColumnValues = new ArrayList<>();
	        	// this is the raw string from the VCF, the last field name is the raw strings of interest, so length-1 is the right index
	        	ArrayList<String> columnValues = VCFsArrayList.get(k).get(
	        			VCFsArrayList.get(k).fieldNames[VCFsArrayList.get(k).fieldNames.length -1]);
	        	
	        	//sanity check
	        	if(columnValues.size() != VCFsArrayList.get(k).nrow()) System.out.println("ERROR 24824776");
	        	
	        	// add a list to store every new field we will be parsing/adding
	        	for(int j = 0; j < newColumnNames.size(); j++) {
	        		newColumnValues.add(new ArrayList<>());
	        	}
	        	
	        	// for every row in the current VCF
	        	for(int j = 0; j <  VCFsArrayList.get(k).nrow(); j++) {
	        		//for every new column
	        		for(int l = 0; l < newColumnNames.size(); l++) {
	        			// to the target column   add the value     in the split[] index
	        			newColumnValues.get(l).add(columnValues.get(j).split(":")[l]);
	        		}
	        	}
	        	
	        	ArrayList<String> columnNamesToAdd = new ArrayList<>();
	        	ArrayList<ArrayList<String>> columnValuesToAdd = new ArrayList<>();
	        	
	        	columnNamesToAdd.add("chr"); 	columnValuesToAdd.add(chr); //System.out.println(chr.size());
	        	columnNamesToAdd.add("start");	columnValuesToAdd.add(start);
	        	columnNamesToAdd.add("end");	columnValuesToAdd.add(end);
	        	columnNamesToAdd.add("name");	columnValuesToAdd.add(name);
	        	columnNamesToAdd.add("sample");	columnValuesToAdd.add(sample);
	        	columnNamesToAdd.addAll(newColumnNames);	columnValuesToAdd.addAll(newColumnValues);
	        	boolean temp_1 = currentVCF.addColumns(columnNamesToAdd, columnValuesToAdd);
	        	if(!temp_1) {
	        		System.out.println("eror 427291233: " + currentClusterVCFsNAME.get(k));
	        	}
	        	//System.out.println(temp_1);
	        	// label dups and dels
	        	ArrayList<String> svtype = new ArrayList<>();
	        	//VCFsArrayList.get(k).nrow() should now be equal to currentVCF?
	        	//System.out.println(currentVCF.columnMapping.toString());
	        	for(int j = 0; j < VCFsArrayList.get(k).nrow(); j++) {
	        		//System.out.println(currentVCF.get("chr").size());
	        		String currChr = currentVCF.get("chr").get(j);
	        		String svt = "NA_j";
	        		if(currChr.equals("chrX")) {
	        			svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 2 ? "DUP" : "DEL";
	        		} else if(currChr.equals("chrY")) {
	        			svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 1 ? "DUP" : "DEL";
	        		} else {
	        			svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 2 ? "DUP" : "DEL";
	        		}
	        		svtype.add(svt);
	        	}
	        	currentVCF.addColumn("svtype", svtype);
	        	
	        	//remove rows where GT=0, this is not the most efficient way to remove indexes from arraylist
	        	ArrayList<Integer> gtIsZero = new ArrayList<>();
	        	for(int j = currentVCF.nrow()-1; j >= 0; j--) {
	        		if(currentVCF.get("GT", j).equals("0")) {
	        			gtIsZero.add(j);
	        		}
	        	}
	        	Collections.sort(gtIsZero, Collections.reverseOrder());
	        	for(int j : gtIsZero) {
	        		currentVCF.df.remove(j);
	        	}
	        	
	        	if(fullClusterBed.df.size() == 0) {
	        		fullClusterBed = currentVCF;
	        	} else {
	        		boolean ctrl = fullClusterBed.rbind(currentVCF);
	        		if(ctrl == false) {
	        			System.out.println("error 23089438  " + k + " " + i);
	        			System.out.println();
	        		}
	        	}
	        }
	        
	        System.out.println(currentCluster  + "\\" + directories[i].getName() + ".bed");
	        fullClusterBed.writeBed(currentCluster  + "\\" + directories[i].getName() + ".java.bed");
	        
	        all_bed_paths.add(currentCluster  + "\\" + directories[i].getName() + ".java.bed");
		}
		
		DataFrame fullMergedBed = new DataFrame(all_bed_paths.get(0), true, "\\t", "@");
		for(int i = 1; i < all_bed_paths.size(); i++) {
			fullMergedBed.rbind(new DataFrame(all_bed_paths.get(i), true, "\\t", "@"));
		}
		System.out.println(wd + "svtk_input_java.bed");
		fullMergedBed.writeFile(wd  + "svtk_input_java.bed", true);
		
		
		
	}
	
	/**
	 * downloads vcfs from cohort and case mode
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void downloadSegmentsVCFs(String entityPath, String membershipPath, String wd) throws IOException, InterruptedException {
//		DataFrame sampleSetMembership = new DataFrame(membershipPath, true, "\\t", "@");	
		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
//		HashMap<String, String> all_bed_paths = new HashMap<>();
		for(int i = 0; i < sampleSetEntity.nrow(); i++) {
			String individualToDownload = sampleSetEntity.get("entity:sample_set_id", i);
			String pathInGoogleBucket = sampleSetEntity.get("segments_vcfs", i).replaceAll("\\[|\\]|\"", "").replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");
			
			if(files.length < 2) {continue;}
			
			File file = new File(wd + individualToDownload + "_temp_files_to_download.txt");
			BufferedWriter output = new BufferedWriter(new FileWriter(file));
			for(String filePath : files) {
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
			
	        //this is only tested on windows, and it works so far. ITS VERY DELICATE. 
	        // download vcfs
	        String[] command;
	        if(System.getProperty("os.name").contains("Windows")) {
	        	String[] dosCommand = {"bash", "-c" ,"'cat", wdUnix+individualToDownload+"_temp_files_to_download.txt", "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix, "'"};
	        	command = dosCommand;
	        } else {
	        	String[] unixCommand = {"cat", wdUnix+individualToDownload+"_temp_files_to_download.txt", "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix};
	        	command = unixCommand;
	        }
	        System.out.println(String.join(" ", command));
			try {
		        new ProcessBuilder(command).inheritIO().start().waitFor();
		    } catch(IOException e) {
		        e.printStackTrace();
		    }
			
			// unzip vcfs
	        if(System.getProperty("os.name").contains("Windows")) {
	        	String[] dosCommand = {"bash", "-c" ,"'gunzip -k ", wdUnix + "*gz'"};
	        	command = dosCommand;
	        } else {
	        	String[] unixCommand = {"gunzip", "-k", wdUnix + "*gz'"};
	        	command = unixCommand;
	        }
	        System.out.println(String.join(" ", command));
			try {
		        new ProcessBuilder(command).inheritIO().start().waitFor();
		    } catch(IOException e) {
		        e.printStackTrace();
		    }
		}
	}
	
	/**
	 * @deprecated - please use parallel version
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @throws IOException
	 */
	public void getCountsMatrix(String sourceFolder, String OUTPUT_PATH) throws IOException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();
		
		System.out.print("finding files.\t");
        Path p = Paths.get(sourceFolder);
        final int maxDepth = 10;
        Stream<Path> matches = Files.find(p, maxDepth, (path, basicFileAttributes) -> String.valueOf(path).endsWith(".tsv"));
        matches.filter(s->s.getFileName().toString().contains("counts")).forEach(barcodeCountsPaths::add);
        matches.close();
        for(Path fp : barcodeCountsPaths) {
        	barcodeCountsFiles.add(fp.toAbsolutePath().toString());
        	sampleNames.add(fp.getFileName().toString().replaceAll(".barcode.counts.tsv", ""));
        }
        System.out.println("found " + sampleNames.size() + " files");
        
        System.out.print("reading files. \t");
        ArrayList<ArrayList<String>> countsArrayList = new ArrayList<>();
        //int iter_c = 1;
        for(String barcodeCountFile : barcodeCountsFiles) {
        	//System.out.println(iter_c + " " + barcodeCountFile); iter_c++;
        	DataFrame countsDF = new DataFrame(barcodeCountFile, true, "\\t", "@");
        	countsArrayList.add(countsDF.get("COUNT"));
        }
        System.out.println("done reading files");
        
        
        
        System.out.print("generating counts matrix. \t");
        ArrayList<String> labels = new ArrayList<>();
        DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
        for(int i = 0; i < countsDF.nrow(); i++) {
        	labels.add(countsDF.get("CONTIG", i) + "_" + countsDF.get("START", i) + "_" + countsDF.get("END", i));
        }
        System.out.println("done generating counts matrix");

        System.out.println("sampleNames.size()\t" + sampleNames.size());
        System.out.println("countsArrayList.size()\t" +  countsArrayList.size());
        System.out.println("countsArrayList.get(0).size()\t" +  countsArrayList.get(0).size());
        System.out.println("labels.size()\t" + labels.size());
        System.out.println("writing counts matrix");

		File file = new File(OUTPUT_PATH);
        BufferedWriter output = new BufferedWriter(new FileWriter(file));
		
		for(int i = 0; i < labels.size(); i++) {
			output.write(labels.get(i));
			if(i != labels.size()-1) {
				output.write("\t");	
			}
		} 
		output.write("\n");
		
		for(int i = 0; i < sampleNames.size(); i++) {
			StringBuilder line = new StringBuilder();
			line.append(sampleNames.get(i) + "\t");
			for(int k = 0; k < countsArrayList.get(i).size(); k++) {
				line.append(countsArrayList.get(i).get(k));
				if(k != countsArrayList.size()-1) {
					line.append("\t");
				}
			}
			if(i != sampleNames.size()-1) {
				line.append("\n");
			}
			output.write(line.toString());
		}
		output.close();
	}
	
	/**
	 * Reads in barcode counts files in parallel, and writes a counts matrix
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void getCountsMatrixParallel(String sourceFolder, String OUTPUT_PATH) throws IOException, InterruptedException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();
		
		System.out.print("finding files.\t");
        Path p = Paths.get(sourceFolder);
        final int maxDepth = 10;
        Stream<Path> matches = Files.find(p, maxDepth, (path, basicFileAttributes) -> String.valueOf(path).endsWith(".tsv"));
        matches.filter(s->s.getFileName().toString().contains("counts")).forEach(barcodeCountsPaths::add);
        matches.close();
        for(Path fp : barcodeCountsPaths) {
        	barcodeCountsFiles.add(fp.toAbsolutePath().toString());
        	sampleNames.add(fp.getFileName().toString().replaceAll(".barcode.counts.tsv", ""));
        }
        System.out.println("found " + sampleNames.size() + " files");
        
        System.out.print("reading files. \t");
        
        List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
        Map<Integer, ArrayList<String>> doneReading = new ConcurrentHashMap<>();
        for(int i = 0; i < sampleNames.size(); i++) {toRead.add(i);}
        int N_THREADS =  Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {
					while(toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading.put(currentFile, countsDF.get("COUNT"));
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

		/*
		 * I could easily-er do the following 
		 * ArrayList<ArrayList<String>> countsArrayList = new ArrayList<>(doneReading.values());
		 * but I want the samples written in read-order for consistency, will be a bit slower  
		 */
		ArrayList<ArrayList<String>> countsArrayList = new ArrayList<>();
		for(int i = 0; i < doneReading.size(); i++) {
			countsArrayList.add(doneReading.get(i));
		}

        
        
        System.out.println("done reading files");
        
        
        
        System.out.print("generating counts matrix. \t");
        ArrayList<String> labels = new ArrayList<>();
        DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
        for(int i = 0; i < countsDF.nrow(); i++) {
        	labels.add(countsDF.get("CONTIG", i) + "_" + countsDF.get("START", i) + "_" + countsDF.get("END", i));
        }
        System.out.println("done generating counts matrix");

        System.out.println("sampleNames.size()\t" + sampleNames.size());
        System.out.println("countsArrayList.size()\t" +  countsArrayList.size());
        System.out.println("countsArrayList.get(0).size()\t" +  countsArrayList.get(0).size());
        System.out.println("labels.size()\t" + labels.size());
        System.out.println("writing counts matrix");

		File file = new File(OUTPUT_PATH);
        BufferedWriter output = new BufferedWriter(new FileWriter(file));
		
		for(int i = 0; i < labels.size(); i++) {
			output.write(labels.get(i));
			if(i != labels.size()-1) {
				output.write("\t");	
			}
		} 
		output.write("\n");
		
		for(int i = 0; i < sampleNames.size(); i++) {
			StringBuilder line = new StringBuilder();
			line.append(sampleNames.get(i) + "\t");
			for(int k = 0; k < countsArrayList.get(i).size(); k++) {
				line.append(countsArrayList.get(i).get(k));
				if(k != countsArrayList.size()-1) {
					line.append("\t");
				}
			}
			if(i != sampleNames.size()-1) {
				line.append("\n");
			}
			output.write(line.toString());
		}
		output.close();
	}

	/**
	 * downloads the barcode counts files from the entity path file
	 * @param entityPath - the full path to the entity file
	 * @param wd - the working directory
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getBarcodeCounts(String entityPath, String wd) throws IOException, InterruptedException {
		DataFrame entityDF = new DataFrame(entityPath, true, "\\t", "@");
		System.out.println(entityDF.nrow());
		System.out.println(entityDF.columnMapping.toString());
		for(int i = 3; i < entityDF.nrow(); i++) {
			String individualToDownload = entityDF.get("entity:sample_set_id", i);
			String pathInGoogleBucket = entityDF.get("output_counts_barcode", i).replaceAll("\\[|\\]|\"", "").replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");
			
			if(files.length < 2) {continue;}
			
			File file = new File(wd + individualToDownload + "_temp_files_to_download.txt");
			BufferedWriter output = new BufferedWriter(new FileWriter(file));
			for(String filePath : files) {
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
			
	        //this is only tested on windows, and it works so far. ITS VERY DELICATE. 
	        String[] command;
	        if(System.getProperty("os.name").contains("Windows")) {
	        	String[] dosCommand = {"bash", "-c" ,"'cat", wdUnix+individualToDownload+"_temp_files_to_download.txt", "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix, "'"};
	        	command = dosCommand;
	        } else {
	        	String[] unixCommand = {"cat", wdUnix+individualToDownload+"_temp_files_to_download.txt", "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix};
	        	command = unixCommand;
	        }
	        System.out.println(String.join(" ", command));
	        
			try {
		        new ProcessBuilder(command).inheritIO().start().waitFor();
		    } catch(IOException e) {
		        e.printStackTrace();
		    }
		}
	}
	
	
	
	

}
