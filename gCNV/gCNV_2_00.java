package gCNV;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
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
import java.util.List;
import java.util.Map;
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
 * @since December 2, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski Lab
 * 
 * This tool is designed as a bridge between the various steps of gCNV
 * 
 * Contact Information
 * website: https://talkowski.mgh.harvard.edu/
 * email: iwong@broadinstitute.org
 * github: https://github.com/theisaacwong/
 *
 */
public class gCNV_2_00 {
	
	
	public String[] initializationArgs;
	public String date;
	
	public final String CHR = "CHR";
	public final String START = "START";
	public final String END = "END";
	public static final String VERSION = "2.3.5";
	
	public gCNV_2_00(String[] args) {
		initializationArgs = args;
		date = Calendar.getInstance().getTime().toString();
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		
		gCNV_2_00 g = new gCNV_2_00(args);
		System.out.println(g.toString());
		System.out.println("version " + VERSION);
		try{g.checkVersion();} catch(Exception e) {}
		System.out.println("Java version: " + System.getProperty("java.version"));
		System.out.println("Heap Size: " + getHeapSize());
		
		
		System.out.println();
		if(args.length==0 || args[0].contains("-help") || args[0].contains("-h")) {
			System.out.println("java -jar gCNV_helper.jar [Command] [required argument(s)] {optional arguement(s)}");
			System.out.println();
				System.out.println("\tgetBarcodeCounts [entityPath] [working-directory] {counts-field-name}");
					System.out.println("\t\tDownload the read count files specified in the entity file");
					System.out.println("\t\t[entityPath] - Full path to the entity file. eg: '/home/gCNV/sample_set_entity.tsv'");
					System.out.println("\t\t[working-directory] - Directory to download files to. eg '/home/gCNV/'");
					System.out.println("\t\t{counts-field-name} - optional - The name of the column in the entity file containing the counts paths. eg 'output_counts_barcode'");
					System.out.println();
				System.out.println("\tgetCountsMatrix [sourceFolder] [OUTPUT_PATH] {regex}");
					System.out.println("\t\tRead in all read count files and generate a matrix file");
					System.out.println("\t\t[sourceFolder] - Directory where read count files are located in, files can be in sub-directories.");
					System.out.println("\t\t[OUTPUT_PATH] -  The full output path where the matrix file will be written to");
					System.out.println("\t\t{regex} - optional - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'");
					System.out.println();
				System.out.println("\tgetCountsMatrixBuffered [sourceFolder] [OUTPUT_PATH] [regex] {buffer-size}");
					System.out.println("\t\tRead in all read count files and generate a matrix file");
					System.out.println("\t\t[sourceFolder] - Directory where read count files are located in, files can be in sub-directories.");
					System.out.println("\t\t[OUTPUT_PATH] -  The full output path where the matrix file will be written to");
					System.out.println("\t\t[regex] - The regex suffix used to identify counts files. eg '.barcode.counts.tsv'");
					System.out.println("\t\t{buffer-size} - number of lines to store in memory for each thread before writing");
					System.out.println();
				System.out.println("\tdownloadSegmentsVCFs [entityPath] [working-directory] {column-name}");
					System.out.println("\t\tDownload the VCF file outputs from running the main gCNV algorithm");
					System.out.println("\t\t[entityPath] - Full path to the entity file. eg: '/home/gCNV/sample_set_entity.tsv'");
					System.out.println("\t\t[working-directory] - Directory to download files to. eg '/home/gCNV/'");
					System.out.println("\t\t{column-name} - optional - The name of the column in the entity file containing the counts paths. eg 'segments_vcfs'");
					System.out.println();
				System.out.println("\tconvertVCFsToBEDFormat [working-directory] [output-path] {prefix-regex} {suffix-regex}");
					System.out.println("\t\tConvert the gCNV VCF output files to BED format for svtk bedlcuster input");
					System.out.println("\t\t[working-directory] - Directory where VCF files are located in, files can be in sub-directories.");
					System.out.println("\t\t[ouput-path] - The output path for the final consolidated BED file");
					System.out.println("\t\t{prefix-regex} - prefix to trim from file name, eg 'genotyped-segments-'");
					System.out.println("\t\t{suffix-regex} - suffix used to identify VCF files, used also to trim from file name. eg '.vcf'");
					System.out.println();
				System.out.println("\tsvtkMatch [svtk_input] [svtk_output] [output_path]");
					System.out.println("\t\tMatch up the gCNV meta data with the svtk bedcluster meta data and write to file");
					System.out.println("\t\t[svtk_input] - The BED file that was given to svtk bedcluster");
					System.out.println("\t\t[svtk_output] - The output file from svtk bedcluster");
					System.out.println("\t\t[output_path] - The full path to write the output file to.");
					System.out.println();
		} else if(args[0].equals("getBarcodeCounts")) {
			if(args.length == 3) {
				g.getBarcodeCounts(args[1], args[2]);	
			} else if(args.length == 4) {
				g.getBarcodeCounts(args[1], args[2], args[3]);
			}
		} else if(args[0].equals("getCountsMatrix")) {
			if(args.length == 3) {
				g.getCountsMatrix(args[1], args[2]);	
			} else if(args.length == 4) {
				g.getCountsMatrix(args[1], args[2], args[3]);
			}
		} else if(args[0].equals("downloadSegmentsVCFs")) {
			if(args.length == 3) {
				g.downloadSegmentsVCFs(args[1], args[2]);	
			} else if(args.length == 4) {
				g.downloadSegmentsVCFs(args[1], args[2], args[3]);
			}
		} else if(args[0].equals("convertVCFsToBEDFormat")) {
			if(args.length == 3) {
				g.convertVCFsToBEDFormat(args[1], args[2]);
			} else if(args.length == 5) {
				g.convertVCFsToBEDFormat(args[1], args[2], args[3], args[4]);
			}
		} else if(args[0].equals("svtkMatch")) {
			if(args.length == 4) {
				g.svtkMatch(args[1], args[2], args[3]);	
			}
		} else if(args[0].equals("getCountsMatrixBuffered")) {
			if(args.length == 4) {
				g.getCountsMatrixBuffered(args[1], args[2], args[3]);
			} else if(args.length == 5) {
				g.getCountsMatrixBuffered(args[1], args[2], args[3], Integer.parseInt(args[4]));
			}
		} 
	}
	
	public void checkVersion() throws IOException {
		//https://stackoverflow.com/questions/9489726/get-raw-text-from-html
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
		if(versionString.equals(VERSION)) {
			System.out.println("(Most current version)");
		} else {
			System.out.println("Waring. You are currently running version " + VERSION + ". The most recent version is " + versionString + ". Updating is strongly recommended.");
		}
	}
	
	public void defrag(String match_output, String defrag_output) {
		/*
		 * 1. give each CNV call its own unique ID number, have a field 
		 */
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
		//svtk will sometimes have two callnames in on index, so duplicate the line for each element in the bin
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
		
		// create arraylists for all the new fields to be added
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
		p.waitFor();
	}
	
	/**
	 * A wrapper to call the full convertVCFsToBEDFormat() using generic arguements
	 * @param wd
	 * @param output
	 * @throws IOException
	 */
	public void convertVCFsToBEDFormat(String wd, String output) throws IOException {
		convertVCFsToBEDFormat(wd, output, "genotyped-segments-", ".vcf");
	}
	
	/**
	 * 
	 * @param wd - the working directory where VCFs are stored. VCFs can be in sub-directories.
	 * @param output - The output path for the final consolidated BED file
	 * @param prefixRegex - prefix to trim from file name, eg "genotyped-segments-"
	 * @param suffixRegex - suffix used to identify VCF files, used also to trim from file name. eg ".vcf"
	 * @throws IOException
	 */
	public void convertVCFsToBEDFormat(String wd, String output, String prefixRegex, String suffixRegex) throws IOException {
		File[] directories = new File(wd).listFiles(File::isDirectory);
		System.out.println(Arrays.toString(directories));
		System.out.println("n directories: " + directories.length);
		
		ArrayList<String> all_bed_paths = new ArrayList<>();
		
		// look through all directories in current working directory
		for(int i = 0; i < directories.length; i++) {
			int cnvNameCounter = 1; // might need to be long
			String currentCluster = directories[i].getAbsolutePath();
			
			ArrayList<Path> currentClusterVCFsPATH = new ArrayList<>();
			ArrayList<String> currentClusterVCFsNAME = new ArrayList<>();
			ArrayList<String> sampleNames = new ArrayList<>();
			
			// get a path to all the vcf files, similar to 'ls wd/*vcf'
	        Path p = Paths.get(currentCluster);
	        final int maxDepth = 1;
	        Stream<Path> matches = Files.find(p, maxDepth, (path, basicFileAttributes) -> String.valueOf(path).endsWith(suffixRegex));
	        matches.filter(s->s.getFileName().toString().endsWith(suffixRegex)).forEach(currentClusterVCFsPATH::add);
	        matches.close();
	        for(Path fp : currentClusterVCFsPATH) {
	        	currentClusterVCFsNAME.add(fp.toAbsolutePath().toString());
	        	sampleNames.add(fp.getFileName().toString().replaceAll(suffixRegex, ""));
	        }
	        
	        // if there are no vcfs in this folder, move on to the next folder
	        if(currentClusterVCFsPATH.size() < 2) {
	        	continue;
	        }
	        
	        // an arraylist containing vcf dataframes
	        ArrayList<DataFrame> VCFsArrayList = new ArrayList<>();
	        for(String vcfFile : currentClusterVCFsNAME) {
	        	VCFsArrayList.add(new DataFrame(vcfFile, true, "\\t", "##"));
	        }
	        

	        // create an array list of to store the new column names, eg "GT, CN, NP, QA" etc
	        String[] newColumns = VCFsArrayList.get(0).get("FORMAT", 0).split(":"); // an array list of the new field names, eg 'GT', 'CN', etc
	        ArrayList<String> newColumnNames = new ArrayList<>();
        	for(int j = 0; j < newColumns.length; j++) {
        		newColumnNames.add(newColumns[j]);
        	}
        	
        	// to store the full merged bed file
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
	        		sample.add(sampleNames.get(k).replace(prefixRegex, "").replace(suffixRegex, ""));
	        	}
	        	
	        	// create a list of lists to store the new column values
	        	ArrayList<ArrayList<String>> newColumnValues = new ArrayList<>();
	        	// this is the raw string from the VCF, the last field name is the raw strings of interest, so length-1 is the right index
	        	ArrayList<String> columnValues = VCFsArrayList.get(k).get(
	        			VCFsArrayList.get(k).fieldNames[VCFsArrayList.get(k).fieldNames.length -1]);
	        	
	        	//sanity check
	        	if(columnValues.size() != VCFsArrayList.get(k).nrow()) System.out.println("ERROR 2482469776");
	        	
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
	        	
	        	columnNamesToAdd.add("chr"); 	columnValuesToAdd.add(chr); 
	        	columnNamesToAdd.add("start");	columnValuesToAdd.add(start);
	        	columnNamesToAdd.add("end");	columnValuesToAdd.add(end);
	        	columnNamesToAdd.add("name");	columnValuesToAdd.add(name);
	        	columnNamesToAdd.add("sample");	columnValuesToAdd.add(sample);
	        	columnNamesToAdd.addAll(newColumnNames);	columnValuesToAdd.addAll(newColumnValues);
	        	boolean temp_1 = currentVCF.addColumns(columnNamesToAdd, columnValuesToAdd);
	        	if(!temp_1) {
	        		System.out.println("eror 42721191233: " + currentClusterVCFsNAME.get(k));
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
	        		if(currChr.contains("X") || currChr.contains("x")) {
	        			svt = Integer.parseInt(currentVCF.get("CN").get(j)) >= 2 ? "DUP" : "DEL";
	        		} else if(currChr.contains("Y") || currChr.contains("y")) {
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
		System.out.println(wd + output);
		fullMergedBed.writeFile(wd  + output, true);
		
		
		
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
	 * @param entityPath - full path to the entity file, eg "sample_set_entity.tsv"
	 * @param wd - working directory, where the files should be downloaded to
	 * @param segments_vcfs_columnName - the name of the column in the entity file containing the segment_vcfs, eg "segments_vcfs"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void downloadSegmentsVCFs(String entityPath, String wd, String segments_vcfs_columnName) throws IOException, InterruptedException {
		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
		for(int i = 0; i < sampleSetEntity.nrow(); i++) {
			String individualToDownload = sampleSetEntity.get("entity:sample_set_id", i);
			String pathInGoogleBucket = sampleSetEntity.get(segments_vcfs_columnName, i).replaceAll("\\[|\\]|\"", "").replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");
			
			if(pathInGoogleBucket.trim().equals("") || (files.length==1 && files[0].trim().equals("")) || files.length==0) {continue;}
			
			String temp_suffix = "_files_to_download.txt";
			
			File file = new File(wd + individualToDownload + temp_suffix);
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
			
	        downloadFiles(wdUnix + individualToDownload + temp_suffix, pathToDownloadToUnix);
	        
//	        //this is only tested on windows, and it works so far. ITS VERY DELICATE. 
	        String[] command;
			
			// unzip vcfs, working as of 12/3/19
	        if(System.getProperty("os.name").contains("Windows")) {
	        	String[] dosCommand = {"bash", "-c" ,"'gunzip", "-k", pathToDownloadToUnix + "/*'"};
	        	command = dosCommand;
	        } else {
	        	String[] unixCommand = {"gunzip", "-k" + pathToDownloadToUnix + "/*'"};
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
	 * this is just a wrapper for the buffered version, which sets the buffer size to 20
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @param countsRegex
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrixBuffered(String sourceFolder, String OUTPUT_PATH, String countsRegex) throws IOException, InterruptedException {
		getCountsMatrixBuffered(sourceFolder, OUTPUT_PATH, countsRegex, 20);
	}
	/**
	 * write buffer after every n files is read per thread, uses less memory as less count data is held in memory at any given time
	 * however, uses more I/O, which will be slower overall
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @param countsRegex
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrixBuffered(String sourceFolder, String OUTPUT_PATH, String countsRegex, int BUFFER_SIZE) throws IOException, InterruptedException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();
		
		System.out.print("finding files.\t");
        Path p = Paths.get(sourceFolder);
        final int maxDepth = 10;
        Stream<Path> matches = Files.find(p, maxDepth, (path, basicFileAttributes) -> String.valueOf(path).endsWith(countsRegex));
        matches.filter(s->s.getFileName().toString().contains(countsRegex)).forEach(barcodeCountsPaths::add);
        matches.close();
        for(Path fp : barcodeCountsPaths) {
        	barcodeCountsFiles.add(fp.toAbsolutePath().toString());
        	sampleNames.add(fp.getFileName().toString().replaceAll(countsRegex, ""));
        }
        System.out.println("found " + sampleNames.size() + " files");
        
        System.out.println("reading files. \t");
        

        
        // toRead functions as a synchronized queue so each thread knows which file to read next
        List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
        for(int i = 0; i < sampleNames.size(); i++) {toRead.add(i);} 
        int N_THREADS =  Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		int totalNFiles = toRead.size();
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run(){
					File file = new File(OUTPUT_PATH + "_" + Thread.currentThread().getId());
			        BufferedWriter output = null;
					try {output = new BufferedWriter(new FileWriter(file));} catch (IOException e1) {e1.printStackTrace();}
			        
					HashMap<String, String> doneReading = new HashMap<>(); // map of sample name to line string
					
			        int bufferCounter = 0;
			        
			        ArrayList<String> labels = new ArrayList<>();
			        DataFrame labelsDF = null;
					try {labelsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");} catch (IOException e1) {e1.printStackTrace();}
			        for(int i = 0; i < labelsDF.nrow(); i++) {
			        	labels.add(labelsDF.get("CONTIG", i) + "_" + labelsDF.get("START", i) + "_" + labelsDF.get("END", i));
			        }
			        labelsDF = null;
			        
			        for(int i = 0; i < labels.size(); i++) {
						try {output.write(labels.get(i));} catch (IOException e) {e.printStackTrace();}
						if(i != labels.size()-1) {
							try {output.write("\t");} catch (IOException e) {e.printStackTrace();}	
						}
					} 
					try {output.write("\n");} catch (IOException e1) {e1.printStackTrace();}
			        
					while(toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading.put(sampleNames.get(currentFile), String.join("\t", countsDF.get("COUNT")));
							progressPercentage(totalNFiles - toRead.size(), totalNFiles, barcodeCountsFiles.get(currentFile));
							countsDF = null; // java gc is a fickle mistress
							
							if(bufferCounter >= BUFFER_SIZE) {
								StringBuilder line = new StringBuilder();
								for(String snKey : doneReading.keySet()) {
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
					for(String snKey : doneReading.keySet()) {
						line.append(snKey); // sample name key
						line.append("\t");
						line.append(doneReading.get(snKey));
						line.append("\n");
					}
					try {output.write(line.toString());} catch (IOException e2) {e2.printStackTrace();}
					doneReading = new HashMap<>();
					bufferCounter = 0;
			
					try {output.close();} catch (IOException e) {e.printStackTrace();}
				}
			});
		}
		exServer.shutdown();
		exServer.awaitTermination(6, TimeUnit.DAYS);

	}
	
	/**
	 * this is just a wrapper for getCountsMatrix() which sets the default regex to ".barcode.counts.tsv"
	 * @param sourceFolder
	 * @param OUTPUT_PATH
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getCountsMatrix(String sourceFolder, String OUTPUT_PATH) throws IOException, InterruptedException {
		getCountsMatrix(sourceFolder, OUTPUT_PATH, ".barcode.counts.tsv");
	}
	/**
	 * This is an improvement upon the previous version as it now more accurately forces GC to 
	 * remove the DF after its been read and only keep the important field
	 * Reads in barcode counts files in parallel, and writes a counts matrix
	 * @param sourceFolder - directory where files are located in, files can be in sub-directories. 
	 * @param OUTPUT_PATH - the full output path where the matrix file will be written to
	 * @param countsRegex - the regex suffix to identify counts files
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void getCountsMatrix(String sourceFolder, String OUTPUT_PATH, String countsRegex) throws IOException, InterruptedException {
		ArrayList<Path> barcodeCountsPaths = new ArrayList<>();
		ArrayList<String> barcodeCountsFiles = new ArrayList<>();
		ArrayList<String> sampleNames = new ArrayList<>();
		
		System.out.print("finding files.\t");
        Path p = Paths.get(sourceFolder);
        final int maxDepth = 10;
        Stream<Path> matches = Files.find(p, maxDepth, (path, basicFileAttributes) -> String.valueOf(path).endsWith(countsRegex));
        matches.filter(s->s.getFileName().toString().contains(countsRegex)).forEach(barcodeCountsPaths::add);
        matches.close();
        for(Path fp : barcodeCountsPaths) {
        	barcodeCountsFiles.add(fp.toAbsolutePath().toString());
        	sampleNames.add(fp.getFileName().toString().replaceAll(countsRegex, ""));
        }
        System.out.println("found " + sampleNames.size() + " files");
        
        System.out.println("reading files. \t");
        
        // toRead functions as a synchronized queue so each thread knows which file to read next
        // each dataframe is then mapped to a unique index number so that the arraylist of dataframes
        // can be assembled again in order even though each one is read out of order
        List<Integer> toRead = Collections.synchronizedList(new ArrayList<Integer>());
        Map<Integer, String> doneReading2 = new ConcurrentHashMap<>();
        for(int i = 0; i < sampleNames.size(); i++) {toRead.add(i);} 
        int N_THREADS =  Runtime.getRuntime().availableProcessors();
		ExecutorService exServer = Executors.newFixedThreadPool(N_THREADS);
		int totalNFiles = toRead.size();
		for (int i = 0; i < N_THREADS; i++) {
			exServer.execute(new Runnable() {
				@Override
				public void run() {
					while(toRead.size() > 0) {
						int currentFile = toRead.remove(0);
						try {
							DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(currentFile), true, "\\t", "@");
							doneReading2.put(currentFile, String.join("\t", countsDF.get("COUNT")));
							progressPercentage(totalNFiles - toRead.size(), totalNFiles, barcodeCountsFiles.get(currentFile));
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
		for(int i = 0; i < doneReading2.size(); i++) {
			countsArrayList.add(doneReading2.get(i));
		}
        
        System.out.println("done reading files");
        
        // 'labels' is the exon labels, / column names
        System.out.print("generating counts matrix. \t");
        ArrayList<String> labels = new ArrayList<>();
        DataFrame countsDF = new DataFrame(barcodeCountsFiles.get(0), true, "\\t", "@");
        for(int i = 0; i < countsDF.nrow(); i++) {
        	labels.add(countsDF.get("CONTIG", i) + "_" + countsDF.get("START", i) + "_" + countsDF.get("END", i));
        }
        System.out.println("done generating counts matrix");

        // for sanity checking, number should all appropriately match
        System.out.println("sampleNames.size()\t" + sampleNames.size());
        System.out.println("countsArrayList.size()\t" +  countsArrayList.size());
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
		
		if(countsArrayList.size() != sampleNames.size()) {
			System.out.println("error 51");
		}
		
		// written in such a way as to be readable in R, adding rownames for the sample ID
		for(int i = 0; i < sampleNames.size(); i++) {
			StringBuilder line = new StringBuilder();
			line.append(sampleNames.get(i) + "\t");
			line.append(countsArrayList.get(i));
			if(i != sampleNames.size()-1) {
				line.append("\n");
			}
			output.write(line.toString());
		}
		output.close();
	}
	
	
	/**
	 * a wrapper for the full getBarcodeCounts() with generic args
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
	 * @param entityPath - full path to the entity file, eg "/documents/sample_set_entity.tsv"
	 * @param wd - the directory to download the counts files to
	 * @param output_counts_barcode_regex - the name of the column containing the counts files, eg "output_counts_barcode"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void getBarcodeCounts(String entityPath, String wd, String output_counts_barcode_regex) throws IOException, InterruptedException {
		DataFrame entityDF = new DataFrame(entityPath, true, "\\t", "@");
		System.out.println(entityDF.nrow());
		System.out.println(entityDF.columnMapping.toString());
		for(int i = 0; i < entityDF.nrow(); i++) {
			String individualToDownload = entityDF.get("entity:sample_set_id", i);
			String pathInGoogleBucket = entityDF.get(output_counts_barcode_regex, i).replaceAll("\\[|\\]|\"", "").replaceAll(",", " ");
			String[] files = pathInGoogleBucket.split(" ");
			
			if(files.length < 2) {continue;}
			
			String temp_suffix = "_temp_files_to_download.txt";
			
			File file = new File(wd + individualToDownload + temp_suffix);
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
			
	        downloadFiles(wdUnix+individualToDownload+temp_suffix, pathToDownloadToUnix);

		}
	}
	
	/**
	 * takes in a file containing a list of files to download. downloads those files. 
	 * this has only been tested on windows, it is VERY DELICATE
	 * @param filesFile
	 * @throws InterruptedException 
	 */
	public static void downloadFiles(String filesFile, String pathToDownloadToUnix) throws InterruptedException {
        String[] command;
        if(System.getProperty("os.name").contains("Windows")) {
        	String[] dosCommand = {"bash", "-c" ,"'cat", filesFile, "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix, "'"};
        	command = dosCommand;
        } else {
        	System.out.println("WARNING: this download functionality has not been tested on macOS, please report any issues");
        	String[] unixCommand = {"cat", filesFile, "|", "gsutil", "-m", "cp", "-I", pathToDownloadToUnix};
        	command = unixCommand;
        }
        System.out.println(String.join(" ", command));
        
		try {
	        new ProcessBuilder(command).inheritIO().start().waitFor();
	    } catch(IOException e) {
	        e.printStackTrace();
	    }
	}
	
	/**
	 * https://stackoverflow.com/questions/852665/command-line-progress-bar-in-java
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
	    System.out.print("\r" + bareDone + bareRemain + " " + remainProcent * 10 + "%  " + remain + " / " + total + "  " + message);
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
		if (bytes < 1024) return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp-1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}
	
	
	
	
	

}
