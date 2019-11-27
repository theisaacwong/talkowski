package gCNV;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

public class gCNV_2_00 {
	
	
	public DataFrame svtk_input;
	public DataFrame svtk_output;
	public DataFrame gCNV_output;
	public String svtk_input_path;
	public String svtk_output_path;
	public String output_file_path;
	
	public gCNV_2_00(String svtk_in, String svtk_out, String out_file_path) throws IOException {
		svtk_input = new DataFrame(svtk_in, true, "\t", "#");
		svtk_output = new DataFrame(svtk_out, true, "\t", "#");
		gCNV_output = new DataFrame();
		
	}
	
	public gCNV_2_00() {
		
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		
//		gCNV_2_00 gCNV = new gCNV_2_00(args[0], args[1], args[2]);
//		gCNV.match();
//		gCNV.writeFile();
		
		gCNV_2_00 g = new gCNV_2_00();
		//g.getBarcodeCounts("C:/Users/iwong/Documents/MGH/CMG_10_29_19/counting_crams_2019_11_12/sample_set_entity.tsv", "C:/Users/iwong/Documents/temp/");
		//g.getCountsMatrix("C:/Users/iwong/Documents/temp/", "C:/Users/iwong/Documents/temp/temp_count_mat.tsv");
		g.getCountsMatrixParallel("C:/Users/iwong/Documents/temp/", "C:/Users/iwong/Documents/temp/temp_count_mat.tsv");
	}
	
	/**
	 * TODO: finish
	 * @param entityPath
	 * @param membershipPath
	 * @param wd
	 * @throws IOException 
	 */
	public void convertVCFsToBEDFormat(String entityPath, String membershipPath, String wd) throws IOException {
		DataFrame sampleSetMembership = new DataFrame(membershipPath, true, "\\t", "@");	
		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
		HashMap<String, String> all_bed_paths = new HashMap<>();
		for(int i = 0; i < sampleSetEntity.nrow(); i++) {
			/*
			 * read in a list of all the vcfs in the individual folder
			 * 
			 */
			
			
		}
	}
	
	/**
	 * downloads vcfs from cohort and case mode
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void downloadSegmentsVCFs(String entityPath, String membershipPath, String wd) throws IOException, InterruptedException {
		DataFrame sampleSetMembership = new DataFrame(membershipPath, true, "\\t", "@");	
		DataFrame sampleSetEntity = new DataFrame(entityPath, true, "\\t", "@");
		HashMap<String, String> all_bed_paths = new HashMap<>();
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
	
	
	
	public void match() {
		HashMap<String, String> qs_map = new HashMap<String, String>();
		HashMap<String, String> np_map = new HashMap<String, String>();
		HashMap<String, String> cn_map = new HashMap<String, String>();
		
		for(int i = 0; i < svtk_input.size(); i++) {
			qs_map.put(svtk_input.get("name", i), svtk_input.get("qs", i));
			np_map.put(svtk_input.get("name", i), svtk_input.get("np", i));
			cn_map.put(svtk_input.get("name", i), svtk_input.get("cn", i));
		}
		
		
		
		for(int i = 0; i < svtk_output.size(); i++) {
			
		}
		
		
	}
	
	
	public void writeFile() {
		
	}
	
	
	

}
