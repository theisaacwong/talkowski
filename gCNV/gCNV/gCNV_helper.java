package gCNV;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.Calendar;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
	public static final String VERSION = "2.25";

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
		start();
		g.run(args);
		stop();
		
	}


	public String toString() {
		return String.join(" ", initializationArgs) + "\n" + date;
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

	public static String getTimeHumanReadable(long time) {
		long seconds = time % 60;
		long minutes = (time / 60) % 60;
		long hours = (time / 3600);

		return "Hours: " + hours + "  Min: " + minutes + "  Sec: " + seconds;
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

	public void printOptionsShort() {
		System.out.println("java -jar gCNV_helper.jar [Command] [required argument(s)] ");
		System.out.println();
		InputStream in = this.getClass().getResourceAsStream("/helpMenuShort"); 
		Scanner sc = new Scanner(in, "UTF-8");
		String line = "";
		while(sc.hasNextLine()) {
			line = sc.nextLine();
			System.out.println(line);
		}
		sc.close();
	}

	public void run0(String[] args) throws IOException, InterruptedException {
		
		gCNVHelperTool gCNVHelperTool;
		System.out.println();
		if (args.length == 0) {
			printOptionsShort();
			return;
		}
		String toolName = args[0];
		switch (toolName) {
		case "addGnomadAnnotations" -> {
			gCNVHelperTool = new AddGnomadAnnotations(args);
		}
		default -> {
			gCNVHelperTool = new gCNVHelperTool(args);
		}
		}
		
		gCNVHelperTool.run();
		
	}
	
	public void run(String[] args) throws IOException, InterruptedException {
		gCNVHelperTool gCNVHelperTool = null;
		
		System.out.println();
		if (args.length == 0) {
			printOptionsShort();
			return;
		}
		String toolName = args[0];
		switch (toolName) {
		case "-help" -> {
			printOptions();
		}
		case "--help" -> {
			printOptions();
		}
		case "-h" -> {
			printOptions();
		}
		case "addGnomadAnnotations" -> {
			gCNVHelperTool = new AddGnomadAnnotations(args);
		}
		case "annotateWithGenes" -> {
			gCNVHelperTool = new AnnotateWithGenes(args);
		}
		case "bedcluster" -> {
			gCNVHelperTool = new Bedcluster(args);
		}
		case "calculateFrequency" -> {
			gCNVHelperTool = new CalculateFrequency(args);
		}
		case "condenseBedtoolsIntersect" -> {
			gCNVHelperTool = new CondenseBedtoolsIntersect(args);
		}
		case "convertCoordinatesToVariantMedianValues" -> {
			gCNVHelperTool = new ConvertCoordinatesToVariantMedianValues(args);
		}
		case "convertToEnsemble" -> {
			gCNVHelperTool = new ConvertToEnsemble(args);
		}
		case "convertVCFsToBEDFormat" -> {
			gCNVHelperTool = new ConvertVCFsToBEDFormat(args);
		}
		case "countExons" -> {
			gCNVHelperTool = new CountExons(args);
		}
		case "defragment" -> {
			gCNVHelperTool = new Defragment(args);
		}
		case "downloadSegmentsVCFs" -> {
			gCNVHelperTool = new DownloadSegmentsVCFs(args);
		}
		case "filter" -> {
			gCNVHelperTool = new Filter(args);
		}
		case "getBarcodeCounts" -> {
			gCNVHelperTool = new GetBarcodeCounts(args);
		}
		case "getBEDtrack" -> {
			gCNVHelperTool = new GetBEDtrack(args);
		}
		case "getCountsMatrix" -> {
			gCNVHelperTool = new GetCountsMatrix(args);
		}
		case "getCountsMatrixBuffered" -> {
			gCNVHelperTool = new GetCountsMatrixBuffered(args);
		}
		case "getPerSampleMetrics" -> {
			gCNVHelperTool = new GetPerSampleMetrics(args);
		}
		case "jointCallVCF" -> {
			gCNVHelperTool = new JointCallVCF(args);
		}
		case "labelMedianQS" -> {
			gCNVHelperTool = new LabelMedianQS(args);
		}
		case "simplify" -> {
			gCNVHelperTool = new Simplify(args);
		}
		case "subsetAnnotations" -> {
			gCNVHelperTool = new SubsetAnnotations(args);
		}
		case "transposeTSV" -> {
			gCNVHelperTool = new TransposeTSV(args);
		}
		case "validateSubsetAnnotations" -> {
			gCNVHelperTool = new ValidateSubsetAnnotations(args);
		}
		default -> {
			System.out.println("unknown command: " + toolName);
			gCNVHelperTool = new gCNVHelperTool(args);
		}

		}

		gCNVHelperTool.run();

	}



}
