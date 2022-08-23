package gCNV;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.Calendar;
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
	public static final String VERSION = "2.30";

	public gCNV_helper(String[] args) {
		initializationArgs = args;
		date = Calendar.getInstance().getTime().toString();
	}

	public static void main(String[] args) throws IOException, InterruptedException {

		gCNV_helper g = new gCNV_helper(args);
		System.out.println(g.toString());
		System.out.println("version " + VERSION);
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
	
	public void run(String[] args) throws IOException, InterruptedException {
		ArgParser userArgs = new ArgParser(args);
		gCNVHelperTool gCNVHelperTool = null;
		
//		if(userArgs.contains("--no-check") == false) {
//			this.checkVersion();
//		}
		
		System.out.println();
		if (args.length == 0) {
			return;
		}
		
		String toolName = userArgs.get("toolname");
		switch (toolName) {
		case "help" -> {
			gCNVHelperTool = new Help(userArgs);
		}
		case "addGnomadAnnotations" -> {
			gCNVHelperTool = new AddGnomadAnnotations(userArgs);
		}
		case "annotateWithGenes" -> {
			gCNVHelperTool = new AnnotateWithGenes(userArgs);
		}
		case "bedcluster" -> {
			gCNVHelperTool = new Bedcluster(userArgs);
		}
		case "bedToVCF" -> {
			gCNVHelperTool = new BedToVCF(userArgs);
		}
		case "calculateFrequency" -> {
			gCNVHelperTool = new CalculateFrequency(userArgs);
		}
		case "condenseBedtoolsIntersect" -> {
			gCNVHelperTool = new CondenseBedtoolsIntersect(userArgs);
		}
		case "convertCoordinatesToVariantMedianValues" -> {
			gCNVHelperTool = new ConvertCoordinatesToVariantMedianValues(userArgs);
		}
		case "convertToEnsemble" -> {
			gCNVHelperTool = new ConvertToEnsemble(userArgs);
		}
		case "convertVCFsToBEDFormat" -> {
			gCNVHelperTool = new ConvertVCFsToBEDFormat(userArgs);
		}
		case "countExons" -> {
			gCNVHelperTool = new CountExons(userArgs);
		}
		case "defragment" -> {
			gCNVHelperTool = new Defragment(userArgs);
		}
		case "downloadFilteredIntervals" -> {
			gCNVHelperTool = new DownloadFilteredIntervals(userArgs);
		}
		case "downloadSegmentsVCFs" -> {
			gCNVHelperTool = new DownloadSegmentsVCFs(userArgs);
		}
		case "filter" -> {
			gCNVHelperTool = new Filter(userArgs);
		}
		case "getBarcodeCounts" -> {
			gCNVHelperTool = new GetBarcodeCounts(userArgs);
		}
		case "getBEDtrack" -> {
			gCNVHelperTool = new GetBEDtrack(userArgs);
		}
		case "getCountsMatrix" -> {
			gCNVHelperTool = new GetCountsMatrix(userArgs);
		}
		case "getCountsMatrixBuffered" -> {
			gCNVHelperTool = new GetCountsMatrixBuffered(userArgs);
		}
		case "getPerSampleMetrics" -> {
			gCNVHelperTool = new GetPerSampleMetrics(userArgs);
		}
		case "jointCallVCF" -> {
			gCNVHelperTool = new JointCallVCF(userArgs);
		}
		case "labelMedianQS" -> {
			gCNVHelperTool = new LabelMedianQS(userArgs);
		}
		case "simplify" -> {
			gCNVHelperTool = new Simplify(userArgs);
		}
		case "subsetAnnotations" -> {
			gCNVHelperTool = new SubsetAnnotations(userArgs);
		}
		case "transposeTSV" -> {
			gCNVHelperTool = new TransposeTSV(userArgs);
		}
		case "updateNP" -> {
			gCNVHelperTool = new updateNP(userArgs);
		}
		case "validateSubsetAnnotations" -> {
			gCNVHelperTool = new ValidateSubsetAnnotations(userArgs);
		}
		default -> {
			System.out.println("unknown command: " + toolName);
//			gCNVHelperTool = new Help(userArgs);
		}

		}

		if(gCNVHelperTool.validate()) {
			gCNVHelperTool.run();	
		}
		

	}



}
