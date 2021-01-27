package editRGSM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class EditRGSM {

	public static void main(String[] args) throws IOException {

		printInfo();
		
		editRGSM(args[0], args[1], args[2]);
		
	}

	public static void printInfo() {
		System.out.println("Java version: " + System.getProperty("java.version"));
		System.out.println("Heap Size: " + getHeapSize());
	}

	public static void editRGSM(String INPUT_PATH, String OUTPUT_PATH, String toChangeTo) throws IOException {

		Boolean hasBeenEdited = false;

		File file = new File(OUTPUT_PATH);
		BufferedWriter output = new BufferedWriter(new FileWriter(file));

		FileInputStream inputStream = new FileInputStream(INPUT_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		String line = "";

		while (!hasBeenEdited && sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.contains("@RG")) {
				output.write("@RG\tID:GATKCopyNumber\tSM:" + toChangeTo + "\n");
				hasBeenEdited = true;
			} else {
				output.write(line + "\n");
			}
		}

		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			output.write(line + "\n");
		}

		output.close();
		sc.close();
	}
	
	public static String getHeapSize() {
		return humanReadableByteCount(Runtime.getRuntime().totalMemory());
	}
	
	public static String humanReadableByteCount(long bytes) {
		if (bytes < 1024) return bytes + " B";
		int exp = (int) (Math.log(bytes) / 6.907755278982137);
		char pre = "KMGTPE".charAt(exp-1);
		return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
	}



}
