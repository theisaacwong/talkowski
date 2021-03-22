package gCNV;

import java.io.IOException;
import java.io.InputStream;
import java.util.Scanner;

public class Help extends gCNVHelperTool {

	public static String[] inputs = new String[] {};
	
	public Help(ArgParser args) {
		super(args, inputs);
	}

	@Override
	public void run() throws IOException, InterruptedException {
		printOptionsShort();
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
	
}
