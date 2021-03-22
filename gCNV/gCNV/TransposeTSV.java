package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class TransposeTSV extends gCNVHelperTool {

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH};
	
	public TransposeTSV(ArgParser args) {
		super(args, inputs);
	}
	
	@Override
	public void run() throws IOException, InterruptedException {
		String INPUT = args.get(INPUT_PATH);
		String OUTPUT = args.get(OUTPUT_PATH);
		this.transposeTSV(INPUT, OUTPUT);
	}
	

	public void transposeTSV(String INPUT, String OUTPUT) throws IOException {
		ArrayList<String[]> rows = new ArrayList<>();


		FileInputStream inputStream = new FileInputStream(INPUT);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		print("reading:");
		while (sc.hasNextLine()) {
			rows.add(sc.nextLine().split("\\t"));
		}
		sc.close();
		print("done reading");

		print("writing");

		BufferedWriter output = null;
		File file = new File(OUTPUT);
		output = new BufferedWriter(new FileWriter(file));

		int nrow = rows.get(0).length; 
		int ncol = rows.size();
		print("nrow: " + nrow);
		print("ncol: " + ncol);
		for(int r = 0; r < nrow; r++) {
			StringBuilder line = new StringBuilder();
			for(int c = 0; c < ncol-1; c++) {
				line.append(rows.get(c)[r]);
				line.append("\t");
			}
			line.append(rows.get(ncol-1)[r]);
			line.append("\n");
			output.write(line.toString());
		}
		output.close();


	}


}
