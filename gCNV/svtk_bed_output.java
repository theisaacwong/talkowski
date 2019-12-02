package gCNV;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class svtk_bed_output {
	public ArrayList<svtk_bed_output_row> df;
	public String INPUT_PATH;
	
	public svtk_bed_output() {
		this.df = new ArrayList<svtk_bed_output_row>();
	}
	
	public svtk_bed_output(String str) throws IOException {
		this.df = new ArrayList<svtk_bed_output_row>();
		this.INPUT_PATH = str;
		readInput();
	}
	
	public String getVal(String columnName, int rowNumber) {
		return "na";
	}
	public int size() {
		return this.df.size();
	}
	
	public void readInput() throws IOException {
		FileInputStream inputStream = new FileInputStream(INPUT_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		long lineNumber = 1;
		
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.charAt(0) == '@' || line.charAt(0) == '#') {continue;}
			String[] linee = line.split("\t");
			if(linee.length != 11) {System.out.println("Error at line " + lineNumber + ", expected 11 columns, got " + linee.length); continue;}
			this.df.add(new svtk_bed_output_row(
					linee[0],
					Integer.parseInt(linee[1]),
					Integer.parseInt(linee[2]),
					linee[3],
					linee[4],
					linee[5],
					linee[6],
					linee[7],
					linee[8],
					linee[9],
					linee[10]
					));
			lineNumber++;
			
		}
		inputStream.close();
		sc.close();
	}

}
