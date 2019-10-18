package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class DataFrame {
	public HashMap<String, Integer> columnMapping;
	public ArrayList<String[]> df;
	public String[] fieldNames;
	public String delimiter;
	public String FILE_PATH;
	public char comment_char;
	public boolean header;
	
	public DataFrame() throws IOException {
		columnMapping = new HashMap<String, Integer>();
		df = new ArrayList<String[]>();
		fieldNames = new String[1];
		
	}
	
	
	/*
	 * for reading in files
	 */
	public DataFrame(String file_input_path, boolean header, String sep, char comment_char) throws IOException {
		this.columnMapping = new HashMap<String, Integer>();
		this.df = new ArrayList<String[]>();
		this.delimiter = sep;
		this.comment_char = comment_char;
		
		this.parseFile();
	}
	public DataFrame(String[] colnames) {
		this.columnMapping = new HashMap<String, Integer>();
		this.df = new ArrayList<String[]>();
	}
	
	
	
	public String[] get(int i) {
		return df.get(i);
	}
	
	public String get(int i, int k) {
		return df.get(i)[k];
	}
	
	public int size() {
		return df.size();
	}
	
	public void add(String[] row) {
		this.df.add(row);
	}
	
	public void writeFile(String OUTPUT_PATH) throws IOException {
		BufferedWriter output = null;
		File file = new File(OUTPUT_PATH);
		output = new BufferedWriter(new FileWriter(file));
		String header = "";
		output.write(header);
		String line = "";
		for(String[] row : df) {
			line = "";
			for(String s : row) {
				line += s + delimiter;
			}
			output.write(line + "\n");
		}
		output.close();
	}
	
	public void parseFile() throws IOException {
		FileInputStream inputStream = new FileInputStream(FILE_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		if(header == true) {
			line = sc.nextLine();
			String[] linee = line.split(delimiter);
			this.fieldNames = linee;
			for(int i = 0; i < linee.length; i++) {
				columnMapping.put(linee[i], i);
			}
		}
		
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.charAt(0) == comment_char) {continue;}
			df.add(line.split(delimiter));
		}
		inputStream.close();
		sc.close();		
	}

	public String get(String columnName, int indexOfValue)
	{
		return df.get(columnMapping.get(columnName))[indexOfValue];
	}

}
