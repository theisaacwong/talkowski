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
	public String comment_char;
	public boolean header;
	
	public DataFrame() throws IOException {
		columnMapping = new HashMap<String, Integer>();
		df = new ArrayList<String[]>();
		fieldNames = new String[1];
		
	}
	
	/*
	 * for reading in files
	 */
	public DataFrame(String file_input_path, boolean header, String sep, String comment_char) throws IOException {
		this.columnMapping = new HashMap<String, Integer>();
		this.df = new ArrayList<String[]>();
		this.delimiter = sep;
		this.comment_char = comment_char;
		this.FILE_PATH = file_input_path;
		this.header = header;
		this.parseFile();
	}
	public DataFrame(String[] colnames) {
		this.columnMapping = new HashMap<String, Integer>();
		this.df = new ArrayList<String[]>();
	}
	
	
	/**
	 * 
	 * @param i: row number
	 * @return String[] row
	 */
	public String[] get(int i) {
		return df.get(i);
	}
	
	/**
	 * 
	 * @param i: row number
	 * @param k: column number
	 * @return String element
	 */
	public String get(int i, int k) {
		return df.get(i)[k];
	}
	
	/**
	 * 
	 * @return nrow()
	 */
	public int size() {
		return df.size();
	}
	
	/**
	 * 
	 * @return number of rows in this data frame
	 */
	public int nrow() {
		return df.size();
	}
	
	/**
	 * 
	 * @return ncol()
	 */
	public int ncol() {
		return df.get(0).length;
	}
	
	/**
	 * adds a row to df
	 * @param row
	 */
	public void add(String[] row) {
		this.df.add(row);
	}
	
	public void writeFile(String OUTPUT_PATH, boolean header) throws IOException {
		BufferedWriter output = null;
		File file = new File(OUTPUT_PATH);
		output = new BufferedWriter(new FileWriter(file));

		StringBuilder sb = new StringBuilder();
		
		if(header==true) {
			for(int i = 0; i < fieldNames.length-1; i++) {
				sb.append(fieldNames[i]);
				sb.append(delimiter);
			}
			sb.append(fieldNames[fieldNames.length-1]);
			sb.append("\n");
			output.write(sb.toString());
			sb = new StringBuilder();
		}
		
		for(String[] row : df) {
			sb = new StringBuilder();
			for(String s : row) {
				sb.append(s);
				sb.append(delimiter);
			}
			sb.append("\n");
			output.write(sb.toString());
		}
		output.close();
	}
	
	/**
	 * reads into memory the file at this object's path
	 * https://stackoverflow.com/questions/1635764/string-parsing-in-java-with-delimiter-tab-t-using-split
	 * @throws IOException
	 */
	public void parseFile() throws IOException {
		FileInputStream inputStream = new FileInputStream(FILE_PATH);
		Scanner sc = new Scanner(inputStream, "UTF-8");

		String line = "";
		
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.indexOf(comment_char)==0) {
				continue;
			} else {
				break;
			}
		}
		
		if(header == true) {
//			line = sc.nextLine();
			String[] linee = line.split(delimiter);
			this.fieldNames = linee;
			for(int i = 0; i < linee.length; i++) {
				columnMapping.put(linee[i], i);
			}
		}
		
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();} catch(Exception e){System.out.println(e);}
			if(line.indexOf(comment_char) == 0) {continue;}
			df.add(line.split(delimiter, -1));
		}
		inputStream.close();
		sc.close();	
	}

	/**
	 * 
	 * @param columnName
	 * @param indexOfValue
	 * @return element at
	 */
	public String get(String columnName, int indexOfValue)
	{
		return df.get(indexOfValue)[columnMapping.get(columnName)];
	}
	
	/**
	 * 
	 * @param columnName
	 * @return returns an arraylist of that column
	 */
	public ArrayList<String> get(String columnName){
		ArrayList<String> rval = new ArrayList<>();
		int colNumber = columnMapping.get(columnName);
		for(String[] rowEntry : this.df) {
			rval.add(rowEntry[colNumber]);
		}
		return rval;
	}
	
	/**
	 * TODO: be able to add multiple columns at once
	 * this is an expensive operation, it has to create a new array for each row
	 * @param columnName
	 * @param columnValues
	 * @return
	 */
	public boolean addColumn(String columnName, ArrayList<String> columnValues) {
		if(columnValues.size() != this.df.size()) {
			return false;
		}
		if(this.columnMapping.containsKey(columnName)) {
			return false;
		}
		
		int nsize = this.fieldNames.length;
		this.columnMapping.put(columnName, nsize);
		String[] newFieldNames = new String[nsize+1];
		for(int i = 0; i < this.fieldNames.length; i++) {
			newFieldNames[i] = this.fieldNames[i];
		}
		newFieldNames[nsize] = columnName;
		this.fieldNames = newFieldNames;
		
		for(int i = 0; i < this.nrow(); i++) {
			String[] newRow = new String[nsize+1];
			for(int k = 0; k < this.get(i).length; k++) {
				newRow[k] = this.get(i,  k);
			}
			newRow[nsize] = columnValues.get(i);
			this.df.set(i, newRow);
		}
		
		return true;
	}
	
	/**
	 * adds multiple new columns at once
	 * @param columnNames, an arraylist of the new column names
	 * @param columnValues, an arraylist of column values, each of which is an array list
	 * @return
	 */
	public boolean addColumns(ArrayList<String> columnNames, ArrayList<ArrayList<String>> columnValues) {
		int firstSize = columnValues.get(0).size();
		for(ArrayList<String> oneCol : columnValues) {
			if(oneCol.size() != firstSize) {
				return false;
			}
		}
		
		if(columnNames.size() != columnValues.size()) {
			return false;
		}
		
		for(String colName : columnNames) {
			if(this.columnMapping.containsKey(colName)) {
				return false;
			}
		}
		
		int nsize = this.fieldNames.length;
		for(int i = 0; i < columnNames.size(); i++) {
			this.columnMapping.put(columnNames.get(i), nsize+i);
		}
		
		String[] newFieldNames = new String[nsize+columnNames.size()];
		for(int i = 0; i < this.fieldNames.length; i++) {
			newFieldNames[i] = this.fieldNames[i];
		}
		for(int i = 0; i < columnNames.size(); i++) {
			newFieldNames[i+nsize] = columnNames.get(i);
		}
		this.fieldNames = newFieldNames;
		
		for(int i = 0; i < this.nrow(); i++) {
			String[] newRow = new String[nsize+columnNames.size()];
			for(int k = 0; k < this.get(i).length; k++) {
				newRow[k] = this.get(i,  k);
			}
			for(int k = 0; k < columnValues.size(); k++) {
				newRow[nsize+k] = columnValues.get(k).get(i);
			}
			this.df.set(i, newRow);
		}
		
		return true;
	}
	

}
