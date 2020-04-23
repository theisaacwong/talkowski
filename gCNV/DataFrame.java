package gCNV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;

import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;

/**
 * 
 * @author Isaac Wong
 * December 2, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski Lab
 * 
 * This class is designed as a bare bones version of R's data.frame
 * 
 * Contact Information
 * website: https://talkowski.mgh.harvard.edu/
 * email: iwong@broadinstitute.org
 * github: https://github.com/theisaacwong/
 *
 */

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
		delimiter = "\t";
		
	}
	
	/*
	 * for reading in files
	 */
	public DataFrame(String file_input_path, boolean header, String sep, String comment_char) throws IOException {
		if(System.getProperty("os.name").contains("indows")) {
			// do nothing, hopefully
		} else {
			file_input_path = file_input_path.replaceAll("\\", "/");
		}
		
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
	
	/**
	 * writes the current data frame to file
	 * @param OUTPUT_PATH
	 * @param header
	 * @throws IOException
	 */
	public void writeFile(String OUTPUT_PATH, boolean header) throws IOException {
		
		if(System.getProperty("os.name").contains("indows")) {
			// do nothing, hopefully
		} else {
			OUTPUT_PATH = OUTPUT_PATH.replaceAll("\\", "/");
		}
		
		BufferedWriter output = null;
		File file = new File(OUTPUT_PATH);
		output = new BufferedWriter(new FileWriter(file));
		
		StringBuilder sb = new StringBuilder();
		if(header == true) {
			for(int i = 0; i < this.fieldNames.length; i++) {
				sb.append(this.fieldNames[i]);
				if(i != this.fieldNames.length) {
					sb.append("\t");
				}
			} sb.append("\n");
			output.write(sb.toString());
			sb = new StringBuilder();
		}
		
		for(int i = 0; i < this.nrow(); i++) {
			output.write(String.join("\t", this.get(i)) + "\n");
		}
		System.out.println(this.nrow());
		output.close();
	}
	
	
    public Reader getReader(String relativePath) throws UnsupportedEncodingException {
    	try {
			return new InputStreamReader(this.getClass().getResourceAsStream(relativePath), "UTF-8");
		} catch (UnsupportedEncodingException e) {
			throw new IllegalStateException("Unable to read input", e);
		}
    }
    public Reader getFileReader(String absolutePath) throws UnsupportedEncodingException, FileNotFoundException {
        return new InputStreamReader(new FileInputStream(new File(absolutePath)), "UTF-8");
    }
	/**
	 * reads into memory the file at this object's path
	 * has errors with the large formats of fireloud's manifest file
	 * https://www.univocity.com/pages/univocity_parsers_tutorial.html#introduction-to-univocity-parsers
	 * @throws IOException
	 */
	public void parseFileU() throws IOException {
		TsvParserSettings settings = new TsvParserSettings();
		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setComment(comment_char.charAt(0));

		// creates a TSV parser
		TsvParser parser = new TsvParser(settings);

		// parses all rows in one go.
		List<String[]> allRows = parser.parseAll(getFileReader(FILE_PATH));
		
		if(header == true) {
			this.fieldNames = allRows.get(0).clone(); // works as intended
			for(int i = 0; i < allRows.get(0).length; i++) {
				columnMapping.put(allRows.get(0).clone()[i], i);
			}
			allRows.remove(0);
		}
		for(int i = 0; i < allRows.size(); i++) {
			df.add(allRows.get(i).clone());
		}
	}
	
	/**
	 * while not technically deprecated, you should use the newer version. There is nothing wrong with this version, but the newer version uses univocity to parse the tsv rather than my homemade tsv parser
	 * using external libraries is better because it offloads bug fixing on someone else
	 * this version still works perfectly fine against my test cases and is sill usable
	 * reads into memory the file at this object's path
	 * ***CORRECTION*** 
	 * use this method, univocity has issues with firecloud's nonstandard tsv peculularities
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
			df.add(line.split(delimiter, -1));
			for(int i = 0; i < df.get(df.size()-1).length; i++) {
				if(df.get(df.size()-1)[i].equals("")) {
					df.get(df.size()-1)[i] = "0";
				}
			}
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
	 * this is an expensive operation, it has to create a new array for each row
	 * @param columnName
	 * @param columnValues
	 * @return false if operation failed, true is succeeded
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
	 * @return false if operation failed, true if success
	 */
	public boolean addColumns(ArrayList<String> columnNames, ArrayList<ArrayList<String>> columnValues) {
		if(this.df.size() == 0) {
			//System.out.println("wefweewfwe");
			this.fieldNames = new String[columnNames.size()];
			for(int i = 0; i < this.fieldNames.length; i++) {
				this.fieldNames[i] = columnNames.get(i);
			}
			for(int i = 0; i < this.fieldNames.length; i++) {
				this.columnMapping.put(this.fieldNames[i], i);
			}

			int firstSize = columnValues.get(0).size();
			for(ArrayList<String> oneCol : columnValues) {
				if(oneCol.size() != firstSize) {
					return false;
				}
			}
			if(columnNames.size() != columnValues.size()) {
				return false;
			}
			
			for(int i = 0; i < columnValues.get(0).size(); i++) {
				this.df.add(new String[columnNames.size()]);
				for(int k = 0; k < columnNames.size(); k++) {
					this.df.get(i)[k] = columnValues.get(k).get(i);
				}
			}
			return true;
			
		}
		
		
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
		System.out.println(Arrays.toString(this.fieldNames));
		System.out.println(Arrays.toString(newFieldNames));
		this.fieldNames = newFieldNames;
		System.out.println(Arrays.toString(this.fieldNames));
		
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
	
	/**
	 * currently performs a shallow copy, deep copy utilities are not quite yet needed as only strings are currently supported
	 */
	public boolean rbind(DataFrame df) {
		//check to make sure columns are equal
		if(!this.columnMapping.equals(df.columnMapping)) {
			return false;
		}
		
		for(int i = 0; i < df.nrow(); i++) {
			this.df.add(df.get(i));
		}
		
		return true;
	}
	
	/**
	 * @throws IOException 
	 * 
	 */
	public void writeBed(String outputPath) throws IOException {
		if(System.getProperty("os.name").contains("indows")) {
			// do nothing, hopefully
		} else {
			outputPath = outputPath.replaceAll("\\", "/");
		}
		
		int chr = this.columnMapping.get("chr");
		int start = this.columnMapping.get("start");
		int end = this.columnMapping.get("end");
		int name = this.columnMapping.get("name");
		int sample = this.columnMapping.get("sample");
		int svtype = this.columnMapping.get("svtype");
		HashSet<Integer> remainginFields = new HashSet<>(this.columnMapping.values());
		remainginFields.remove(chr);
		remainginFields.remove(start);
		remainginFields.remove(end);
		remainginFields.remove(name);
		remainginFields.remove(sample);
		remainginFields.remove(svtype);
		ArrayList<Integer> remainingFields = new ArrayList<>(remainginFields);
		
		BufferedWriter output = null;
		File file = new File(outputPath);
		output = new BufferedWriter(new FileWriter(file));
		StringBuilder sb = new StringBuilder();
		
		//write header
		sb.append("#chr\t");
		sb.append("start\t");
		sb.append("end\t");
		sb.append("name\t");
		sb.append("sample\t");
		sb.append("svtype\t");
		for(int i = 0; i < remainingFields.size(); i++) {
			sb.append(fieldNames[remainingFields.get(i)]);
			if(i != remainingFields.size()) {
				sb.append("\t");
			}
		}
		sb.append("\n");
		output.write(sb.toString());
		sb = new StringBuilder();
		
		for(int i = 0; i < this.nrow(); i++) {
			sb = new StringBuilder();
			sb.append(this.get("chr", i)); sb.append("\t");
			sb.append(this.get("start", i)); sb.append("\t");
			sb.append(this.get("end", i)); sb.append("\t");
			sb.append(this.get("name", i)); sb.append("\t");
			sb.append(this.get("sample", i)); sb.append("\t");
			sb.append(this.get("svtype", i)); sb.append("\t");
			for(int k = 0; k < remainingFields.size(); k++) {
				sb.append(this.get(i, remainingFields.get(k)));
				
				if(k != remainingFields.size()) {
					sb.append("\t");
				}
			}
			sb.append("\n");
			output.write(sb.toString());
		}
		output.close();
		
		
	}
	

}
