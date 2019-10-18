package gCNV;

//import java.io.BufferedWriter;
//import java.io.File;
//import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class gCNV_2_00 {
	
	
	public DataFrame svtk_input;
	public DataFrame svtk_output;
	public DataFrame gCNV_output;
	public String svtk_input_path;
	public String svtk_output_path;
	public String output_file_path;
	
	public gCNV_2_00(String svtk_in, String svtk_out, String out_file_path) throws IOException {
		svtk_input = new DataFrame(svtk_in, true, "\t", '#');
		svtk_output = new DataFrame(svtk_out, true, "\t", '#');
		gCNV_output = new DataFrame();
		
	}
	
	
	public static void main_not(String[] args) throws IOException {
		
		gCNV_2_00 gCNV = new gCNV_2_00(args[0], args[1], args[2]);
		gCNV.match();
		gCNV.writeFile();
		
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
