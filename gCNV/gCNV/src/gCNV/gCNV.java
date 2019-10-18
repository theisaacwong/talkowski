package gCNV;

/**
 * Author: Isaac Wong
 * October 9, 2019
 * Massachusetts General Hospital
 * Center for Genomic Medicine
 * Talkowski lab
 * 
 * This takes as input the merged bed file from gcnv and the bed file from svtk bedcluster. 
 * It then matches the qs/np/cn scores from the gcnv output to the svtk output
 * this is much faster than doing it in R
 * 
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class gCNV {

	public svtk_bed_input input;
	public svtk_bed_output output;
	public String input_string;
	public String output_string;
	public String OUTPUT_PATH;
	
	public gCNV() {
		input = null;
		output = null;
		input_string = null;
		output_string = null;
		OUTPUT_PATH = null;
	}
	
	public static void main(String[] args) throws IOException {
		
//		args = new String[]{"C:/Users/iwong/Documents/MGH/CMG/sample_set_10_15_19/full_merged_bed.bed", 
//				"C:/Users/iwong/Documents/MGH/CMG/sample_set_10_15_19/full_merged_bed_svtk_output.bed",
//				"C:/Users/iwong/Documents/MGH/CMG/sample_set_10_15_19/full_merged_svtk_java_out_temp.bed"};
		
		if(args.length == 1) {
			System.out.println("usage: gCNV [svtk_input] [svtk_output] [output_file]");
		} else {
			gCNV gcnv = new gCNV();
			gcnv.input_string = args[0];
			gcnv.output_string = args[1];
			gcnv.OUTPUT_PATH = args[2];
			gcnv.svtk_output_merge();	
		}
		
		
	}
	
	
	public void svtk_output_merge() throws IOException {
		input = new svtk_bed_input(input_string);
		output = new svtk_bed_output(output_string);
		
		System.out.println(input.size());
		System.out.println(output.size());
		
		HashMap<String, Integer> m_qs = new HashMap<String, Integer>();
		HashMap<String, Integer> m_np = new HashMap<String, Integer>();
		HashMap<String, Integer> m_cn = new HashMap<String, Integer>();
		
		for(svtk_bed_input_row r : input.df) {
			m_qs.put(r.name, r.qs);
			m_np.put(r.name, r.np);
			m_cn.put(r.name, r.copy_num);
		}

		svtk_merged rval = new svtk_merged(output);
		
		System.out.println(rval.df.size());
		
		for(svtk_merged_row r : rval.df) {
			//System.out.println(r.call_name);
			r.qs = m_qs.get(r.call_name);
			r.np = m_np.get(r.call_name);
			r.cn = m_cn.get(r.call_name);
		}	
		
		BufferedWriter output = null;
		File file = new File(OUTPUT_PATH);
		output = new BufferedWriter(new FileWriter(file));
		String header = "chr\tstart\tend\tname\tsvtype\tsample\tcall_name\tvaf\tvac\tpre_rmsstd\tpost_rmsstd\tnp\tcn\tqs\n";
		output.write(header);
		for(svtk_merged_row r : rval.df) {
			output.write(r.toString() + "\n");
		}
		output.close();
		
		
		
		
	}
	
	
	
	
	
}
