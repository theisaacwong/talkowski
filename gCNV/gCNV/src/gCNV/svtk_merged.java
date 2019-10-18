package gCNV;

import java.util.ArrayList;

public class svtk_merged {

	public ArrayList<svtk_merged_row> df;

	public svtk_merged(ArrayList<svtk_merged_row> df) {
		this.df = df;
	}
	
	public svtk_merged(svtk_bed_output output) {
		this.df = new ArrayList<svtk_merged_row>();
		for(svtk_bed_output_row row : output.df) {
			if(row.call_name.contains(",")) {
				String[] splits = row.call_name.split(",");
				for(String s : splits) {
					this.df.add(new svtk_merged_row(
							row.chr,
							row.start,
							row.end,
							row.name,
							row.svtype,
							row.sample,
							s,
							row.vaf,
							row.vac,
							row.pre_rmsstd,
							row.post_rmsstd,
							-1,
							-1,
							-1
							));		
				}
			} else {
				this.df.add(new svtk_merged_row(
						row.chr,
						row.start,
						row.end,
						row.name,
						row.svtype,
						row.sample,
						row.call_name,
						row.vaf,
						row.vac,
						row.pre_rmsstd,
						row.post_rmsstd,
						-1,
						-1,
						-1
						));				
			}

		}
	}
	
	public int size() {
		return this.df.size();
	}
	
	

}
