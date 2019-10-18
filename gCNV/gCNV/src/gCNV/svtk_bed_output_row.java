package gCNV;

public class svtk_bed_output_row {
	public String chr;		// 1
	public int start;		// 2
	public int end;			// 3
	public String name;		// 4
	public String svtype;	// 5
	public String sample;	// 6
	public String call_name;	// 7
	public String vaf;			// 8 I might want ot make these be strings
	public String vac;			// 9
	public String pre_rmsstd;	// 10
	public String post_rmsstd;	/// 11
	
	public svtk_bed_output_row() {
		
	}
	
	public svtk_bed_output_row(String chr, int start, int end, String name, String svtype, String sample, String call_name, String vaf, String vac, String pre_rmsstd, String post_rmsstd) {
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.name = name;
		this.svtype = svtype;
		this.sample = sample;
		this.call_name = call_name;
		this.vaf = vaf;
		this.vac = vac;
		this.pre_rmsstd = pre_rmsstd;
		this.post_rmsstd = post_rmsstd;
	}
	
	
	
}
