package gCNV;

public class svtk_bed_input_row {
	public String chr;		// 1
	public int start;		// 2
	public int end;			// 3
	public String name;		// 4
	public String sample;	// 5
	public String cn;		// 6
	public int qs;			// 7
	public int size;		// 8
	public int copy_num;	// 9
	public int np;			// 10
	public String group;	/// 11
	
	public svtk_bed_input_row() {
		
	}
	
	public svtk_bed_input_row(String chr, 
			int start, 
			int end, 
			String name, 
			String sample, 
			String cn, 
			int qs, 
			int size, 
			int copy_num, 
			int np, 
			String group) {
		
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.name = name;
		this.sample = sample;
		this.cn = cn;
		this.qs = qs;
		this.size = size;
		this.copy_num = copy_num;
		this.np = np;
		this.group = group;
	}
	
	
	
}
