package gCNV;

import java.util.ArrayList;

public class Gene implements Comparable<Gene>{
	public final String name;
	public final String chr;
	public final ArrayList<Integer> starts;
	public final ArrayList<Integer> ends;
	public final int minStart;
	public final int maxEnd;
	public final String ID;
	public final int nexons;
	public final int nbases;
	
	public Gene(String name, String chr, ArrayList<Integer> rstarts, ArrayList<Integer> rends) {
		this.name = name;
		this.chr = chr;
		
		//condense
		this.starts = new ArrayList<>();
		this.ends = new ArrayList<>();
		int i = 0;
		while(i < rstarts.size()) {
			int currentEnd = rends.get(i);
			int j = i+1;
			while(j < rstarts.size() && currentEnd > rstarts.get(j)) {
				currentEnd = Math.max(currentEnd, rends.get(j));
				j++;
			}
			this.starts.add(rstarts.get(i));
			this.ends.add(currentEnd);
			i=j;
		}
		
		this.nexons = this.starts.size();
		this.minStart = min(this.starts);
		this.maxEnd = max(this.ends);
		this.ID = chr + ":" + minStart + "-" + maxEnd + "_" + name;
		
		int sum = 0;
		for(int k = 0; k < this.nexons; k++) {
			sum += ends.get(k) - starts.get(k);
		}
		this.nbases = sum;
		
		if(starts.size() != ends.size()) {
			System.out.println("error");
		}
	}
	
	/**
	 * https://stackoverflow.com/questions/16691524/calculating-the-overlap-distance-of-two-1d-line-segments
	 * @param min1
	 * @param max1
	 * @param min2
	 * @param max2
	 * @return
	 */
	public int getNOverlap(int min1, int max1, int min2, int max2) {
		    return Math.max(0, Math.min(max1, max2) - Math.max(min1, min2));

	}
	
	/**
	 * calculate percentage overlap of a CNV, assuming chr is equal
	 * @param cnvStart
	 * @param cnvEnd
	 * @return
	 */
	public double calculateOverlapPercent(int cnvStart, int cnvEnd) {
		int cumulativeOverlap = 0;
		for(int i = 0; i < nexons; i++) {
			cumulativeOverlap += getNOverlap(cnvStart, cnvEnd, starts.get(i), ends.get(i)); 
		}
		return (double)cumulativeOverlap/ (double)nbases;
	}
	
	public void condense() {
		ArrayList<Integer> newStarts = new ArrayList<>();
		ArrayList<Integer> newEnds = new ArrayList<>();
		
		int i = 0;
		while(i < starts.size()) {
			int currentEnd = ends.get(i);
			int j = i+1;
			while(currentEnd > starts.get(j)) {
				currentEnd = Math.max(currentEnd, ends.get(j));
				j++;
			}
			newStarts.add(starts.get(i));
			newEnds.add(currentEnd);
			i=j;
		}
	}
	
	public String toString() {
		return this.ID;
	}
	
	public static int min(ArrayList<Integer> arr) {
		int min = arr.get(0);
		for(int i : arr) {
			min = Math.min(min, i);
		}
		return min;
	}
	
	public static int max(ArrayList<Integer> arr) {
		int max = arr.get(0);
		for(int i : arr) {
			max = Math.max(max, i);
		}
		return max;
	}
	
	@Override
	public int compareTo(Gene otherGene) {
		if(this.name.equals(otherGene.name)) {
			return 0;
		} else if(this.chr.equals(otherGene.chr)) {
			if(this.starts.get(0) > otherGene.starts.get(0)) {
				return 1;
			} else if(this.starts.get(0) < otherGene.starts.get(0)) {
				return -1;
			} else {
				return 0;
			}
		} else {
			return this.chr.compareTo(otherGene.chr);
		}
	}
		
	
	
}