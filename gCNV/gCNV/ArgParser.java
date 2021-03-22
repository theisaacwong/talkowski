package gCNV;

import java.util.HashMap;

public class ArgParser {

	public HashMap<String, String> parsedArgs;
	public String[] args;
	
	public ArgParser(String[] args) {
		this.parsedArgs = new HashMap<>();
		this.args = args;
		
	}
	
	public HashMap<String, String> parse(){
		for(int i = 0; i < this.parsedArgs.size(); i++) {
			if(args[i].contains("--")) {
				parsedArgs.put(args[i], "<empty>");
			} else if(args[i].contains("-")) {
				parsedArgs.put(args[i], args[i+1]);
			} else {
				parsedArgs.put("toolname", args[i]);
			}
		}
		
		return this.parsedArgs;
	}
	
	public boolean contains(String str) {
		for(String arg : parsedArgs.keySet()) {
			if(arg.equals(str)) {
				return(true);
			}
		}
		return(false);
	}
	
	public String get(String key) {
		return(parsedArgs.get(key));
	}
		
		
		
}
