package gCNV;

import java.util.Arrays;
import java.util.HashMap;

public class ArgParser {

	public HashMap<String, String> parsedArgs;
	public String[] args;
	
	public ArgParser(String[] args) {
		this.parsedArgs = new HashMap<>();
		this.args = args;
		this.parse();
	}
	
	public HashMap<String, String> parse(){
		for(int i = 0; i < this.parsedArgs.size(); i++) {
			if(args[i].contains("--")) {
				parsedArgs.put(args[i], "<empty>");
			} else if(args[i].contains("-")) {
				parsedArgs.put(args[i], args[i+1]);
				i++;
			} else if(parsedArgs.containsKey("toolname") == false){
				parsedArgs.put("toolname", args[i]);
			} else {
				if(parsedArgs.containsKey("extra") == false) {
					parsedArgs.put("extra", args[i]);
				} else {
					parsedArgs.put("extra", parsedArgs.get("extra") + "," + args[i]);
				}
				
			}
		}
		
		return this.parsedArgs;
	}
	
	public int size() {
		return parsedArgs.size();
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
	
	public boolean validate(String [] requiredArgs) {
		
		if(parsedArgs.keySet().containsAll(Arrays.asList(requiredArgs))) {
			return true;
		}
		
		return false;
	}
		
		
		
}
