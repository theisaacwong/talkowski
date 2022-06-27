package gCNV;

import java.util.Arrays;
import java.util.HashMap;

public class ArgParser {

	public HashMap<String, String> parsedArgs;
	public String[] args;
	public boolean isValid;
	
	public ArgParser(String[] args) {
		System.out.println("The code to parse arguements just went through a major rework, please report all issues to Isaac");
		
		this.parsedArgs = new HashMap<>();
		this.args = args;
		this.isValid = true;
		this.parse();
	}
	
	public HashMap<String, String> parse(){
		for(int i = 0; i < this.args.length; i++) {
			if(args[i].startsWith("--")) {
				parsedArgs.put(args[i], "<empty>");
			} else if(args[i].startsWith("-")) {
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
		
		if(parsedArgs.containsKey("toolname") == false) {
			parsedArgs.put("toolname", "help");
			parsedArgs.put("--help", "help");
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
		
		for(String k : requiredArgs) {
			System.out.println("required: " + k);
		}
		for(String k : parsedArgs.keySet()) {
			System.out.println("parsed: " + k + " : " + parsedArgs.get(k));
		}
		
		if(parsedArgs.keySet().containsAll(Arrays.asList(requiredArgs))) {
			return true;
		}
		
		for(String requiredArg : requiredArgs) {
			if(parsedArgs.keySet().contains(requiredArg) == false)
				System.out.println("missing: " + requiredArg);
		}
		this.isValid = false;
		
		return false;
	}
		
		
		
}
