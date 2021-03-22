package gCNV;

import java.io.IOException;

public class GetPerSampleMetrics extends gCNVHelperTool {

	public GetPerSampleMetrics(String[] args) {
		super(args);
	}

	@Override
	public void run() throws IOException, InterruptedException {
		String INPUT_PATH = args[1];
		String fieldColumn = args[2];
		String OUTPUT_PATH = args[3];
		boolean writeFile = true;
		this.getPerSampleMetrics(INPUT_PATH, fieldColumn, OUTPUT_PATH, writeFile);
	}

	
}
