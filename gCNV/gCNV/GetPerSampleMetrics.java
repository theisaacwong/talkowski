package gCNV;

import java.io.IOException;

public class GetPerSampleMetrics extends gCNVHelperTool {

	public static String[] inputs = new String[] {INPUT_PATH, OUTPUT_PATH,COLUMN_NAME};

	public GetPerSampleMetrics(ArgParser args) {
		super(args, inputs);
	}

	@Override
	public void run() throws IOException, InterruptedException {
		String INPUT = args.get(INPUT_PATH);
		String fieldColumn = args.get(COLUMN_NAME);
		String OUTPUT = args.get(OUTPUT_PATH);
		boolean writeFile = true;
		this.getPerSampleMetrics(INPUT, fieldColumn, OUTPUT, writeFile);
	}

	
}
