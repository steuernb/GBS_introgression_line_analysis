

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;
import java.util.zip.GZIPInputStream;


public class GBSSNPs {
	Hashtable<String,HashSet<Long>> coverage;
	Hashtable<String,HashSet<Long>> unusable;
	Hashtable<String,HashSet<Long>> usable;
	
	
	public GBSSNPs(File pileupBackground, File pileupDonor)throws IOException{
		
		 this.coverage = new Hashtable<String,HashSet<Long>>();
		 this.unusable = new Hashtable<String,HashSet<Long>>();
		 this.usable = new Hashtable<String,HashSet<Long>>();
		 
		 this.loadBackGround(pileupBackground);
		 this.loadDonor(pileupDonor); 
		 
		 
		
		 
	}
	
	
	
	public static void main(String[] args) {
		CLI cli = new CLI();
		cli.parseOptions(args);
		
		
		try {
			
			if( !cli.hasOption("b")) {
				throw new CLIParseException("No background file specified. Parameter -b missing.");
			}
			if( !cli.hasOption("d")) {
				throw new CLIParseException("No donor file specified. Parameter -d missing.");
			}
			if( !cli.hasOption("i")) {
				throw new CLIParseException("No introgression file specified. Parameter -i missing.");
			}
			if( !cli.hasOption("o")) {
				throw new CLIParseException("No output file specified. Parameter -o missing.");
			}
			if( !cli.hasOption("p")) {
				throw new CLIParseException("No bed file parts2chromosomes specified. Parameter -p missing.");
			}
			if( !cli.hasOption("l")) {
				throw new CLIParseException("No length of interval specified. Parameter -l missing.");
			}
			
			
			File backgroundPileup = new File(cli.getArg("b"));
			File donorPileup = new File(cli.getArg("d"));
			File introgressionPileup = new File(cli.getArg("i"));
			File outputFile = new File(cli.getArg("o"));
			File parts2chromosomes = new File(cli.getArg("p"));
			int intervalLength = Integer.parseInt(cli.getArg("l"));
			
			
			GBSSNPs gbs = new GBSSNPs(backgroundPileup, donorPileup);
			gbs.checkSNPs(parts2chromosomes, introgressionPileup, outputFile, intervalLength);
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		
	}
	
	
	
	
	
	
	public   void checkSNPs(File cs_parts2chromosomes, File pileupIntrogression,  File outputFile, int interval)throws IOException{
		
		
		
		
		
		Hashtable<String,Long> offset = new Hashtable<String,Long>();
		BufferedReader in = new BufferedReader(new FileReader(cs_parts2chromosomes));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			offset.put(split[0], Long.parseLong(split[4]));
		}

		in.close();
		
	
		Hashtable<String, HashSet<Long>> mySNPs = new Hashtable<String, HashSet<Long>>();
		
		BufferedReader in4;
		FileInputStream fis = new FileInputStream(pileupIntrogression);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();

		
		if(gzip){
			in4=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pileupIntrogression))));
		}else{
			in4 = new BufferedReader(new FileReader(pileupIntrogression));
		}
		
		
		 
		
		
		
		for (String inputline = in4.readLine(); inputline!= null; inputline = in4.readLine()){
			MPileupLine line = new MPileupLine(inputline);
			if( line.getReferenceAlleleFrequency() <0.1){
				long pos =  line.getPosition() ;
				String chr = line.getChromosome();
				
				if( (coverage.containsKey(chr) && coverage.get(chr).contains(pos)  )&&  ( !unusable.containsKey(chr) || !unusable.get(chr).contains(pos)) && (usable.containsKey(chr) &&usable.get(chr).contains(pos))){
					
				if(!mySNPs.containsKey(chr)){
						mySNPs.put(chr, new HashSet<Long>());
					}
					mySNPs.get(chr).add(pos);
				}
				
			}
			
		}
		
		in4.close();
		
		int countPositions = 0;
		for(Enumeration<String> myenum = mySNPs.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			countPositions = countPositions + mySNPs.get(key).size();
		}
		System.out.println(pileupIntrogression + " Positions found: " + countPositions);
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Enumeration<String> myenum1 = usable.keys(); myenum1.hasMoreElements();){
			String chr = myenum1.nextElement();
			Vector<Long> v = new Vector<Long>();
			for(Iterator<Long> iterator = usable.get(chr).iterator(); iterator.hasNext();){
				v.add(iterator.next());
			}
			Collections.sort(v);
			
			int countDonor =0;
			int count = 0;   //count SNPs in introgression line that are the same as in donor
			int base = 0;
			
			
			for(Enumeration<Long> myenum2 = v.elements(); myenum2.hasMoreElements();){
				long pos = myenum2.nextElement();
				
				if( pos > base + interval){
					out.write( chr.split("_")[0] + "\t" + (base + interval/2 + offset.get(chr).longValue()) +"\t" + countDonor +"\t" + count  );
					out.newLine();
					base = base + interval;
					count = 0;
					countDonor = 0;
				}
				
				countDonor ++;
				if(mySNPs.containsKey(chr) && mySNPs.get(chr).contains(pos)){
					count ++;
				}
			}
			
			out.write( chr.split("_")[0] + "\t" + (base + interval/2 +offset.get(chr).longValue()) +"\t" + countDonor +"\t" + count  );
			out.newLine();
		}
		
		
		out.close();
		
		
		
		
		
	}
	
	
	
	
	private void loadBackGround(File pileupBackground)throws IOException{
		
			
		BufferedReader in2;
		FileInputStream fis = new FileInputStream(pileupBackground);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();

		
		if(gzip){
			in2=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pileupBackground))));
		}else{
			in2 = new BufferedReader(new FileReader(pileupBackground));
		}
		
			
			
		

		for(String inputline = in2.readLine(); inputline!= null; inputline = in2.readLine()){
			if( inputline.startsWith("#")){
				continue;
			}
			
			MPileupLine line = new MPileupLine(inputline);
			if( !coverage.containsKey(line.getChromosome())){
				coverage.put(line.getChromosome(), new HashSet<Long>());
				unusable.put(line.getChromosome(), new HashSet<Long>());
				System.out.print(".");
			}
			coverage.get(line.getChromosome()).add((long)line.getPosition());
			if( line.getReferenceAlleleFrequency() < 1){
				unusable.get(line.getChromosome()).add((long)line.getPosition());
				
			}
			
		}
		
		in2.close();
		
		
		
			
		
		
		long covered = 0;
		long unusabl = 0;
		for(Enumeration<String> myenum = coverage.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			covered = covered + coverage.get(key).size();
			unusabl = unusabl + unusable.get(key).size();
		}
		
		
		
		System.out.println("Positions covered: " + covered + "; positions unusable: " +  unusabl);
		
		
		
	}
	
	
	private void loadDonor(File pileupDonor)throws IOException{

		
		
		BufferedReader in3;
		FileInputStream fis = new FileInputStream(pileupDonor);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();

		
		if(gzip){
			in3=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pileupDonor))));
		}else{
			in3 = new BufferedReader(new FileReader(pileupDonor));
		}
		
		
		
		for(String inputline = in3.readLine(); inputline!= null; inputline = in3.readLine()){
			if( inputline.startsWith("#")){
				continue;
			}
			
			MPileupLine line = new MPileupLine(inputline);
			
			String chr = line.getChromosome();
			long pos = line.getPosition();
			
			
			if( !unusable.get(chr).contains(pos)){
				if( !usable.containsKey(line.getChromosome())){
					usable.put(line.getChromosome(), new HashSet<Long>());
					System.out.print(".");
				}
				if( line.getReferenceAlleleFrequency() <1){
					usable.get(chr).add(pos);
				}
			}
			
		}
		
		in3.close();
		
		
		long usabl = 0;
		
		for(Enumeration<String> myenum = usable.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			
			usabl = usabl + usable.get(key).size();
		}
		
		System.out.println("Positions usable: " +  usabl);
		
		
		
	}
	
	
	
	
	
	
	
	
	


	
	
	
	
}
