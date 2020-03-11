/*
 * This program formats ActiveSequences.txt and InactiveSequences.txt and prints it to TrainingSequences.txt
 * It uses FeatureCalculator to calculate the peptide characteristics
 */
import java.io.*;
import java.util.*;


public class FilesReader 
{
	// read a file for HIV sequences, false sequences
	// later, classify sampled sequences
	public static BufferedReader activebr;
	public static BufferedReader uniprotbr; // reads data from uniprot
	public static PrintWriter inactiveout; // prints uniprot sequences to InactiveSequences
	public static BufferedReader inactivebr;
	public static PrintWriter inactiveout2; //prints selected 1000 sequences to inactiveSequences
	public static PrintWriter out;
	
	
	public static void main (String[] args) throws IOException
	{
		activebr = new BufferedReader(new FileReader("ActiveSequences.txt"));
		uniprotbr = new BufferedReader(new FileReader("uniprot_all_data.txt"));
		inactivebr = new BufferedReader(new FileReader("InactiveSequencesTotal.txt"));
		inactiveout = new PrintWriter(new BufferedWriter(new FileWriter("InactiveSequencesTotal.txt")));
		out = new PrintWriter(new BufferedWriter(new FileWriter("TrainingSequences.txt")));
		inactiveout2 = new PrintWriter(new BufferedWriter(new FileWriter("InactiveSequences.txt")));
		
		// formatting tag for the NN program
		out.println("Sequence,Length,InterfaceHydrophobicity,OctanolHydrophobicity,GRAVY," + 
		"TotalCharge,TotalPositiveCharge,TotalNegativeCharge,AveragePositivePosition,AverageNegativePosition," + 
				"Weight,Antimicrobial");
		
		// first read active sequences
		// don't print active sequences with length <5
		String line;
		int[] lengths = new int[63]; // keep track of how many peptides have each length 
        while( (line = activebr.readLine()) != null ){
        	if (line.length() >= 5) {
        		String s = FeatureCalculator.getFeaturesString(line, true);
        		out.println(s);
        		lengths[line.length()]++;
        	}
        }
        for (int i = 0; i<lengths.length; i++)
        	System.out.println(i + ": " + lengths[i]);
        
        ArrayList<String>[] inactiveSequences = new ArrayList[175]; //store inactive sequences by length (length is index)
        for (int i = 0; i<inactiveSequences.length; i++) //initialize the array
        	inactiveSequences[i] = new ArrayList<String>();
        // this loop formats uniprot_all_data.txt, print the sequences to inactiveSequences.txt
        line = uniprotbr.readLine(); // skip the first line
        String s = "";
        int x = 0;
        while( x<6072922 && (line = uniprotbr.readLine()) != null ){ //read 6 million sequences
        	if (line.startsWith(">"))
        	{
        		inactiveout.println(s); // print the sequence to inactiveSequences
        		inactiveSequences[s.length()].add(s);
        		s = "";
        		x++;
        	}
        	else // it's part of a sequence
        	{
        		s += line.trim();
        	}
        }
        for (int i = 0; i< inactiveSequences.length; i++)
        	System.out.println(i + ": " + inactiveSequences[i].size());
        inactiveout.close();
        
		// then write inactive sequences to TrainingSequences.txt randomly
        
        int totalInactiveData = 0;
        for (int i = 5; i <= 62; i++) // for each specific length
        {
        	System.out.println(i);
        	for (int j = 0; j < lengths[i] && inactiveSequences[i].size()>0; j++) // for the peptides you need of a length
        	{
	        	int n = (int)(Math.random()*inactiveSequences[i].size()); // random index
	        	s = FeatureCalculator.getFeaturesString(inactiveSequences[i].get(n), false); //print to TrainingSequences.txt
	        	inactiveout2.println(inactiveSequences[i].get(n));
	        	out.println(s);
	        	System.out.println(s);
	        	inactiveSequences[i].remove(n); // take out the used
	        	totalInactiveData++;
        	}
        }
        System.out.println(totalInactiveData);
		// close writer!
        inactiveout2.close();
        out.close();
	}
}
