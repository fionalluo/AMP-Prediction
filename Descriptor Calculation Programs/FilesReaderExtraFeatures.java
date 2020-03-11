/*
 * FilesReaderExtraFeatures adds more information to TrainingSequences.txt for the anti-HIV peptides by adding on 
 * additional features from the ExtraFeatures text files. 
 */
import java.io.*;
import java.util.*;


public class FilesReaderExtraFeatures 
{
	// read a file for HIV sequences, false sequences
	// later, classify sampled sequences
	public static BufferedReader activebr;
	public static BufferedReader uniprotbr; // reads data from uniprot
	public static PrintWriter inactiveout; // prints uniprot sequences to InactiveSequences
	public static BufferedReader inactivebr;
	public static BufferedReader inactivefixedbr; // the ~1000 fixed inactive sequences, randomly generated
	public static PrintWriter inactiveout2; //prints selected 1000 sequences to inactiveSequences
	public static PrintWriter out;
	public static BufferedReader extraFeaturesActivebr; // read in the hardnesses of active sequences from HardnessHIVActive.txt
	public static BufferedReader extraFeaturesInactivebr; // read in hardnesses of the InactiveSequencesFixed file
	
	public static void main (String[] args) throws IOException
	{
		activebr = new BufferedReader(new FileReader("ActiveSequences.txt"));
		uniprotbr = new BufferedReader(new FileReader("uniprot_all_data.txt"));
		inactivebr = new BufferedReader(new FileReader("InactiveSequencesTotal.txt"));
		inactiveout = new PrintWriter(new BufferedWriter(new FileWriter("InactiveSequencesTotal.txt")));
		out = new PrintWriter(new BufferedWriter(new FileWriter("TrainingSequences.txt")));
		inactiveout2 = new PrintWriter(new BufferedWriter(new FileWriter("InactiveSequences.txt")));
		extraFeaturesActivebr = new BufferedReader(new FileReader("ExtraFeaturesHIVActive.txt"));
		extraFeaturesInactivebr = new BufferedReader(new FileReader("ExtraFeaturesHIVInactive.txt"));
		inactivefixedbr = new BufferedReader(new FileReader("InactiveSequencesFixed.txt"));
		
		// formatting tag for the NN program
		String tags = "Sequence,Length,InterfaceHydrophobicity,OctanolHydrophobicity,GRAVY," + 
		"TotalCharge,TotalPositiveCharge,TotalNegativeCharge,AveragePositivePosition,AverageNegativePosition," + 
				"Weight,Antimicrobial";
		//tags += ",hardness,eo,sigma";
		tags += ",hmol,hmol_pos,hmol_neg,smol_pos,smol_neg,ave_smol,ave_hmol_pos,ave_hmol_neg,ave_smol_pos,ave_smol_neg,"+
				"hmol_pos_largest,hmol_neg_largest,hmol_pos_smallest,hmol_neg_smallest,eo_pos,sigma,";
		tags += "stericmol_largest,stericmol_smallest,stericatom_largest,stericatom_smallest";
		out.println(tags);
		
		
		// first read active sequences
		// don't print active sequences with length <5
		String line;
		int[] lengths = new int[63]; // keep track of how many peptides have each length 
		int x = 0;
        while( (line = activebr.readLine()) != null ){
        	String hardness = extraFeaturesActivebr.readLine().trim();
        	if (line.length() >= 5) {
        		String s = FeatureCalculator.getFeaturesString(line, true);
        		out.println(s + ", " + hardness);
        		lengths[line.length()]++;
        	}
        	x++;
        	System.out.println(x);
        }
        for (int i = 0; i<lengths.length; i++)
        	System.out.println(i + ": " + lengths[i]);
        
        while( (line = inactivefixedbr.readLine()) != null ){
        	String hardness = extraFeaturesInactivebr.readLine().trim();
        	if (line.length() >= 5) {
        		String s = FeatureCalculator.getFeaturesString(line, false);
        		out.println(s + ", " + hardness);
        	}
        }      
        
		// close writer!
        inactiveout2.close();
        out.close();
	}
}
