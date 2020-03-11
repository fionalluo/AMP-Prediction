/*
 * SampleReader.java is used to read, format, and print out the descriptors for the LSTM generated
 * sequences. 
 */
import java.io.*;
import java.util.*;


public class SampleReader 
{
	public static BufferedReader br;
	public static PrintWriter out;
	public static BufferedReader extraFeaturesbr; // read in the hardnesses of active sequences from HardnessHIVActive.tx
	
	public static void main (String[] args) throws IOException
	{
		br = new BufferedReader(new FileReader("sampled_sequences_AMP.txt"));
		out = new PrintWriter(new BufferedWriter(new FileWriter("samples_AMP_formatted.txt")));
		extraFeaturesbr = new BufferedReader(new FileReader("ExtraFeaturesAMPSamples.txt"));
		
		// formatting tag for the NN program
		String tags = "Sequence,Length,InterfaceHydrophobicity,OctanolHydrophobicity,GRAVY," + 
		"TotalCharge,TotalPositiveCharge,TotalNegativeCharge,AveragePositivePosition,AverageNegativePosition," + 
				"Weight,Antimicrobial";
		//tags += ",hardness,eo,sigma";
		tags += ",hmol,hmol_pos,hmol_neg,smol_pos,smol_neg,ave_smol,ave_hmol_pos,ave_hmol_neg,ave_smol_pos,ave_smol_neg,"+
				"hmol_pos_largest,hmol_neg_largest,hmol_pos_smallest,hmol_neg_smallest,eo_pos,sigma,";
		tags += "stericmol_largest,stericmol_smallest,stericatom_largest,stericatom_smallest";
		out.println(tags);
		
		
		// don't print active sequences with length <5
		String line;
		while( (line = br.readLine()) != null ){
        	String hardness = extraFeaturesbr.readLine().trim();
        	if (line.length() >= 5) {
        		String s = FeatureCalculator.getFeaturesString(line, true);
        		out.println(s + ", " + hardness);
        	}
        }
		// close writer!
        out.close();
	}
}
