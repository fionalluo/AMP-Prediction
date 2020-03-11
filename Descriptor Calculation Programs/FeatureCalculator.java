/**

This program calculates the features of a peptide given a sequence and prints it to a csv file
Key:

Sequence,Length,InterfaceHydrophobicity,OctanolHydrophobicity,GRAVY,
TotalCharge,TotalPositiveCharge,TotalNegativeCharge,AveragePositivePosition,AverageNegativePosition,
Weight,Antimicrobial

 */
import java.io.*;

public class FeatureCalculator 
{	
	public static void main(String[] args) throws FileNotFoundException
	{
		//System.out.println(getFeaturesString("GLFDIIKKIAESF", true));
		System.out.println(getFeaturesString("ACDEFGHIKLMNPQRSTVWY", true));
	}
	
	// activity is 1 if it's active, 0 if inactive
	public static String getFeaturesString(String seq, boolean activity)
	{
		String s = seq;
		s += ", " + getLength(seq);
		s += ", " + getInterfaceHydrophobicity(seq);
		s += ", " + getOctanolHydrophobicity(seq);
		s += ", " + getGRAVY(seq);
		s += ", " + getTotalCharge(seq);
		s += ", " + getTotalPositiveCharge(seq);
		s += ", " + getTotalNegativeCharge(seq);
		s += ", " + getAveragePositivePosition(seq);
		s += ", " + getAverageNegativePosition(seq);
		s += ", " + getWeight(seq);
		if (activity == true) // if active, 1
			s += ", " + 1;
		else
			s += ", " + 0; // if inactive, 0
		return s;
	}
	
	// if NO ACTIVITY IS GIVEN
	public static String getFeaturesString(String seq)
	{
		String s = seq;
		s += ", " + getLength(seq);
		s += ", " + getInterfaceHydrophobicity(seq);
		s += ", " + getOctanolHydrophobicity(seq);
		s += ", " + getGRAVY(seq);
		s += ", " + getTotalCharge(seq);
		s += ", " + getTotalPositiveCharge(seq);
		s += ", " + getTotalNegativeCharge(seq);
		s += ", " + getAveragePositivePosition(seq);
		s += ", " + getAverageNegativePosition(seq);
		s += ", " + getWeight(seq);
		return s;
	}
	
	public static int getLength(String sequence)
	{
		return sequence.length();
	}
	
	// calculate using Wimbley-White hydrophobicity scale, experimentally determined values
	// source: https://en.wikipedia.org/wiki/Hydrophobicity_scales#Wimley%E2%80%93White_whole_residue_hydrophobicity_scales
	// source: http://www.tulane.edu/~biochem/faculty/facfigs/NSB.pdf
	public static double getInterfaceHydrophobicity(String sequence)
	{
		double hydrophobicity = 0;
		for (int i = 0; i<sequence.length(); i++)
		{
			char c = sequence.charAt(i);
			switch (c) 
			{
				case 'A': hydrophobicity += 0.17; break;	
				case 'R': hydrophobicity += 0.81; break;	
				case 'N': hydrophobicity += 0.42; break;	
				case 'D': hydrophobicity += 1.23; break;	 //asp- use this
				//case 'D': hydrophobicity += -0.07; break;	 //asp0 
				case 'C': hydrophobicity += -0.24; break;	
				case 'Q': hydrophobicity += 0.58; break;	 
				case 'E': hydrophobicity += 2.02; break;	 //glu- use this
				//case 'E': hydrophobicity += -0.01; break;	 //glu0
				case 'G': hydrophobicity += 0.01; break;	
				//case 'H': hydrophobicity += 0.96; break;	 //his+ 
				case 'H': hydrophobicity += 0.17; break;	 //his0 use this FOR NOW!!!!!!
				case 'I': hydrophobicity += -0.31; break;	
				case 'L': hydrophobicity += -0.56; break;	
				case 'K': hydrophobicity += 0.99; break;	
				case 'M': hydrophobicity += -0.23; break;	
				case 'F': hydrophobicity += -1.13; break;	
				case 'P': hydrophobicity += 0.45; break;	
				case 'S': hydrophobicity += 0.13; break;	
				case 'T': hydrophobicity += 0.14; break;	
				case 'W': hydrophobicity += -1.85; break;	
				case 'Y': hydrophobicity += -0.94; break;	
				case 'V': hydrophobicity +=	0.07; break;	
			}
		}
		hydrophobicity = (int)(hydrophobicity*100)/100.0;
		return hydrophobicity;
	}
	
	// calculate using Wimbley-White hydrophobicity scale, experimentally determined values
	// source: https://en.wikipedia.org/wiki/Hydrophobicity_scales#Wimley%E2%80%93White_whole_residue_hydrophobicity_scales
	// source: http://www.tulane.edu/~biochem/faculty/facfigs/NSB.pdf
	public static double getOctanolHydrophobicity(String sequence)
	{
		double hydrophobicity = 0;
		for (int i = 0; i<sequence.length(); i++)
		{
			char c = sequence.charAt(i);
			switch (c) 
			{
				case 'A': hydrophobicity += 0.50; break;	
				case 'R': hydrophobicity += 1.81; break;	
				case 'N': hydrophobicity += 0.85; break;	
				case 'D': hydrophobicity += 3.64; break;	 //asp- use this
				//case 'D': hydrophobicity += 0.43; break;	 //asp0 
				case 'C': hydrophobicity += -0.02; break;	
				case 'Q': hydrophobicity += 0.77; break;	 
				case 'E': hydrophobicity += 3.63; break;	 //glu- use this
				//case 'E': hydrophobicity += 0.11; break;	 //glu0
				case 'G': hydrophobicity += 1.15; break;	
				//case 'H': hydrophobicity += 2.33; break;	 //his+ 
				case 'H': hydrophobicity += 0.11; break;	 //his0 use this FOR NOW
				case 'I': hydrophobicity += -1.12; break;	
				case 'L': hydrophobicity += -1.25; break;	
				case 'K': hydrophobicity += 2.80; break;	
				case 'M': hydrophobicity += -0.67; break;	
				case 'F': hydrophobicity += -1.71; break;	
				case 'P': hydrophobicity += 0.14; break;	
				case 'S': hydrophobicity += 0.46; break;	
				case 'T': hydrophobicity += 0.25; break;	
				case 'W': hydrophobicity += -2.09; break;	
				case 'Y': hydrophobicity += -0.71; break;	
				case 'V': hydrophobicity +=	-0.46; break;	
			}
		}
		hydrophobicity = (int)(hydrophobicity*100)/100.0;
		return hydrophobicity;
	}
	
	// Formula: sum of hydropathy values divided by number of residues in sequence
	public static double getGRAVY(String s)
	{
		double hydrophobicity = 0;
		for (int i = 0; i<s.length(); i++)
		{
			char c = s.charAt(i);
			switch (c) 
			{
				case 'A': hydrophobicity += 1.8; break;	
				case 'R': hydrophobicity += -4.5; break;	
				case 'N': hydrophobicity += -3.5; break;	
				case 'D': hydrophobicity += -3.5; break;	 //asp
				case 'C': hydrophobicity += 2.5; break;	
				case 'Q': hydrophobicity += -3.5; break;	 
				case 'E': hydrophobicity += -3.5; break;	 //glu
				case 'G': hydrophobicity += -0.4; break;	
				case 'H': hydrophobicity += -3.2; break;	 //his
				case 'I': hydrophobicity += 4.5; break;	
				case 'L': hydrophobicity += 3.8; break;	
				case 'K': hydrophobicity += -3.9; break;	
				case 'M': hydrophobicity += 1.9; break;	
				case 'F': hydrophobicity += 2.8; break;	
				case 'P': hydrophobicity += -1.6; break;	
				case 'S': hydrophobicity += -0.8; break;	
				case 'T': hydrophobicity += -0.7; break;	
				case 'W': hydrophobicity += -0.9; break;	
				case 'Y': hydrophobicity += -1.3; break;	
				case 'V': hydrophobicity +=	4.2; break;	
			}
		}
		return hydrophobicity/s.length();
	}
	
	
	// K and R are positive, E and D are negative
	// Histinine has pka of 6-6.5, so counted it as having 0.1 positive charge [see source below]
	// source: https://pepcalc.com/ppc.php
	// Limitations: Discounts polarities
	public static double getTotalCharge(String seq)
	{
		String pos = "KR";
		String neg = "ED";
		double charge = 0;
		for (int i = 0; i<seq.length(); i++)
		{
			if (pos.indexOf(seq.charAt(i)) != -1)
				charge++;
			else if (neg.indexOf(seq.charAt(i)) != -1)
				charge --;
			else if (seq.charAt(i) == 'H')
				charge += 0.1;
		}
		charge = (int)(charge*10)/10.0; //fix rounding errors
		return charge;
	}
	
	// total positive charge
	public static double getTotalPositiveCharge(String s)
	{
		String pos = "KR";
		double charge = 0.0;
		for (int i = 0; i<s.length(); i++)
		{
			if (pos.indexOf(s.charAt(i)) != -1)
				charge++;
			else if (s.charAt(i) == 'H')
				charge += 0.1;
		}
		charge = (int)(charge*10)/10.0; //fix rounding errors
		return charge;
	}
	
	// total negative charge
	public static double getTotalNegativeCharge(String s)
	{
		String neg = "ED";
		double charge = 0;
		for (int i = 0; i<s.length(); i++)
		{
			if (neg.indexOf(s.charAt(i)) != -1)
				charge --;
		}
		charge = (int)(charge*10)/10.0; //fix rounding errors
		return charge;
	}
	
	// average position of positive charge
	// formula: (xw1+xw2+xw3)/(w1+w2+w3) <- weight is magnitude of charge
	// if no positive charges, return the midpoint
	public static double getAveragePositivePosition(String s)
	{
		String pos = "KR";
		double weightedSum = 0.0;
		double totalCharge = 0.0;
		for (int i = 0; i<s.length(); i++)
		{
			if (pos.indexOf(s.charAt(i)) != -1)
			{
				weightedSum += i; // add the position * 1
				totalCharge ++;
			}
			else if (s.charAt(i) == 'H')
			{
				weightedSum += i * 0.1;
				totalCharge += 0.1;
			}
		}
		if (totalCharge != 0)
			return weightedSum/totalCharge;
		return s.length()/2.0;
	}
	
	// average position of negative charge
	// formula: (xw1+xw2+xw3)/(w1+w2+w3) <- weight is magnitude of charge
	public static double getAverageNegativePosition(String s)
	{
		String neg = "ED";
		double weightedSum = 0.0;
		double totalCharge = 0.0;
		for (int i = 0; i<s.length(); i++)
		{
			if (neg.indexOf(s.charAt(i)) != -1)
			{
				weightedSum += i; // add the position * 1
				totalCharge ++;
			}
		}
		if (totalCharge != 0)
			return weightedSum/totalCharge;
		return s.length()/2.0;
	}
	
	// weight measured in DA's, aka g/mol
	// Formula: Add weight of each amino acid, subtract (N-1)*17.99 <- the weight of each H2O bond
	// Disregards weight of terminus
	// source: https://www.selleckchem.com/peptide-calculator.html (molecular weight calculator)
	// source2: https://pepcalc.com/notes.php?mw
	public static double getWeight(String sequence)
	{
		double weight = 0;
		for (int i = 0; i<sequence.length(); i++)
		{
			char c = sequence.charAt(i);
			switch (c) 
			{
				case 'A': weight +=	89.09; break;
				case 'C': weight +=	121.16; break;
				case 'D': weight += 133.10; break;
				case 'E': weight += 147.13; break;
				case 'F': weight += 165.19; break;
				case 'G': weight += 75.07; break;
				case 'H': weight += 155.16; break;
				case 'I': weight += 131.18; break;
				case 'K': weight += 146.19; break;
				case 'L': weight += 131.18; break;
				case 'M': weight += 149.21; break;
				case 'N': weight += 132.12; break;
				case 'P': weight += 115.13; break;
				case 'Q': weight += 146.15; break;
				case 'R': weight += 174.20; break;
				case 'S': weight += 105.09; break;
				case 'T': weight += 119.12; break;
				case 'V': weight += 117.15; break;
				case 'W': weight += 204.23; break;
				case 'Y': weight += 181.19; break;
			}
		}
		weight -= (sequence.length()-1)*17.99; // subtract weight of N-1 H2O bonds
		// round to hundredths place to avoid rounding errors, and keep 2 decimal places
		weight = (int)(weight*100)/100.0;
		return weight;
	}
	
	public static int getAcidic(String sequence)
	{
		return 0;
	}
	public static double getBasic(String sequence)
	{
		return 0.0;
	}
}







