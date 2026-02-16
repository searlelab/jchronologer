Cartographer (Python) is a PyTorch ML model similar in structure to Chronologer, except that instead of predicting retention time, it predicts fragmentation patterns for peptides sequences, also given an NCE and a charge state. On this computer, Cartographer exists as part of the Encyclopydia project (/Users/searle.brian/Documents/projects/encyclopydia). The repo has been changed slightly to work with the Wilhelmlab datasets (described in PARQUET_DATA_STRUCTURES.md) where the changes to training are noted in TRAINING_CHANGES.md. The main reason for this was to support a wider degree modifications, with details in PARQUET_MODIFICATIONS.md. 

The new Parquet data required some changes to inputs and outputs, which the above docs discuss. First, for inputs, modifications are described using Unimod, rather than mass values (PARQUET_MODIFICATIONS.md gives a conversion table). For output predictions, a spectrum is represented by a 174-dimensional vector (y/b ions, 3 charges, 29 fragment ions) and ordered as follows: y1 (1+), y1 (2+), y1 (+3), b1 (1+), b1 (2+), b1 (2+), y2 (1+) and so on. Fragment ion intensity values at impossible dimensions (that is, y20 for a 7-mer) are set to −1. This will require reshaping the vector to produce y-ion and b-ion ladders.

I would like to have JChronologer (this project) provide a similar adaption Java layer for this modified version of Cartographer. Input should still allow for Chronologer-styled mass encoding, but also allow for Unimod-based encoding as well. Converters can adjust the sequences as necessary. Note, the Unimod approach can use compound modifications, which aren't cleanly supported by the mass encoding approach. For example, "[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]" shows a TMT modified peptide where the n-terminus, and both lysines (K) are modified. With the mass encoding approach (peptideModSeq), this peptide might be encoded as "K[458.325864]PGLAITFAK[229.162932]" or "[229.162932]K[229.162932]PGLAITFAK[229.162932]". Internally, let's adjust to Unimod encoding since it is more explicit, but we need converters for peptideModSeq for now for backwards compatability.

We want to use both Chronologer and Cartographer models to be used in tandem to predict a complete peptide library entry given a peptide sequence, NCE, and charge as input. The input should allow for a range of NCE and charge pairs so that multiple forms of the fragmentation can be predicted in the same query (reusing the retention time prediction for the peptide sequence). For now, let's not adjust the RT prediction CLI frontend, but instead build sufficient real-world unit-tests to make sure the system is performing correctly.

Here is a skeleton for how ChronologerLibraryEntry should be designed:
``` java
public class ChronologerLibraryEntry {
	private final String unimodPeptideSequence;
	private final byte precursorCharge;
	private final double precursorNCE;
	
	private final double precursorMZ;
	private final float retentionTimeInSeconds;
	private final double[] massArray;
	private final float[] intensityArray;
	private final String[] ionTypeArray;
	
	// constructors, getters, etc...
	
	public String getPeptideModSeq() {
	   // return a converted sequence...
	}
}
```
where ionTypeArray is a simple String that indicates the ion type and index (e.g., "2+y20" for the 2+ form of the y20 ion).

Please ask me any questions you may have as you plan so we can get this right on the first go.