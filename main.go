package main

import (
	"fmt"
	"strconv"
	//"github.com/fogleman/gg"
)

func main() {

	// Read the parameters file
	parameters := ReadParameters("CFparameters.txt")

	// Read the AAIndex file
	aaIndexMap := ReadAAIndexMap("aa_index_map.txt")

	// Determine input. Then either translate or read from FASTA.
	fmt.Println("Reading from 'FASTA'? Or 'DNA'?")
	var fileType string
	fmt.Scanln(&fileType)
	var proteins []Protein
	// Read the fasta file
	if fileType == "FASTA" {
		reader := GenerateFASTAReader("rcsb_pdb_1UBQ.fasta")
		proteins = ReadProteinsFASTA(reader)
		reader.file.Close() // closing the file after calling the ReadProteins function
	} else if fileType == "DNA" {
		reader := GenerateDNAReader("dna.txt")
		proteins = append(proteins, TranslateDNA(reader))
		reader.file.Close() // closing the file after calling the ReadProteins function
	} else {
		panic("Incorrect reading source entered.")
	}

	// windowsize initialised
	// windowSize := 0
	// made a slice of slice of CFScore, for each protein in fasta file
	ProteinPredArray := make([][]CFScore, len(proteins))
	for itr, protein := range proteins {
		ProteinPredArray[itr] = ChouFasman(protein, parameters, aaIndexMap)
		//fmt.Println("\n\n\nPrediction Array:", ProteinPredArray[itr])

		// predict helices
		helices := IdentifyHelicies(ProteinPredArray[itr])
		fmt.Println("\n\nFound Helicies:", helices)
		// predict beta sheets
		betaSheets := IdentifySheets2(ProteinPredArray[itr])
		fmt.Println("\n\nFound Beta Sheets:", betaSheets)
		// reassign the helices and sheets as appropriate
		reassignedABHelixSheet := AHelicalBSheetAssignment(append(helices, betaSheets...))
		reassignedABHelixSheet = FillGapsInSequence(len(ProteinPredArray), reassignedABHelixSheet)

		fmt.Println("\n\nFound after Reassignment:", reassignedABHelixSheet)
		fmt.Println("\n\nFound Beta bends:, ", IdentifyTurns(protein, parameters, aaIndexMap))

		// VISUALIZATION CODE BELOW

		// make int slice, which would store the information about the secondary structure of each position
		// (1 = helix, 2 = sheet, 3 = loop).
		aaSecStruct := make([]int, len(protein.Sequence))

		// It should be initialized with 3s (as if not helix or sheet, then loop) and we assign helix and sheet below
		for i := 0; i < len(aaSecStruct); i++ {
			aaSecStruct[i] = 3
		}

		// convert the reassignedABHelixSheet slice to aaSecStruct slice
		aaSecStruct = ConvertABHelixSheetToAASecStruct(reassignedABHelixSheet, aaSecStruct)

		// 2D VISUALIZATION
		Make2DPlot(aaSecStruct, "outputs/2d_plots/2DPlot_"+strconv.Itoa(itr)+".png")

		// 3D VISUALIZATION
		Make3DPlot(aaSecStruct, "3d_visualization_resources/1ubq.pdb", "outputs/3d_htmls/1ubq.html")
	}

}
