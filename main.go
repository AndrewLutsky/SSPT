package main

import (
	"fmt"
	//"github.com/fogleman/gg"
)

func main() {

	// Read the parameters file
	parameters := ReadParameters("CFparameters.txt")
	// print what was read
	fmt.Println(parameters)

	// Read the AAIndex file
	aaIndexMap := ReadAAIndexMap("aa_index_map.txt")
	// print what was read
	fmt.Println(aaIndexMap)

	//seq := "MKTLLLTLTLACAGTTAATN"

	// generate
	reader := GenerateFASTAReader("sample_fasta_file.fasta")
	proteins := ReadProteins(reader)
	reader.file.Close() // closing the file after calling the ReadProteins function

	// windowsize initialised
	windowSize := 1
	// made a slice of slice of CFScore, for each protein in fasta file
	ProteinPredArray := make([][]CFScore, len(proteins))
	for itr, protein := range proteins {
		ProteinPredArray[itr] = ChouFasman(protein, windowSize, parameters, aaIndexMap)
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
	}

	// Testing Visualization
	aaSecStruct := TestVisualization(proteins[0], ProteinPredArray[0])
	// call the visualization function (2D)
	Make2DPlot(aaSecStruct, "2d_plots/2DPlot.png")
	// call the visualization function
	// coordinates := Visualize(aaSecStruct, 0.9, 1.3)
	// write the coordinates to a pdb file
	// err := WriteCoordinatesToPDB(coordinates, "pdb_files/testLessPitchMoreDist.pdb")

}
