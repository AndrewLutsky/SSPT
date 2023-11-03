package main

import (
	"fmt"
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
	windowSize := 7
	// made a slice of slice of CFScore, for each protein in fasta file
	ProteinPredArray := make([][]CFScore, len(proteins))
	for itr, protein := range proteins {
		ProteinPredArray[itr] = ChouFasman(protein, windowSize, parameters, aaIndexMap)
		fmt.Println("\n\n\nPrediction Array:", ProteinPredArray[itr])
		fmt.Println("\n\nFound Helicies:", IdentifyHelicies(ProteinPredArray[itr]))
		fmt.Println("\n\nFound Beta bends:, ", IdentifyTurns(protein, parameters, aaIndexMap))
	}

	// Testing Visualization
	aaSecStruct := TestVisualization(proteins[0], ProteinPredArray[0])
	// call the visualization function
	coordinates := Visualize(aaSecStruct, 0.9, 1.3)
	// write the coordinates to a pdb file
	err := WriteCoordinatesToPDB(coordinates, "pdb_files/testLessPitchMoreDist.pdb")
	if err != nil {
		fmt.Println(err)
	}
}
