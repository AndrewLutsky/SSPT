package main

import (
	"fmt"
)

func main() {

	// Read the parameters file
	parameters := ReadParameters("dummy_parameters.txt")
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

	// windowsize initialised
	windowSize := 7

	for _, protein := range proteins {
		predArray := ChouFasman(protein, windowSize, parameters, aaIndexMap)
		fmt.Println("\n\n\nPrediction Array:", predArray)
		fmt.Println("\n\nFound Helicies:", IdentifyHelicies(predArray))
	}

}
