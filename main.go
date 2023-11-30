package main

import (
	"fmt"
)

func main() {
	/*cifReader := GenerateCIFReader("cif_files/1fez.cif")
	test, _ := ReadCIFToFasta(cifReader) */

	//fmt.Println(GORMethodInput("GorParams/alpha_helix.txt"))

	// Read the parameters file
	parameters := ReadParameters("CFparameters.txt")

	// Read the AAIndex file
	aaIndexMap := ReadAAIndexMap("aa_index_map.txt")

	// Determine input. Then either translate or read from FASTA.
	// There will three types of input: an array of uniprot IDs, a FASTA file, or a DNA sequence.
	fmt.Println("Are you i)supplying an 'array' of Uniprots OR ii)Reading from 'FASTA' OR iii)Reading from 'DNA'?")
	var fileType string
	fmt.Scanln(&fileType)
	var proteins []Protein
	var outputNames []string
	var uniprotIDs []string

	// NOTE: THE 3D VISUALIZATION IS ONLY AVAILABLE FOR THE ARRAY INPUT, AS ONLY THEN CAN WE GET THE PDB FILE

	// Make the proteins slice according to the input
	if fileType == "array" {
		fmt.Println("Enter the Ensembl IDs separated by commas:")
		var ensemblInput string

		// ask the user to enter the uniprot IDs
		fmt.Scanln(&ensemblInput)

		// make a slice of uniprot IDs using the uniportInput string
		ensemblIDs := MakeEnsemblArray(ensemblInput)

		// make a slice of proteins using the ensembl IDs
		proteins = EnsemblToProteinSlice(ensemblIDs)

		// convert the ensembl IDs to uniprot IDs
		uniprotIDs = EnsemblToUniprotSlice(ensemblIDs)

		// need to download the pdb files for the uniprot IDs
		DownloadPDBFiles(uniprotIDs)

		// also make a outputNames array which will store the outputNames of the output files
		outputNames = ensemblIDs

	} else if fileType == "FASTA" {
		// ask the user to enter the name of the fasta file
		fmt.Println("Enter the name of the FASTA file:")
		var fastaInput string
		fmt.Scanln(&fastaInput)

		// make a slice of proteins using the fasta file
		reader := GenerateFASTAReader(fastaInput)
		proteins = ReadProteinsFASTA(reader)
		reader.file.Close() // closing the file after calling the ReadProteins function

		// make a outputNames array which will store the outputNames of the output files
		outputNames = MakeOutputNames(proteins)

	} else if fileType == "DNA" {
		// ask the user to enter the name of the DNA file
		fmt.Println("Enter the name of the DNA file:")
		var dnaInput string
		fmt.Scanln(&dnaInput)

		// make a slice of proteins using the dna file
		reader := GenerateDNAReader(dnaInput)
		proteins = append(proteins, TranslateDNA(reader))
		reader.file.Close() // closing the file after calling the ReadProteins function

		// make a outputNames array (specific to DNA type run) which will store the outputNames of the output files
		outputNames = MakeDNAOutputNames(proteins)

	} else {
		panic("Incorrect reading source entered.")
	}

	// made a slice of slice of CFScore, for each protein in fasta file
	ProteinPredArray := make([][]CFScore, len(proteins))
	//ProteinPredGorArray := make([][]CFScore, len(proteins))
	for itr, protein := range proteins {
		ProteinPredArray[itr] = ChouFasman(protein, parameters, aaIndexMap)
		fmt.Println(protein, parameters, aaIndexMap)

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
		// reassignedABHelixSheet = IdentifyTurns(protein, parameters, aaIndexMap, reassignedABHelixSheet)
		// fmt.Println("Completed secondary structure assignment using CF!")

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
		// name according to the outputNames array
		Make2DPlot(aaSecStruct, "outputs/2d_plots/2DPlot_"+outputNames[itr]+".png")

		// 3D VISUALIZATION (only if filetype is 'array')
		if fileType == "array" {
			// name output according to the outputNames array and also the pdb will be existing in the
			// 3d_visualization_resources folder with the same name
			Make3DPlot(aaSecStruct, "3d_visualization_resources/"+uniprotIDs[itr]+".pdb", "outputs/3d_htmls/3DPlot_"+outputNames[itr]+".html")

		}
	}
	/*
		for itr, protein := range proteins {
			ProteinPredGorArray[itr] = GORPrediction(protein, parameters, aaIndexMap)
			fmt.Println(protein, parameters, aaIndexMap)
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
			reassignedABHelixSheet = IdentifyTurns(protein, parameters, aaIndexMap, reassignedABHelixSheet)
			fmt.Println("Completed secondary structure assignment using CF!")

		}
	*/
}
