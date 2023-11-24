// Rohit Nandakumar

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
)

func GenerateFASTAReader(fileName string) FASTAReader {
	file, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}

	// (Shashank) below line needed to be commented out to make the file non empty. Else, the ReadProteins
	// function will return an empty slice
	// defer file.Close()
	return FASTAReader{file: file}
}

// ReadProteins reads proteins from given FASTA file
func ReadProteinsFASTA(fastaReader FASTAReader) []Protein {
	var proteins []Protein
	var currProtein Protein
	scanner := bufio.NewScanner(fastaReader.file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if currProtein.Identifier != "" {
				proteins = append(proteins, currProtein)
				currProtein = Protein{}
			}
			currProtein.Identifier = strings.TrimPrefix(line, ">")

		} else if line != "" {
			currProtein.Sequence += line
		} else if line == "" {
			proteins = append(proteins, currProtein)
			currProtein = Protein{}
		}
	}

	proteins = append(proteins, currProtein)
	currProtein = Protein{}

	return proteins
}

func GenerateCIFReader(fileName string) *CIFReader {
	file, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}

	// (Shashank) below line needed to be commented out to make the file non empty. Else, the ReadProteins
	// function will return an empty slice
	// defer file.Close()
	return &CIFReader{file: file}
}

// ReadCIFToFasta reads protein sequences from a CIF file and writes them to a FASTA format.
func ReadCIFToFasta(cifReader *CIFReader) (Protein, error) {
	// Define the identifier to look for in the CIF file.
	identifier := "_entity_poly.pdbx_seq_one_letter_code"

	protein := Protein{}

	// Create a scanner to read the file line by line.
	scanner := bufio.NewScanner(cifReader.file)

	// Flags to keep track of whether the sequence has been foundand whether we've encountered the second semicolon.
	sequenceFound := false
	secondSemicolon := false

	// Variable to store the extracted sequence.
	var sequence string

	// Loop over each line in the file.
	for scanner.Scan() {
		line := scanner.Text()

		// Parse the protein entry id into the protein
		if strings.Contains(line, "_entry.id") {
			parsedString := strings.TrimPrefix(line, "_entry.id")
			parsedString = strings.TrimSpace(parsedString)
			protein.Identifier = parsedString
		}

		// Check if the line contains the identifier.
		// If it does, start capturing the sequence.
		if strings.Contains(line, identifier) {
			sequenceFound = true
		}

		// If we've found the sequence start and not yet encountered the second semicolon,
		// trim any leading semicolon and append the line to the sequence.
		if strings.HasPrefix(line, ";") && sequenceFound && !secondSemicolon {
			parsedString := strings.TrimPrefix(line, ";")
			sequence += strings.TrimSpace(parsedString)
		}

		// If we encounter a semicolon and have already started capturing the sequence,
		// set the flag indicating that we've reached the end of the sequence.
		if strings.HasPrefix(line, ";") && sequence != "" {
			secondSemicolon = true
		}
	}

	// Put sequence into protein
	protein.Sequence = sequence

	// Check for any errors that occurred during the file reading.
	if err := scanner.Err(); err != nil {
		return Protein{}, err
	}

	// Return the protein which will have an identifier and a protein.
	return protein, nil
}

// IdentifyBetaSheet is a function that takes as its input a prediction array of CFScores and
// returns an array of BetaSheet objects identified from the sequence
func IdentifyBetaSheet(predArray PredArray) []ABHelixSheet {
	betaSheets := make([]ABHelixSheet, 0)
	n := len(predArray)
	// scan through array until viable starting points are found
	for i := 0; i < n; i++ { // can't range here, have to manually delcare so that we can enforce window skipping
		predPoint := predArray[i]
		if predPoint.Sheet > 1.05 && i < n-3 { //Valid starting position (note to Jon, check if i < n - 3 is right)
			j := i
			helixSum := predPoint.Helix // Start tracking for average
			sheetSum := predPoint.Sheet
			for j = i + 1; j < len(predArray); j++ { //Scan forward until stopping point is found or end of sequence is reached.
				helixSum += predArray[j].Helix
				sheetSum += predArray[j].Sheet
				if predPoint.Sheet < 1.00 { // End of the potential sheet
					break
				}
			}
			// At this point, we've either reached the end of the potential sheet or hit the end of the sequence.
			// Now, validate the potential sheet.
			length := j - i
			avgHelixScore := helixSum / float64(length)
			avgSheetScore := sheetSum / float64(length)
			if avgSheetScore > avgHelixScore && avgSheetScore > 1.05 { // We have a valid sheet! Create the sheet and add it to the array, then skip to the end of the sheet.
				var newBetaSheet ABHelixSheet
				newBetaSheet.StartIndex = i
				newBetaSheet.EndIndex = j
				newBetaSheet.Score = avgSheetScore
				newBetaSheet.typeAB = "sheet"
				//return betaSheet
				betaSheets = append(betaSheets, newBetaSheet)
			}
			i = j // ALWAYS jump the window once you have found a helix. Chou-Fasman says that you must also extend the search window backwards when calculating helicies
			// so by always jumping the window we never have to consider this.
		}
	}
	return betaSheets
}

// AHelicalBSheetAssignment identifies overlapping regions of alpha-helices and beta-sheets,
// and removes the region with the lower score.
func AHelicalBSheetAssignment(allItems []ABHelixSheet) []ABHelixSheet {
	//fmt.Println("All Items for Overlap Deletion:", allItems)
	var toDelete []int // Slice to keep track of indices of items to be deleted

	// Loop through all items EXCEPT last one
	for i := 0; i < len(allItems)-1; i++ {
		element1 := allItems[i]
		// Loop through all items starting from i + 1 to the very end
		for j := 0; j < len(allItems); j++ {
			element2 := allItems[j]
			// checks to see if each ABHelixSheet is of a different type (helix or sheet)
			if element1 != element2 {
				// Determine if there is overlap between the two regions
				// Note to group: Just make sure overlap looks right please
				overlap := (element1.EndIndex >= element2.StartIndex && element1.EndIndex <= element2.EndIndex) || (element2.EndIndex >= element1.StartIndex && element2.EndIndex <= element1.EndIndex)
				//fmt.Println("Comparing Elem1:\n", element1, "\nAnd Elem2:\n", element2, "\nIs there overlap?", overlap)
				// in this case, there is overlap...
				if overlap {
					// If the score at position i is greater than at position j, mark j for deletion
					if element1.Score > element2.Score {
						toDelete = append(toDelete, j)
					} else {
						// If the score at position j is greater or equal, mark i for deletion
						toDelete = append(toDelete, i)
					}
				}
			}
		}
	}

	// Call a separate function to actually perform the deletion of marked items
	// We do this because any deletions inside the nested for loop that we had
	// would cause all sorts of problems
	allItems = DeleteItemsFromABHelixSheet(allItems, toDelete)

	return allItems
}

// DeleteItemsFromABHelixSheet creates a new slice excluding the items marked for deletion.
func DeleteItemsFromABHelixSheet(allItems []ABHelixSheet, toDelete []int) []ABHelixSheet {
	var newItems []ABHelixSheet // Slice to hold the items not marked for deletion

	// SUGGESTION: ADD "DELETEME" FIELD TO THE ITEMS, THEN JUST DELETE OBJECTS IN THE SLICE WITH THAT FLAG
	// Loop through all items in the original slice
	for i, item := range allItems {
		delete := false
		// Check if the current index is in the toDelete slice
		for _, deleteIndex := range toDelete {
			if i == deleteIndex {
				delete = true // Mark the item for deletion
				break
			}
		}
		// If the item is not marked for deletion, append it to the newItems slice
		if !delete {
			newItems = append(newItems, item)
		}
	}
	return newItems
}

// function to open a file named output, and write the contents of the toEdit slice to it
func WriteToFileString(output string, toEdit string) {
	// Open the file for writing
	file, err := os.Create(output)
	if err != nil {
		return
	}
	defer file.Close()

	// Write the contents of the toEdit slice to the file
	fmt.Fprintln(file, toEdit)
}
