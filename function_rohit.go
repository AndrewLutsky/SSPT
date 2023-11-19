// Rohit Nandakumar

package main

import (
	"bufio"
	"log"
	"os"
	"strings"
)

func GenerateFASTAReader(fileName string) *FASTAReader {
	file, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}

	// (Shashank) below line needed to be commented out to make the file non empty. Else, the ReadProteins
	// function will return an empty slice
	// defer file.Close()
	return &FASTAReader{file: file}
}

// ReadProteins reads proteins from given FASTA file
func ReadProteins(fastaReader *FASTAReader) []Protein {
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
		}
	}

	// Add the last protein after the loop finishes
	if currProtein.Identifier != "" {
		proteins = append(proteins, currProtein)
	}
	return proteins
}

// IdentifyBetaSheet is a function that takes as its input a prediction array of CFScores and
// returns an array of BetaSheet objects identified from the sequence
func IdentifyBetaSheet(predArray PredArray) []ABHelixSheet {
	betaSheets := make([]ABHelixSheet, 0)
	n := len(predArray)
	// scan through array until viable starting points are found
	for i := 0; i < n; i++ { // can't range here, have to manually delcare so that we can enforce window skipping
		predPoint := predArray[i]
		if predPoint.Sheet > 1.00 && i < n-3 { //Valid starting position (note to Jon, check if i < n - 3 is right)
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
	var toDelete []int // Slice to keep track of indices of items to be deleted

	// Loop through all items EXCEPT last one
	for i := 0; i < len(allItems)-1; i++ {
		// Loop through all items starting from i + 1 to the very end
		for j := i + 1; j < len(allItems); j++ {
			// checks to see if each ABHelixSheet is of a different type (helix or sheet)
			if allItems[i].typeAB != allItems[j].typeAB {
				// Determine if there is overlap between the two regions
				// Note to group: Just make sure overlap looks right please
				overlap := (allItems[i].EndIndex > allItems[j].StartIndex && allItems[i].EndIndex < allItems[j].EndIndex) || (allItems[j].EndIndex > allItems[i].StartIndex && allItems[j].EndIndex < allItems[i].EndIndex)

				// in this case, there is overlap...
				if overlap {
					// If the score at position i is greater than at position j, mark j for deletion
					if allItems[i].Score > allItems[j].Score {
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
