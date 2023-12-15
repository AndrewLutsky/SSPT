package main

import (
	"bufio"
	"fmt"
	"io/fs"
	"log"
	"os"
	"strconv"
)

//Function to identify bends in the protein.

/*
To identify a bend at residue number j, calculate the following value
p(t) = f(j)f(j+1)f(j+2)f(j+3) where the f(j+1) value for the j+1 residue is used,
the f(j+2) value for the j+2 residue is used and the f(j+3) value for the j+3 residue is used.
If: (1) p(t) > 0.000075; (2) the average value for P(turn) > 1.00 in the tetrapeptide;
and (3) the averages for the tetrapeptide obey the inequality P(a-helix) < P(turn) > P(b-sheet),
then a beta-turn is predicted at that location.
*/

// Function that identifies turns given a particular protein structure, the parameters, an index map,
// and a given ABHelix Sheet structure. It only finds turns in beta sheets.
// Outputs: A new ABHelix Sheet structure with modified Beta Sheets that are interrupted with turns.
// For example given a Sheet that goes from 1 to 10 and a turn that has been identified at 5, the
// returned structure would be a sheet from 1 to 4, a turn at 5, and a sheet from 6 to 10.
func IdentifyTurns(protein Protein, parameters [][]float64, aaIndexMap map[rune]int, secStruc []ABHelixSheet) []ABHelixSheet {
	//Create a slice of turn data types
	turns := make([]Turn, 0)

	//Range through the protein sequence up until the last three(no parameters for these variables)
	for i := 0; i < len(protein.Sequence)-3; i++ {
		//Step 1 calculate the product of f(j) * f(j+1) * f(j+2) * f(j+3)

		//Get the parameters for the next four amino acids

		//J'th amino acid index.
		aaIndex := aaIndexMap[rune(protein.Sequence[i])]
		//J + 1 amino acid index.
		aaIndexPlusOne := aaIndexMap[rune(protein.Sequence[i+1])]
		//J+2 amino acid index.
		aaIndexPlusTwo := aaIndexMap[rune(protein.Sequence[i+2])]
		//J+3 amino acid index.
		aaIndexPlusThree := aaIndexMap[rune(protein.Sequence[i+3])]

		//Create the parameter value for those indeces using the parameters look up table.
		p_t := parameters[aaIndex][3]
		p_t *= parameters[aaIndexPlusOne][4]
		p_t *= parameters[aaIndexPlusTwo][5]
		p_t *= parameters[aaIndexPlusThree][6]

		//Create a boolean to see if it satisfies the first condition.
		satisfiesCond1 := p_t > 0.000075

		//Step 2 calculate the average of P(turn) for the tetrapeptide.
		//Sum all the parameter values for the tetrapeptide.
		avgTurn := parameters[aaIndex][2]
		avgTurn += parameters[aaIndexPlusOne][2]
		avgTurn += parameters[aaIndexPlusTwo][2]
		avgTurn += parameters[aaIndexPlusThree][2]

		//Create another boolean to see if it satisfies the second condition that the
		//average is greater than one.
		satisfiesCond2 := avgTurn/4.0 > 1.0

		//Step 3 - Check Inequality

		//Finds the P(alpha) for the tetrapeptide.
		avgAlpha := parameters[aaIndex][0]
		avgAlpha += parameters[aaIndexPlusOne][0]
		avgAlpha += parameters[aaIndexPlusTwo][0]
		avgAlpha += parameters[aaIndexPlusThree][0]
		avgAlpha /= 4.0

		//Find sthe P(beta) for the tetrapeptide.
		avgBeta := parameters[aaIndex][1]
		avgBeta += parameters[aaIndexPlusOne][1]
		avgBeta += parameters[aaIndexPlusTwo][1]
		avgBeta += parameters[aaIndexPlusThree][1]
		avgBeta /= 4.0

		//Creates a boolean to see if it satisfies the inequality.
		satisfiesCond3 := avgTurn > avgAlpha && avgTurn > avgBeta

		//Checks to see if all three conditions are met.
		if satisfiesCond1 && satisfiesCond2 && satisfiesCond3 {
			//Creates a new turn data type.
			//Sets the index of the turn.
			newTurn := new(Turn)
			newTurn.Index = i

			//Appends the new turn to the slice of turns.
			turns = append(turns, *newTurn)
		}

	}

	//Now that we have created a list of
	//Checks to see if the given index is a beta sheet ->
	return MergeTurns(turns, secStruc)
}

func MergeTurns(identifiedTurns []Turn, secStruc []ABHelixSheet) []ABHelixSheet {

	for i := range identifiedTurns {
		flagBreak := false
		indTurn := identifiedTurns[i].Index
	out:
		for j, structure := range secStruc {
			if indTurn > structure.EndIndex {
				continue
			}
			if structure.typeAB != "sheet" {
				continue
			}

			newTurn := new(ABHelixSheet)
			newTurn.StartIndex = indTurn
			newTurn.EndIndex = indTurn
			newTurn.typeAB = "Turn"

			if inSecondStructure(indTurn, structure) {
				//Three cases
				//index of the turn is at start index
				if indTurn == structure.StartIndex {
					fmt.Println("Start Index")
					if j == 0 {
						secStruc[j].StartIndex += 1
						secStruc = append([]ABHelixSheet{*newTurn}, secStruc...)
					} else {
						upToTurn := make([]ABHelixSheet, len(secStruc[:j]))
						afterTurn := make([]ABHelixSheet, len(secStruc[j:]))
						upToTurn = CopySecondaryStructure(secStruc[:j])
						upToTurn[len(upToTurn)-1].EndIndex -= 1
						afterTurn = CopySecondaryStructure(secStruc[j:])
						afterTurn[0].StartIndex += 1
						secStruc = append(upToTurn, *newTurn)
						secStruc = append(secStruc, afterTurn...)
					}
				} else if indTurn == structure.EndIndex {
					fmt.Println("End index", j)
					upToTurn := make([]ABHelixSheet, 0)
					afterTurn := make([]ABHelixSheet, 0)
					if j != 0 {

						upToTurn = CopySecondaryStructure(secStruc[:j+1])
						fmt.Println(secStruc[:j])
						upToTurn[j].EndIndex -= 1
					} else {
						upToTurn = []ABHelixSheet{secStruc[j]}
						upToTurn[j].EndIndex -= 1
					}
					afterTurn = CopySecondaryStructure(secStruc[j+1:])

					secStruc = append(upToTurn, *newTurn)
					secStruc = append(secStruc, afterTurn...)
				} else {
					//split the secondary structure into two
					fmt.Println("Split")
					firstStruc := *new(ABHelixSheet)
					secondStruc := *new(ABHelixSheet)

					firstStruc = CopyHelixSheet(secStruc[j])
					secondStruc = CopyHelixSheet(firstStruc)

					firstStruc.EndIndex = indTurn - 1
					secondStruc.StartIndex = indTurn + 1

					if secondStruc.StartIndex > secondStruc.EndIndex {
						panic("This shouldn't happen! Start index is bigger!")
					}

					upToTurn := make([]ABHelixSheet, len(secStruc[:j]))
					afterTurn := make([]ABHelixSheet, len(secStruc[j:]))
					upToTurn = CopySecondaryStructure(secStruc[:j])
					afterTurn = CopySecondaryStructure(secStruc[j+1:])

					secStruc = append(upToTurn, firstStruc)
					secStruc = append(secStruc, *newTurn)
					secStruc = append(secStruc, secondStruc)
					secStruc = append(secStruc, afterTurn...)

				}
				flagBreak = true
				fmt.Println("exit turn")
			}
			if flagBreak {
				break out
			}
		}
		fmt.Println(secStruc)
	}

	return secStruc

}

func CopySecondaryStructure(secStruc []ABHelixSheet) []ABHelixSheet {
	newSecStruc := make([]ABHelixSheet, len(secStruc))
	for i := range secStruc {
		newSecStruc[i] = CopyHelixSheet(secStruc[i])
	}

	return newSecStruc
}

func CopyHelixSheet(helixOrSheet ABHelixSheet) ABHelixSheet {
	newHelixSheet := new(ABHelixSheet)
	newHelixSheet.typeAB = helixOrSheet.typeAB
	newHelixSheet.Score = helixOrSheet.Score
	newHelixSheet.StartIndex = helixOrSheet.StartIndex
	newHelixSheet.EndIndex = helixOrSheet.EndIndex

	return *newHelixSheet
}

func inSecondStructure(index int, betaSheet ABHelixSheet) bool {
	startInd, endInd := betaSheet.StartIndex, betaSheet.EndIndex
	if startInd <= index && index <= endInd {
		return true
	}
	return false
}

// Function to assess accuracy of the program.
func assessAccuracy(predSS, realSS string) (float64, float64, float64, float64) {
	if len(predSS) != len(realSS) {
		print(len(predSS), len(realSS))
		panic("Not the same length secondary structure!")
	}

	var count float64
	var countCorrectHelices, totalHelices float64
	var countCorrectSheets, totalSheets float64
	var countCorrectCoils, totalCoils float64

	for i := range predSS {
		fmt.Println(predSS[i], realSS[i])
		if predSS[i] == realSS[i] {
			count++
			if string(predSS[i]) == "H" {
				countCorrectHelices++
			} else if string(predSS[i]) == "S" {
				countCorrectSheets++
			} else {
				countCorrectCoils++
			}
		}

		if string(realSS[i]) == "H" {
			totalHelices++
		} else if string(realSS[i]) == "S" {
			totalSheets++
		} else {
			totalCoils++
		}

	}

	acc := count / float64(len(predSS))

	accH := countCorrectHelices / totalHelices
	accS := countCorrectSheets / totalSheets
	accC := countCorrectCoils / totalCoils

	return acc, accH, accS, accC
}

func convertAASecStrucToString(aaSecStruc []int) string {
	totalString := ""
	for _, val := range aaSecStruc {
		if val == 2 {
			totalString += "S"
		} else if val == 1 {
			totalString += "H"
		} else {
			totalString += "L"
		}
	}

	return totalString
}

func TestGorAndFasman(dir string, parameters [][]float64, aaIndexMap map[rune]int) {
	//Range through all of the results validation tests.
	files := ReadDirectory(dir)
	ProteinPredArray := make([][]CFScore, len(files))
	ProteinPredGorArray := make([][]CFScore, len(files))

	alphaParam := GORMethodInput("GorParams/alpha_helix.txt")
	sheetParam := GORMethodInput("GorParams/beta_strand.txt")
	turnParam := GORMethodInput("GorParams/beta_turn.txt")
	coilParam := GORMethodInput("GorParams/coil.txt")

	file, err := os.Create("results/Output/results.csv")
	if err != nil {
		panic(err)
	}
	defer file.Close()
	fmt.Fprintf(file, "Protein Name, accCF, accCFHelices, accCFSheets, accCFCoils, accGor, accGorHelices, accGorSheets, accGorCoils")
	for i, inputFile := range files {
		//Do Chou-Fasman on this secondary structure
		//create a protein from this input file
		var newProt Protein
		var actualStructure string
		newProt.Sequence, actualStructure = GORReadTests(dir + "/" + inputFile.Name())
		newProt.Identifier = inputFile.Name()[0:5]
		actualStructure = mapSSToLocalStructure(actualStructure)
		ProteinPredArray[i] = ChouFasman(newProt, parameters, aaIndexMap)
		//Predict helices
		helices := IdentifyHelicies(ProteinPredArray[i])

		// predict beta sheets
		betaSheets := IdentifySheets2(ProteinPredArray[i])

		// reassign the helices and sheets as appropriate
		reassignedABHelixSheet := AHelicalBSheetAssignment(append(helices, betaSheets...))
		reassignedABHelixSheet = FillGapsInSequence(len(ProteinPredArray), reassignedABHelixSheet)

		//Reassign structures for turns
		//reassignedABHelixSheet = IdentifyTurns(newProt, parameters, aaIndexMap, reassignedABHelixSheet)

		aaSecStructCFSingle := make([]int, len(newProt.Sequence))

		// make int slice, which would store the information about the secondary structure of each position
		// (1 = helix, 2 = sheet, 3 = loop).

		// It should be initialized with 3s (as if not helix or sheet, then loop) and we assign helix and sheet below
		for i := 0; i < len(aaSecStructCFSingle); i++ {
			aaSecStructCFSingle[i] = 3
		}

		// convert the reassignedABHelixSheet slice to aaSecStruct slice
		aaSecStructSingle := ConvertABHelixSheetToAASecStruct(reassignedABHelixSheet, aaSecStructCFSingle)

		cfstr := convertAASecStrucToString(aaSecStructSingle)
		accCF, accCFHelices, accCFSheets, accCFCoils := assessAccuracy(cfstr, actualStructure)

		//Do GOR on this file
		params := make([][][]float64, 4)
		params[0], params[1], params[2], params[3] = alphaParam, sheetParam, turnParam, coilParam
		aaIndexMap := ReadAAIndexMap("aa_index_map_gor.txt")
		ProteinPredGorArray[i] = GORPrediction(newProt, params, aaIndexMap)
		ssStruc := GorPredictionConv(ConvertPredToArr(ProteinPredGorArray[i]))
		// make int slice, which would store the information about the secondary structure of each position
		// (1 = helix, 2 = sheet, 3 = loop).
		aaSecStructGORSingle := make([]int, len(newProt.Sequence))

		// It should be initialized with 3s (as if not helix or sheet, then loop) and we assign helix and sheet below
		for i := 0; i < len(aaSecStructGORSingle); i++ {
			aaSecStructGORSingle[i] = 3
		}

		// convert the reassignedABHelixSheet slice to aaSecStruct slice
		aaSecStructSingle = ConvertABHelixSheetToAASecStruct(ssStruc, aaSecStructGORSingle)

		gorstr := convertAASecStrucToString(aaSecStructSingle)
		accGor, accGorHelices, accGorSheets, accGorCoils := assessAccuracy(gorstr, actualStructure)

		//fmt.Println(accCF, ",", accCFHelices, ",", accCFSheets, ",", accCFCoils, ",", accGor, ",", accGorHelices, ",", accGorSheets, ",", accGorCoils)

		accCFStr := strconv.FormatFloat(accCF, 'f', -1, 64)
		accCFHelicesStr := strconv.FormatFloat(accCFHelices, 'f', -1, 64)
		accCFSheetsStr := strconv.FormatFloat(accCFSheets, 'f', -1, 64)
		accCFCoilsStr := strconv.FormatFloat(accCFCoils, 'f', -1, 64)
		accGorStr := strconv.FormatFloat(accGor, 'f', -1, 64)
		accGorHelicesStr := strconv.FormatFloat(accGorHelices, 'f', -1, 64)
		accGorSheetsStr := strconv.FormatFloat(accGorSheets, 'f', -1, 64)
		accGorCoilsStr := strconv.FormatFloat(accGorCoils, 'f', -1, 64)

		fmt.Fprintf(file, "\n")
		fmt.Fprintf(file, newProt.Identifier+","+accCFStr+","+accCFHelicesStr+","+accCFSheetsStr+","+accCFCoilsStr+","+accGorStr+","+accGorHelicesStr+","+accGorSheetsStr+","+accGorCoilsStr)
	}

}

// This functtion takes in the inputs from the GorParams folder
func GORReadTests(file string) (string, string) {
	f, err := os.Open(file)
	if err != nil {
		log.Panicf("Cannot open file: %s", err)
	}
	var ss string
	var seq string
	scanner := bufio.NewScanner(f)

	var counter int
	for scanner.Scan() {
		line := scanner.Text()
		if counter == 0 {
			seq = line
		}
		if counter == 1 {
			ss = line
		}
		counter++
	}

	defer f.Close()

	return seq, ss
}

// ReadDirectory reads in a directory and returns a slice of fs.DirEntry
// objects containing file info for the directory.
func ReadDirectory(dir string) []fs.DirEntry {
	//read in all files in the given directory
	files, err := os.ReadDir(dir)
	//Checking if there is an error.
	if err != nil {
		panic(err)
	}
	//Return the files from the directory.
	return files
}

func mapSSToLocalStructure(sequenceFromSS string) string {
	var newSequence string
	for i := range sequenceFromSS {
		if string(sequenceFromSS[i]) == "C" {
			newSequence += "L"
		} else if string(sequenceFromSS[i]) == "E" {
			newSequence += "S"
		} else {
			newSequence += "H"
		}
	}
	return newSequence
}
