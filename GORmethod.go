package main

// The GOR (Garnier-Osguthorpe-Robson) method for predicting protein secondary structure is based on information theory
// and uses a statistical approach.
func GORPrediction(protein Protein, parameters [][][]float64, aaIndexMap map[rune]int) PredArray {
	windowSize := 17 // Typical window size used in GOR
	halfWindowSize := windowSize / 2
	sequence := protein.Sequence
	prediction := make(PredArray, len(sequence))

	// Loop through the protein sequence
	for i := 0; i < len(sequence); i++ {
		var scoreHelix, scoreStrand, scoreTurn, scoreCoil float64
		var diff int
		//Need to reimplement
		//at 0 it window starts at zero
		//at 8 it starts at 0
		// Calculate the start and end of the sliding window
		windowStart := i - halfWindowSize
		if windowStart < 0 {
			windowStart = 0
		}
		windowEnd := i + halfWindowSize
		if windowEnd >= len(sequence) {
			windowEnd = len(sequence) - 1
		}

		// Loop through the window
		// at i = 0 [0:8]
		// at i = 1 [0:9]
		// at i = 2 [0:10]
		//at i = 3 [0:11]
		//... until [0:17]
		//at i = 9 [1:18]
		if i >= 9 {
			diff = i - 8
		}
		for j := windowStart; j <= windowEnd; j++ {

			aaIndex := aaIndexMap[rune(sequence[j])]
			scoreHelix += parameters[0][aaIndex][j-diff]
			scoreStrand += parameters[1][aaIndex][j-diff]
			scoreTurn += parameters[2][aaIndex][j-diff]
			scoreCoil += parameters[3][aaIndex][j-diff]
		}

		// Predict the structure based on the highest score
		if scoreHelix > scoreStrand && scoreHelix > scoreCoil && scoreHelix > scoreCoil {
			prediction[i] = CFScore{Helix: scoreHelix} // Helix
		} else if scoreStrand > scoreHelix && scoreStrand > scoreCoil && scoreStrand > scoreTurn {
			prediction[i] = CFScore{Sheet: scoreStrand} // Sheet
		} else if scoreTurn > scoreHelix && scoreTurn > scoreStrand && scoreTurn > scoreTurn {
			prediction[i] = CFScore{Turn: scoreTurn} // Loop
		} else {
			prediction[i] = CFScore{Loop: scoreCoil}
		}
	}

	return prediction
}

func GorPredictionConv(intermediateArray []int) []ABHelixSheet {
	ssStruc := make([]ABHelixSheet, 0)
	//intermediateArray := ConvertPredToArr(predArr)
	lastNum := intermediateArray[0]
	currentSecStruc := new(ABHelixSheet)
	currentSecStruc.typeAB = ConvertIntToType(intermediateArray[0])
	for i := 0; i < len(intermediateArray); i++ {
		currentNum := intermediateArray[i]

		//make new ABHelixSheet
		if currentNum != lastNum {
			//append to secondary structure slice
			ssStruc = append(ssStruc, *currentSecStruc)
			newSecStruc := new(ABHelixSheet)
			newSecStruc.StartIndex = i
			newSecStruc.EndIndex = i
			newSecStruc.typeAB = ConvertIntToType(intermediateArray[i])
			currentSecStruc = newSecStruc
		} else {
			//increment the end index
			currentSecStruc.EndIndex = i
		}

		lastNum = intermediateArray[i]
	}

	return ssStruc
}

func ConvertIntToType(val int) string {
	if val == 0 {
		return "helix"
	} else if val == 1 {
		return "sheet"
	} else if val == 2 {
		return "loop"
	} else {
		return "coil"
	}
}

func ConvertPredToArr(predArr PredArray) []int {
	ssStruc := make([]int, len(predArr))
	for i := range predArr {
		ssValues := make([]float64, 4)
		ssValues[0] = predArr[i].Helix
		ssValues[1] = predArr[i].Sheet
		ssValues[2] = predArr[i].Loop
		ind := max(ssValues)
		ssStruc[i] = ind
	}

	return ssStruc
}

func max(vals []float64) int {
	max := vals[0]
	index := 0
	for i, val := range vals {
		if val > max {
			max = val
			index = i
		}
	}

	return index
}
