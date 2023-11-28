package main

// The GOR (Garnier-Osguthorpe-Robson) method for predicting protein secondary structure is based on information theory
// and uses a statistical approach.
func GORPrediction(protein Protein, parameters [][]float64, aaIndexMap map[rune]int) PredArray {
	windowSize := 17 // Typical window size used in GOR
	halfWindowSize := windowSize / 2
	sequence := protein.Sequence
	prediction := make(PredArray, len(sequence))

	// Loop through the protein sequence
	for i := 0; i < len(sequence); i++ {
		var scoreHelix, scoreStrand, scoreCoil float64

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
		//Probably needs to be fixed
		for j := windowStart; j <= windowEnd; j++ {
			aaIndex := aaIndexMap[rune(sequence[j])]
			scoreHelix += parameters[aaIndex][0]
			scoreStrand += parameters[aaIndex][1]
			scoreCoil += parameters[aaIndex][2]
		}

		// Predict the structure based on the highest score
		if scoreHelix > scoreStrand && scoreHelix > scoreCoil {
			prediction[i] = CFScore{Helix: 1} // Helix
		} else if scoreStrand > scoreHelix && scoreStrand > scoreCoil {
			prediction[i] = CFScore{Sheet: 1} // Sheet
		} else {
			prediction[i] = CFScore{Loop: 1} // Loop
		}
	}

	return prediction
}

func GorPredictionConv(predArr PredArray) []ABHelixSheet {
	ssStruc := make([]ABHelixSheet, 0)
	//intermediateArr := ConvertPredToArr(predArr)

	return ssStruc
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
