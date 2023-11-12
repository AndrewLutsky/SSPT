package main

// BetterChouFasmanWindow takes a window sequence and returns a CFScore of each of the parameters
func BetterChouFasmanWindow(window string, parameters [][]float64, aaIndexMap map[rune]int) CFScore {

	alphaHelix := float64(0)
	betaSheet := float64(0)
	loop := float64(0)

	// Calculate the score for the window
	for i := 0; i < len(window); i++ {
		aaIndex := aaIndexMap[rune(window[i])]
		alphaHelix += parameters[aaIndex][0]
		betaSheet += parameters[aaIndex][1]
		loop += parameters[aaIndex][2]
	}

	var score CFScore
	score.Helix = alphaHelix
	score.Sheet = betaSheet
	score.Loop = loop

	return score
}

// IdentifyHelicies is a function that takes as its input a prediction array of CFScores and
// returns an array of Helix objects identified from the sequence
func IdentifyHelicies(predArray PredArray) []ABHelixSheet {
	helicies := make([]ABHelixSheet, 0)
	n := len(predArray)
	// scan through array until viable starting points are found
	for i := 0; i < n-5; i++ { // can't range here, have to manually delcare so that we can enforce window skipping
		startPoint := i
		endPoint := i + 5
		predStartArray := predArray[startPoint:endPoint]
		validHStart, _ := ValidStartingPosition(predStartArray, 4)
		//numHits, _ := NumHits(predStartArray, 4)

		if validHStart { // Valid starting array
			//fmt.Println("Found the starting array for a helix! NumHits:", numHits, "Start Point:", startPoint, "End Point:", endPoint, ". Extending backwards.")
			// Start extending backwards
			backFlag := true
			j := startPoint
			for backFlag && j > 0 {
				j -= 1
				backArray := predArray[j : j+3]
				backArrayHScore, _ := ArrayAverage(backArray)
				backFlag = backArrayHScore >= 1.00
				//fmt.Println("Searching backwards... j=", j, "backArrayHScore:", backArrayHScore)
			}
			startPoint = j
			//fmt.Println("Found backmost point at", startPoint, "! Searching forwards.")
			// Start extending forwards
			forFlag := true
			k := endPoint
			for forFlag && k < n-3 {
				k += 1
				forArray := predArray[k : k+3]
				forArrayHScore, _ := ArrayAverage(forArray)
				forFlag = forArrayHScore >= 1.00
				//fmt.Println("Searching forwards... k=", k, "forArrayHScore:", forArrayHScore)
			}
			endPoint = k
			//fmt.Println("Found foremost point at", endPoint, "!")
		}

		avgHelixScore, avgSheetScore := ArrayAverage(predArray[startPoint:endPoint])
		if avgHelixScore > avgSheetScore { // We have a valid helix! Create the helix and add it to the array, then skip to the end of the helix.
			var newHelix ABHelixSheet
			newHelix.StartIndex = startPoint
			newHelix.EndIndex = endPoint
			newHelix.Score = avgHelixScore
			newHelix.typeAB = "helix"
			helicies = append(helicies, newHelix)
		}
		i = endPoint
	}
	return helicies
}

// IdentifySheets2 is a function that takes as its input a prediction array of CFScores and
// returns an array of Sheet objects identified from the sequence
func IdentifySheets2(predArray PredArray) []ABHelixSheet {
	sheets := make([]ABHelixSheet, 0)
	n := len(predArray)
	// scan through array until viable starting points are found
	for i := 0; i < n-4; i++ { // can't range here, have to manually delcare so that we can enforce window skipping
		startPoint := i
		endPoint := i + 4
		predStartArray := predArray[startPoint:endPoint]
		_, validSStart := ValidStartingPosition(predStartArray, 3)

		if validSStart { // Valid starting array
			// Start extending backwards
			backFlag := true
			j := startPoint
			for backFlag && j > 0 {
				j -= 1
				backArray := predArray[j : j+3]
				_, backArraySScore := ArrayAverage(backArray)
				backFlag = backArraySScore >= 1.00
			}
			startPoint = j
			// Start extending forwards
			forFlag := true
			k := endPoint
			for forFlag && k < n-3 {
				k += 1
				forArray := predArray[k : k+3]
				_, forArraySScore := ArrayAverage(forArray)
				forFlag = forArraySScore >= 1.00
			}
			endPoint = k
		}

		avgHelixScore, avgSheetScore := ArrayAverage(predArray[startPoint:endPoint])
		if avgHelixScore > avgSheetScore && avgSheetScore > 1.05 { // We have a valid helix! Create the helix and add it to the array, then skip to the end of the helix.
			var newSheet ABHelixSheet
			newSheet.StartIndex = startPoint
			newSheet.EndIndex = endPoint
			newSheet.Score = avgSheetScore
			newSheet.typeAB = "sheet"
			sheets = append(sheets, newSheet)
		}
		i = endPoint
	}
	return sheets
}

// ArrayAverage takes in a PredArray subslice and returns the average helix and sheet scores, respectively.
func ArrayAverage(a PredArray) (float64, float64) {
	n := float64(len(a))
	hscore := 0.0
	sscore := 0.0
	for _, item := range a {
		hscore += item.Helix
		sscore += item.Sheet
	}
	return hscore / n, sscore / n
}

// ValidStartingPosition takes in a PredArray subslice and an int of necessary hits numHits, and returns bools of if the slice has valid helix or sheet.
func ValidStartingPosition(a PredArray, numHits int) (bool, bool) {
	helixHits := 0
	sheetHits := 0
	for _, item := range a {
		if item.Helix > 1.0 {
			helixHits++
		}
		if item.Sheet > 1.0 {
			sheetHits++
		}
	}
	return helixHits >= numHits, sheetHits >= numHits
}

// ValidStartingPosition takes in a PredArray subslice and an int of necessary hits numHits, and returns bools of if the slice has valid helix or sheet.
func NumHits(a PredArray, numHits int) (int, int) {
	helixHits := 0
	sheetHits := 0
	for _, item := range a {
		if item.Helix > 1.0 {
			helixHits++
		}
		if item.Sheet > 1.0 {
			sheetHits++
		}
	}
	return helixHits, sheetHits
}
