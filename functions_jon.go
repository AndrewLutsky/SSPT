package main

// BetterChouFasmanWindow takes a window sequence and returns a CFScore of each of the parameters
func BetterChouFasmanWindow(window string, parameters [][]float64, aaIndexMap map[rune]int) CFScore {

	alphaHelix := float64(0)
	betaSheet := float64(0)
	loop := float64(0)
	//Needed toa dd new variables for ith residue, i+1, i+2...
	f_i := float64(0)
	f_i1 := float64(0)
	f_i2 := float64(0)
	f_i3 := float64(0)

	// Calculate the score for the window
	for i := 0; i < len(window); i++ {
		aaIndex := aaIndexMap[rune(window[i])]
		alphaHelix += parameters[aaIndex][0]
		betaSheet += parameters[aaIndex][1]
		loop += parameters[aaIndex][2]
		f_i += parameters[aaIndex][3]
		f_i1 += parameters[aaIndex][4]
		f_i2 += parameters[aaIndex][5]
		f_i3 += parameters[aaIndex][6]
	}

	var score CFScore
	score.Helix = alphaHelix
	score.Sheet = betaSheet
	score.Loop = loop
	score.F_i = f_i
	score.F_i1 = f_i1
	score.F_i2 = f_i2
	score.F_i3 = f_i3

	return score
}

// IdentifyHelicies is a function that takes as its input a prediction array of CFScores and
// returns an array of Helix objects identified from the sequence
func IdentifyHelicies(predArray PredArray) []Helix {
	helicies := make([]Helix, 0)
	n := len(predArray)
	// scan through array until viable starting points are found
	for i := 0; i < n; i++ { // can't range here, have to manually delcare so that we can enforce window skipping
		predPoint := predArray[i]
		if predPoint.Helix > 1.00 && i < n-4 { //Valid starting position
			j := i
			helixSum := predPoint.Helix // Start tracking for average
			sheetSum := predPoint.Sheet
			for j = i + 1; j < len(predArray); j++ { //Scan forward until stopping point is found or end of sequence is reached.
				helixSum += predArray[j].Helix
				sheetSum += predArray[j].Sheet
				if predPoint.Helix < 1.00 { // End of the potential helix
					break
				}
			}
			// At this point, we've either reached the end of the potential helix or hit the end of the sequence.
			// Now, validate the potential helix.
			length := j - i
			if length > 5 { // long enough to be a helix, calculate avg scores
				avgHelixScore := helixSum / float64(length)
				avgSheetScore := sheetSum / float64(length)
				if avgHelixScore > avgSheetScore { // We have a valid helix! Create the helix and add it to the array, then skip to the end of the helix.
					var newHelix Helix
					newHelix.StartIndex = i
					newHelix.EndIndex = j
					newHelix.Score = avgHelixScore
					helicies = append(helicies, newHelix)
				}
			}
			i = j // ALWAYS jump the window once you have found a helix. Chou-Fasman says that you must also extend the search window backwards when calculating helicies
			// so by always jumping the window we never have to consider this.
		}
	}
	return helicies
}
