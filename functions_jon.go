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
					var newHelix ABHelixSheet
					newHelix.StartIndex = i
					newHelix.EndIndex = j
					newHelix.Score = avgHelixScore
					newHelix.typeAB = "helix"
					helicies = append(helicies, newHelix)
				}
			}
			i = j // ALWAYS jump the window once you have found a helix. Chou-Fasman says that you must also extend the search window backwards when calculating helicies
			// so by always jumping the window we never have to consider this.
		}
	}
	return helicies
}

// func IdentifyHelicies(predArray PredArray) []ABHelixSheet {
// 	helicies := make([]ABHelixSheet, 0)
// 	n := len(predArray)

// 	averagePHelix := func(scores []CFScore) float64 {
// 		sum := 0.0
// 		for _, score := range scores {
// 			sum += score.Helix
// 		}
// 		return sum / float64(len(scores))
// 	}

// 	averagePSheet := func(scores []CFScore) float64 {
// 		sum := 0.0
// 		for _, score := range scores {
// 			sum += score.Sheet
// 		}
// 		return sum / float64(len(scores))
// 	}

// 	for i := 0; i <= n-6; i++ {
// 		// Look for 4 out of 6 residues in the window having P(a-helix) > 1.00
// 		count := 0
// 		for j := i; j < i+6; j++ {
// 			if predArray[j].Helix > 1.00 {
// 				count++
// 			}
// 		}

// 		if count >= 4 {
// 			// Potential helix start found, now extend it
// 			start := i
// 			end := i + 5

// 			// Extend forward
// 			for {
// 				if end+4 >= n { // Reached the end of array, stop extending
// 					break
// 				}
// 				// Look at the next four residues to decide whether to extend
// 				if averagePHelix(predArray[end+1:end+5]) < 1.00 {
// 					break
// 				}
// 				end++
// 			}

// 			// Extend backward
// 			for {
// 				if start-1 < 0 { // Reached the beginning of array, stop extending
// 					break
// 				}
// 				// Look at the previous four residues to decide whether to extend
// 				if averagePHelix(predArray[start-4:start]) < 1.00 {
// 					break
// 				}
// 				start--
// 			}

// 			// Calculate averages for the helix
// 			avgHelixScore := averagePHelix(predArray[start : end+1])
// 			avgSheetScore := averagePSheet(predArray[start : end+1])

// 			// Check final criteria for assigning as helix
// 			if (end-start) >= 5 && avgHelixScore > avgSheetScore {
// 				helicies = append(helicies, ABHelixSheet{
// 					StartIndex: start,
// 					EndIndex:   end,
// 					Score:      avgHelixScore,
// 					typeAB:     "helix",
// 				})
// 				i = end // Jump the window to the end of the helix
// 			}
// 		}
// 	}

// 	return helicies
// }
