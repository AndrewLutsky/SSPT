package main

//Function to identify bends in the protein.

/*
To identify a bend at residue number j, calculate the following value

p(t) = f(j)f(j+1)f(j+2)f(j+3) where the f(j+1) value for the j+1 residue is used,
the f(j+2) value for the j+2 residue is used and the f(j+3) value for the j+3 residue is used.
If: (1) p(t) > 0.000075; (2) the average value for P(turn) > 1.00 in the tetrapeptide;
and (3) the averages for the tetrapeptide obey the inequality P(a-helix) < P(turn) > P(b-sheet),
then a beta-turn is predicted at that location.
*/
func IdentifyTurns(protein Protein, parameters [][]float64, aaIndexMap map[rune]int, secStruc []ABHelixSheet) []ABHelixSheet {
	turns := make([]Turn, 0)
	for i := 0; i < len(protein.Sequence)-3; i++ {
		//Step 1 calculate the product of f(j) * f(j+1) * f(j+2) * f(j+3)
		aaIndex := aaIndexMap[rune(protein.Sequence[i])]
		aaIndexPlusOne := aaIndexMap[rune(protein.Sequence[i+1])]
		aaIndexPlusTwo := aaIndexMap[rune(protein.Sequence[i+2])]
		aaIndexPlusThree := aaIndexMap[rune(protein.Sequence[i+3])]
		p_t := parameters[aaIndex][3] * parameters[aaIndexPlusOne][4] * parameters[aaIndexPlusTwo][5] * parameters[aaIndexPlusThree][6]
		satisfiesCond1 := p_t > 0.000075

		//Step 2 calculate the average of P(turn) for the tetrapeptide.
		//sum all the P(turns)
		avgTurn := parameters[aaIndex][2]
		avgTurn += parameters[aaIndexPlusOne][2]
		avgTurn += parameters[aaIndexPlusTwo][2]
		avgTurn += parameters[aaIndexPlusThree][2]

		satisfiesCond2 := avgTurn/4.0 > 1.0

		//Step 3 check inequality

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

		//Checks the inequality.
		satisfiesCond3 := avgTurn > avgAlpha && avgTurn > avgBeta

		//Checks to see if all three conditions are met.
		if satisfiesCond1 && satisfiesCond2 && satisfiesCond3 {
			newTurn := new(Turn)
			newTurn.Index = i
			turns = append(turns, *newTurn)
		}

	}

	//Checks to see if the given index is a beta sheet ->
	for i := range turns {
		indTurn := turns[i].Index
		for j, structure := range secStruc {
			if indTurn > structure.EndIndex {
				continue
			}
			newTurn := new(ABHelixSheet)
			newTurn.StartIndex = indTurn
			newTurn.EndIndex = indTurn
			newTurn.typeAB = "Turn"
			newSecStruc := make([]ABHelixSheet, len(secStruc))

			if inSecondStructure(indTurn, structure) {
				//Three cases
				//index of the turn is at start index
				if indTurn == structure.StartIndex {
					newSecStruc = append(secStruc[:j-1], *newTurn)
					newSecStruc = append(newSecStruc, secStruc[j:]...)
					newSecStruc[j].StartIndex += 1
					if j != 0 {
						newSecStruc[j-1].EndIndex -= 1
					}
					return newSecStruc
				} else if indTurn == structure.EndIndex {
					newSecStruc = append(secStruc[:j], *newTurn)
					newSecStruc = append(newSecStruc, secStruc[j:]...)
					newSecStruc[j].EndIndex -= 1
					if !(j >= len(secStruc)) {
						newSecStruc[j+1].StartIndex += 1
					}
				} else {
					//split the secondary structure into two
					firstStruc := new(ABHelixSheet)
					secondStruc := new(ABHelixSheet)

					firstStruc = &secStruc[j]
					secondStruc = firstStruc

					firstStruc.EndIndex = indTurn - 1
					secondStruc.StartIndex = indTurn + 1

					if secondStruc.StartIndex > secondStruc.EndIndex {
						panic("This shouldn't happen! Start index is bigger!")
					}

					newSecStruc = append(secStruc[:j-1], *firstStruc)
					newSecStruc = append(newSecStruc, *newTurn)
					newSecStruc = append(newSecStruc, *secondStruc)
					newSecStruc = append(newSecStruc, secStruc[j+1:]...)
					return newSecStruc

				}
			}
		}
	}

	return secStruc
}

func inSecondStructure(index int, betaSheet ABHelixSheet) bool {
	startInd, endInd := betaSheet.StartIndex, betaSheet.EndIndex
	if startInd <= index && index <= endInd {
		return true
	}
	return false
}
