package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

// function to take a window sequence and return the average chou-fasman score of the window
func ChouFasmanWindow(window string, parameters [][]float64, aaIndexMap map[rune]int) int {

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

	// CAN ALSO RETURN THE SCORES FOR EACH OF THE THREE (LATER ON)

	// Return the number corresponding to the highest score of the three
	// return helix
	if alphaHelix > betaSheet && alphaHelix > loop {
		return 1
	}
	// return sheet
	if betaSheet > alphaHelix && betaSheet > loop {
		return 2
	}
	// return loop
	return 3
}

// function to read a text file containing map of amino acids to their indices
func ReadAAIndexMap(file string) map[rune]int {

	// Define the map
	aaIndexMap := make(map[rune]int)

	// Open the file
	f, err := os.Open(file)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()

	// Read the file line by line
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		// Split the line into a slice of strings
		lineSplit := strings.Split(line, " ")

		// make the first character as key and the second as value
		// convert the first character to rune
		key := rune(lineSplit[0][0])
		// convert the second part to int
		val, _ := strconv.Atoi(lineSplit[1])

		aaIndexMap[key] = val // [0], as [1] is the error term
	}

	// Return the parameters array
	return aaIndexMap
}

// function to read a text file containing the parameters for the Chou-Fasman algorithm
func ReadParameters(file string) [][]float64 {

	// Define the parameters array
	parameters := make([][]float64, 0)

	// Open the file
	f, err := os.Open(file)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()

	// Read the file line by line
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		// Split the line into a slice of strings
		lineSplit := strings.Split(line, " ")
		// Convert the slice of strings into a slice of floats
		lineFloat := make([]float64, len(lineSplit))
		for i := 0; i < len(lineSplit); i++ {
			lineFloat[i], _ = strconv.ParseFloat(lineSplit[i], 64)
		}
		// Add the slice of floats to the parameters array
		parameters = append(parameters, lineFloat)
	}

	// Return the parameters array
	return parameters
}

// Function to call the Chou-Fasman algorithm on a window of a sequence and return the prediction array
// func ChouFasman(seq string, parameters [][]float64, aaIndexMap map[rune]int) []int {
func ChouFasman(protein Protein, parameters [][]float64, aaIndexMap map[rune]int) PredArray {

	// Define the prediction array
	//predictionArray := make([]int, len(seq)) // 0: not determined, 1: helix, 2: sheet, 3: loop
	predictionArray := make(PredArray, len(protein.Sequence))

	for i := 0; i < len(protein.Sequence); i++ {
		//predictionArray[leftPointer+windowSize/2] = ChouFasmanWindow(seq[leftPointer:rightPointer+1], parameters, aaIndexMap)
		predictionArray[i] = BetterChouFasmanWindow(protein.Sequence[i:i+1], parameters, aaIndexMap)
	}

	//fmt.Println(predictionArray)
	return predictionArray
}

// read in parameters as a dictionary - INCOMPLETE
// Input: file string
// Output: Map from amino acid character to a score given by that amino acid
func ReadParametersDict(file string) map[string]float64 {
	//Create map
	param := make(map[string]float64)
	namesToChars := NameToChar("Names.txt")
	//open file
	f, err := os.Open(file)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		// Split the line into a slice of strings
		lineSplit := strings.Split(line, " ")

		//First value is amino acid name
		aminoChar := lineSplit[0]
		aminoChar = namesToChars[aminoChar]

		//Second value is P(a).
		aminoScore, err := strconv.ParseInt(lineSplit[1], 10, 63)

		//If there is an error panic.
		if err != nil {
			panic(err)

		}

		//Third Value is P(b)

		//Fourth value is P(turn)

		//Fifth value is P(i+1)

		//Sixth value is P(i+2)

		//Seventh Value is P(i+3)

		//Map the amino acid character to the amino acid score.
		param[aminoChar] = float64(aminoScore)

	}

	// Return the parameters hashmap.
	return param
}

// Function that reads amino acid string name to a character. See Names.txt to examine file mappings from names to chars
// Input: A file name that contains mapping from amino acid name to amino acid char. (I.e. Isoleucine -> I)
// Output: A hashmap that maps a amino acid name to amino acid character.
func NameToChar(file string) map[string]string {
	namesToChars := make(map[string]string)
	//open file
	f, err := os.Open(file)
	if err != nil {
		fmt.Println(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		// Split the line into a slice of strings
		lineSplit := strings.Split(line, ",")

		//First value is amino acid character
		aminoName := lineSplit[0]

		//Third value is one letter code.
		aminoChar := lineSplit[2]

		//Map the amino acid character to the amino acid score.
		namesToChars[aminoName] = aminoChar

	}
	return namesToChars
}

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

func FillGapsInSequence(n int, sequence []ABHelixSheet) []ABHelixSheet {
	// Sort the input slice based on StartIndex
	sort.Slice(sequence, func(i, j int) bool {
		return sequence[i].StartIndex < sequence[j].StartIndex
	})

	var result []ABHelixSheet

	currentIndex := 0
	for _, element := range sequence {
		// Check for a gap before the current element
		if element.StartIndex > currentIndex {
			result = append(result, ABHelixSheet{
				StartIndex: currentIndex,
				EndIndex:   element.StartIndex - 1,
				Score:      0,
				typeAB:     "neither",
			})
		}
		result = append(result, element)
		currentIndex = element.EndIndex + 1
	}

	// Check for a gap at the end
	if currentIndex < n {
		result = append(result, ABHelixSheet{
			StartIndex: currentIndex,
			EndIndex:   n - 1,
			Score:      0,
			typeAB:     "neither",
		})
	}

	return result
}

func GenerateDNAReader(fileName string) *DNAReader {
	file, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}

	// (Shashank) below line needed to be commented out to make the file non empty. Else, the ReadProteins
	// function will return an empty slice
	// defer file.Close()
	return &DNAReader{file: file}
}

// Translates DNA read from given DNA file into a slice of Protein. Checks 3 ORFs and returns the largest one.
func TranslateDNA(dnaReader *DNAReader) Protein {
	var allDNA []string
	var potentialProteins []Protein
	scanner := bufio.NewScanner(dnaReader.file)
	for scanner.Scan() {
		line := scanner.Text()
		allDNA = append(allDNA, line)
	}
	dnaString := strings.Join(allDNA, "")
	fmt.Println("Imported DNA sequence as", dnaString)
	rnaString := TranscribeDNA(dnaString)
	fmt.Println("Transcribed DNA to RNA Sequench:", rnaString)
	id := 0
	for j := 0; j < 3; j++ {
		var workingAcids []AminoAcid
		for i := j; i < len(rnaString)-2; i += 3 {
			codon := rnaString[i : i+3]
			nextAcid := TranslateCodon(codon)
			if nextAcid.Identifier != "*" {
				workingAcids = append(workingAcids, nextAcid)
			} else {
				nextProtein := AcidsToProtein(workingAcids, id)
				id++
				potentialProteins = append(potentialProteins, nextProtein)
				workingAcids = make([]AminoAcid, 0)
			}
		}
	}
	protein := potentialProteins[0]
	n := len(potentialProteins[0].Sequence)
	for _, newProtein := range potentialProteins {
		if len(newProtein.Sequence) > n {
			protein = newProtein
			n = len(newProtein.Sequence)
		}
	}
	fmt.Println("Translated sequence with longest ORF:", protein)
	return protein
}

// TranscribeDNA transcribes a DNA string into RNA
func TranscribeDNA(dnaString string) string {
	rnaString := strings.Replace(dnaString, "T", "U", -1)
	rnaString = strings.Replace(dnaString, "t", "u", -1)
	return rnaString
}

// TranslateCodon translates a given codon into its corresponding protein
func TranslateCodon(codon string) AminoAcid {
	codon = strings.ToUpper(codon)
	codonToAminoAcid := map[string]string{
		"UUU": "F", "UUC": "F",
		"UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
		"AUU": "I", "AUC": "I", "AUA": "I",
		"AUG": "M",
		"GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
		"UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
		"CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
		"ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
		"GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
		"UAU": "Y", "UAC": "Y",
		"UAA": "*", "UAG": "*", "UGA": "*",
		"CAU": "H", "CAC": "H",
		"CAA": "Q", "CAG": "Q",
		"AAU": "N", "AAC": "N",
		"AAA": "K", "AAG": "K",
		"GAU": "D", "GAC": "D",
		"GAA": "E", "GAG": "E",
		"UGU": "C", "UGC": "C",
		"UGG": "W",
		"CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
		"AGU": "S", "AGC": "S",
		"GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
	}

	aminoAcid, exists := codonToAminoAcid[codon]
	if !exists {
		return AminoAcid{"unknown"}
	}

	return AminoAcid{Identifier: aminoAcid}
}

// AcidsToProtein takes in an array of amino acides and joins them into a protein object.
func AcidsToProtein(acids []AminoAcid, id int) Protein {
	var acidStrings []string
	for _, acid := range acids {
		acidStrings = append(acidStrings, acid.Identifier)
	}
	aaString := strings.Join(acidStrings, "")
	return Protein{Identifier: fmt.Sprint(id), Sequence: aaString}
}
