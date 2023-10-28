package main

import (
	"bufio"
	"fmt"
	"os"
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
// func ChouFasman(seq string, windowSize int, parameters [][]float64, aaIndexMap map[rune]int) []int {
func ChouFasman(protein Protein, windowSize int, parameters [][]float64, aaIndexMap map[rune]int) PredArray {

	// Define the prediction array
	//predictionArray := make([]int, len(seq)) // 0: not determined, 1: helix, 2: sheet, 3: loop
	predictionArray := make(PredArray, len(protein.Sequence))

	if len(protein.Sequence) < windowSize+1 {
		fmt.Println("Sequence too short")
		return predictionArray
	}

	leftPointer := 0
	rightPointer := windowSize - 1

	for rightPointer < len(protein.Sequence) {
		//predictionArray[leftPointer+windowSize/2] = ChouFasmanWindow(seq[leftPointer:rightPointer+1], parameters, aaIndexMap)
		predictionArray[leftPointer+windowSize/2] = BetterChouFasmanWindow(protein.Sequence[leftPointer:rightPointer+1], parameters, aaIndexMap)
		leftPointer += 1
		rightPointer += 1
	}

	fmt.Println(predictionArray)
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
