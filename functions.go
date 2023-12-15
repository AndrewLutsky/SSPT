package main

import (
	"bufio"
	"fmt"
	"io"
	"io/fs"
	"log"
	"math"
	"net/http"
	"os"
	"os/exec"
	"sort"
	"strconv"
	"strings"

	"github.com/fogleman/gg"
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

// ____________FIRST 3D VISUALIZATION ATTEMPT FUNCTIONS BELOW_______________

// DISCLAIMER: I have assumed that there are exactly 4 amino acids in a helix turn and hence, each on lies on the
// origin-XZ plane or the origin-YZ plane

func Visualize(predSeq []int, pitch, aaDist float64) []Coordinate {

	// REMEMBER, OUR PROTEIN IS ALONG THE Z AXIS (AND OUR COORDINATE SYSTEM IS Y-UP)
	// an empty slice of coordinate objects, which will store the coordinates of each amino acid.
	coordinates := make([]Coordinate, len(predSeq))
	coordinates[0] = Coordinate{0, 0, 0}
	// looping through the prediction sequence, and finding the coordinates of each amino acid
	for char := 1; char < len(predSeq); char++ {
		if predSeq[char] == 1 {
			// helix
			coordinates[char] = HelixCoordinate(predSeq, pitch, aaDist, char, coordinates[char-1])
		} else if predSeq[char] == 2 {
			// sheet
			// temporally using the LoopCoordinates function instead of SheetCoordinates
			coordinates[char] = LoopCoordinates(predSeq, pitch, aaDist, char, coordinates[char-1])
		} else if predSeq[char] == 3 {
			// loop
			coordinates[char] = LoopCoordinates(predSeq, pitch, aaDist, char, coordinates[char-1])
		} else {
			panic("bad prediction (not helix/sheet/loop)")
		}
	}
	return coordinates
}

// function to find the coordinates of an amino acid, which has been predicted to be a part of a sheet
// we need to first check the number of bends in the sheet sequence and then decide the coordinates
// func SheetCoordinates(predSeq []int, pitch, aaDist float64, char int, prevCoor Coordinate) Coordinate {

// 	var coor Coordinate
// 	var newSheet bool
// 	if predSeq[char-1] != 2 {
// 		newSheet = true
// 	} else {
// 		newSheet = false
// 	}

// 	return coor
// }

// function to find the coordinates of an amino acid, which has been predicted to be a part of a loop
func LoopCoordinates(predSeq []int, pitch, aaDist float64, char int, prevCoor Coordinate) Coordinate {

	var coor Coordinate

	coor.X = 0
	coor.Y = 0
	coor.Z = prevCoor.Z + aaDist // incrementing by not the pitch but the distance between two amino acids

	return coor
}

// this is a function to find out the coordinate of an amino acid, which has been predicted to be a part
// of an alpha helix. it assumes that there are exactly 4 amino acids per turn in the helix
// and it follows left hand corkscrew rule.
func HelixCoordinate(predSeq []int, pitch, aaDist float64, char int, prevCoor Coordinate) Coordinate {

	// This assumes that the previous amino acid was either from a loop (hence on z axis)
	// or it was a part of alpha helix
	var coor Coordinate
	// location of the previous aa. Basically on which axis it lies.
	prevLoc := FindXYQuadrant(prevCoor)
	if prevLoc == 0 {
		// previous aa was on the z axis, i.e. loop
		coor.X = math.Sqrt(aaDist*aaDist - pitch*pitch) // this is how far off the z axis the aa will be
		coor.Y = 0
		coor.Z = prevCoor.Z + aaDist // incrementing by the distance because otherwise, PYMOL will make unncessary bond
	} else if prevLoc == 1 {
		// previous aa was a part of alpha helix
		coor.X = 0
		coor.Y = math.Sqrt(aaDist*aaDist - pitch*pitch)
		coor.Z = prevCoor.Z + pitch
	} else if prevLoc == 2 {
		// previous aa was a part of alpha helix
		coor.X = -math.Sqrt(aaDist*aaDist - pitch*pitch)
		coor.Y = 0
		coor.Z = prevCoor.Z + pitch
	} else if prevLoc == 3 {
		// previous aa was a part of alpha helix
		coor.X = 0
		coor.Y = -math.Sqrt(aaDist*aaDist - pitch*pitch)
		coor.Z = prevCoor.Z + pitch
	} else if prevLoc == 4 {
		// previous aa was a part of alpha helix
		coor.X = math.Sqrt(aaDist*aaDist - pitch*pitch)
		coor.Y = 0
		coor.Z = prevCoor.Z + pitch
	} else {
		panic("FindXYQuadrant function is returning a value it shouldn't")
	}

	// the distances between every two amino acids is not necessarily equal to aaDist
	return coor
}

// function initially made with alpha helices in mind to categorize points
// according to the axis they lie on
func FindXYQuadrant(coor Coordinate) int {

	// looking along z axis to find the quadrant
	if coor.X == 0 && coor.Y == 0 {
		// on the z axis
		return 0
	} else if coor.X > 0 && coor.Y == 0 {
		// quad 1 and 4 boundary (angle 0)
		return 1
	} else if coor.X == 0 && coor.Y > 0 {
		// quad 1 and 2 boundary (angle 90)
		return 2
	} else if coor.X < 0 && coor.Y == 0 {
		// quad 2 and 3 boundary (angle 180)
		return 3
	} else if coor.X == 0 && coor.Y < 0 {
		// quad 4 and 5 boundary (angle 270)
		return 4
	}
	return -1 // shouldnt happen
}

// ____________SECOND 3D VISUALIZATION ATTEMPT (AND 2D VISUALIZATION) FUNCTIONS BELOW____________

// function to write the coordinates to a PDB file
func WriteCoordinatesToPDB(coordinates []Coordinate, fileName string) error {
	file, err := os.Create(fileName)
	if err != nil {
		return err
	}
	defer file.Close()

	for i, coord := range coordinates {
		// Format follows the ATOM record type format for PDB
		fmt.Fprintf(file, "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n",
			// show as water molecule instead of carbon atom
			// fmt.Fprintf(file, "ATOM  %5d  O   HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n",

			i+1, i+1, coord.X, coord.Y, coord.Z)
	}
	return nil
}

// function to make a proxy secondary structure prediction array for testing visualization
func TestVisualization(protein Protein, predArray []CFScore) []int {
	// Testing visualization
	// initialize the aaSecStruct int slice with size of the protein sequence
	aaSecStruct := make([]int, len(protein.Sequence))
	// for each tuple in the predArray, if the first element is max, then append 1 to the aaSecStruct slice
	// if the second element is max, then append 2 to the aaSecStruct slice and so on
	for itr, tuple := range predArray {
		if tuple.Helix > tuple.Sheet && tuple.Helix > tuple.Loop {
			aaSecStruct[itr] = 1
		} else if tuple.Sheet > tuple.Helix && tuple.Sheet > tuple.Loop {
			aaSecStruct[itr] = 2
		} else {
			aaSecStruct[itr] = 3
		}
	}
	return aaSecStruct
}

// function to convert the reassignedABHelixSheet slice to aaSecStruct slice
func ConvertABHelixSheetToAASecStruct(reassignedABHelixSheet []ABHelixSheet, aaSecStruct []int) []int {
	// iterate over the reassignedABHelixSheet slice
	for _, abHelixSheet := range reassignedABHelixSheet {
		// iterate over the positions of the helix or sheet
		for i := abHelixSheet.StartIndex; i <= abHelixSheet.EndIndex; i++ {
			// assign the appropriate value to the aaSecStruct slice
			if abHelixSheet.typeAB == "helix" {
				aaSecStruct[i] = 1
			} else if abHelixSheet.typeAB == "sheet" {
				aaSecStruct[i] = 2
			}
		}
	}
	return aaSecStruct
}

// function to make a 2d plot for the prediction output
func Make2DPlot(predSeq []int, output string) {
	l := len(predSeq)
	lineWidth := 5                            // width of the vertical lines in the figure
	spaceBuffer := 50                         // to have buffer space on left and right
	figWidth := 2*spaceBuffer + 2*lineWidth*l // keeping space between lines same as linewidth

	// initialize the plot
	dc := gg.NewContext(figWidth, 200)
	// set the background color
	dc.SetRGB255(255, 255, 255) // make the whole fig white initially
	dc.Clear()

	// set the font
	dc.LoadFontFace("/Library/Fonts/Arial.ttf", 96)

	// set the line width for all lines
	dc.SetLineWidth(float64(lineWidth))

	for index := range predSeq {
		// value of xCoor, given the values of index, figWidth, buffer, lineWidth
		xCoor := spaceBuffer + index*lineWidth*2

		// draw one vertical line
		MakeVerticalLine(predSeq, index, xCoor, spaceBuffer, dc)
		// print every 10th line's coordinates (or the end coordinate)
		if index == l-1 || index%10 == 0 {
			dc.SetRGB255(0, 0, 0)
			// subtracted linewidth from the xCoor to shift text slightly to left (better this way)
			dc.DrawString(strconv.Itoa(index), float64(xCoor-lineWidth), float64(dc.Height())-float64(spaceBuffer)/2)
		}
	}

	// print the legend:
	dc.DrawString("helix is red, sheet is green, loop is blue", float64(spaceBuffer), float64(spaceBuffer)/2)

	// save the plot
	dc.SavePNG(output)

	// print success message
	fmt.Println("Successfully created the png (2D-visualization) file/s")
}

// function to make a vertical line for one amino acid in the prediction output of the appropriate
// color according to the secondary structure
func MakeVerticalLine(predSeq []int, index int, xCoor int, spaceBuffer int, dc *gg.Context) {
	// set the color for the plot
	if predSeq[index] == 1 {
		// red for helix
		dc.SetRGB255(255, 0, 0)
	} else if predSeq[index] == 2 {
		// green for sheet
		dc.SetRGB255(0, 255, 0)
	} else if predSeq[index] == 3 {
		// blue for loop
		dc.SetRGB255(0, 0, 255)
	} else {
		panic("bad prediction (not helix/sheet/loop)")
	}

	// draw the line
	dc.DrawLine(float64(xCoor), float64(0+spaceBuffer), float64(xCoor), float64(dc.Height()-spaceBuffer))
	dc.Stroke()
}

// 3D VISUALIZATION FUNCTIONS BELOW

// function to make a 3D plot for the prediction output.
// input is the input pdb file's name. output is the output html file's name
func Make3DPlot(predSeq []int, inputPDB string, outputHTML string) {

	// loop through the predSeq and assign values to colr, colg, colb according to its type
	// and fill up the colors slice
	colorInput := make([]string, 0)
	for index := range predSeq {
		var colr int
		var colg int
		var colb int
		if predSeq[index] == 1 {
			// red for helix
			colr = 255
			colg = 0
			colb = 0
		} else if predSeq[index] == 2 {
			// green for sheet
			colr = 0
			colg = 255
			colb = 0
		} else if predSeq[index] == 3 {
			// blue for loop
			colr = 0
			colg = 0
			colb = 255
		} else {
			panic("bad prediction (not helix/sheet/loop)")
		}
		// append to the colorInput slice, something of the following type:
		// {{start_residue_number: {}, end_residue_number: {}, color:{{r:{},g:{},b:{}}}}}
		colorInput = append(colorInput, fmt.Sprintf("{start_residue_number: %d, end_residue_number: %d, color:{r:%d,g:%d,b:%d}}", index, index, colr, colg, colb))
	}
	// convert the colorInput slice to a string
	result := "[" + strings.Join(colorInput, ", ") + "]"

	// read the html file into a toEdit slice
	toEdit := ReadFile("3d_visualization_resources/final_html_template.html")

	// append the resultant color string to the toEdit slice (which is the raw html file)
	toEdit[144] = result + "})}) \n"

	// read the pdb file into a slice of strings
	pdbFile := ReadFile(inputPDB)

	// writing to a pdb file
	line := 205 // where to begin writing pdb info within the html file

	// Insert pdbArray elements into toEdit starting at line 205
	for i := range pdbFile {
		toEdit = append(toEdit[:line+i], append([]string{pdbFile[i]}, toEdit[line+i:]...)...)
	}

	// write the toEdit slice to a new html file
	WriteToFile(outputHTML, toEdit)

	// print success message
	fmt.Println("Successfully created the html (3D-visualization) file/s")
}

// function to open a file named output, and write the contents of the toEdit slice to it
func WriteToFile(output string, toEdit []string) {
	// Open the file for writing
	file, err := os.Create(output)
	if err != nil {
		fmt.Println("Error creating file:", err)
		return
	}
	defer file.Close()

	// Write the contents of the toEdit slice to the file
	for _, line := range toEdit {
		fmt.Fprintln(file, line)
	}
}

// function to take an input string, and read the whole file line by line and store in a slice of strings
func ReadFile(input string) []string {
	// Open the file for reading
	file, err := os.Open(input)
	if err != nil {
		fmt.Println("Error opening file:", err)
		return nil
	}
	defer file.Close()

	// Read the file line by line
	var toEdit []string
	scanner := bufio.NewScanner(file)
	// increase the max token size to accomodate the long lines in the html file
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 4096*1024)
	for scanner.Scan() {
		toEdit = append(toEdit, scanner.Text())
	}

	// Check for errors during scanning
	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading file:", err)
		return nil
	}
	return toEdit
}

// MakeEnsemblArray is a function that takes as its input a string of Ensembl IDs separated by spaces
// and returns a slice of strings containing the Ensembl IDs
func MakeEnsemblArray(ensemblInput string) []string {
	// split the input string into a slice of strings
	ensemblIDs := strings.Split(ensemblInput, ",")
	return ensemblIDs
}

// EnsemblToProteinSlice is a function that takes as its input a slice of strings containing Ensembl IDs
// and returns a slice of Protein objects
func EnsemblToProteinSlice(ensemblIDs []string) []Protein {
	// make a slice of proteins
	proteins := make([]Protein, 0)
	// iterate over the ensembl IDs
	for _, ensemblID := range ensemblIDs {
		// make a protein object using the ensembl ID
		proteins = append(proteins, EnsemblToProtein(ensemblID))
	}
	return proteins
}

// EnsemblToProtein is a function that takes as its input a string containing an Ensembl ID
// and returns a Protein object
func EnsemblToProtein(ensemblID string) Protein {
	// make a protein object
	var protein Protein
	// assign the ensembl ID to the protein object
	protein.Identifier = ensemblID
	protein.Sequence = EnsemblToSequence(ensemblID)
	// return the protein object
	return protein
}

// EnsemblToSequence is a function that takes as its input a string containing an Ensembl ID
// and returns a string containing the sequence of the protein
func EnsemblToSequence(ensemblID string) string {
	// in the console, run gget command with the translate flag as true
	// this command basically: gget seq --translate ENSG00000130234 > ensemblID.fasta
	// this will give you the sequence of the protein in the console, we need to save it in a file named ensemblID.fasta
	cmd := exec.Command("gget", "seq", "--translate", ensemblID)
	output, err := cmd.Output()
	if err != nil {
		fmt.Println(err)
	}

	err = os.WriteFile("outputs/fastas/"+ensemblID+".fasta", output, 0644)
	if err != nil {
		fmt.Println(err)
	}

	// make a reader for the file
	reader := GenerateFASTAReader("outputs/fastas/" + ensemblID + ".fasta")
	// read the protein sequence from the reader
	sequence := ReadProteinsFASTA(reader)

	// return the sequence of the first protein in the slice
	return sequence[0].Sequence
}

// EnsemblToUniprotSlice is a function that takes as its input a slice of strings containing Ensembl IDs
// and returns a slice of strings containing Uniprot IDs
func EnsemblToUniprotSlice(ensemblIDs []string) []string {
	// make a slice of strings
	uniprotIDs := make([]string, 0)
	// iterate over the ensembl IDs
	for _, ensemblID := range ensemblIDs {
		// make a uniprot ID using the ensembl ID
		uniprotIDs = append(uniprotIDs, EnsemblToUniprot(ensemblID))
	}
	return uniprotIDs
}

// EnsemblToUniprot is a function that takes as its input a string containing an Ensembl ID
// and returns a string containing a Uniprot ID
func EnsemblToUniprot(ensemblID string) string {
	// a file named ensemblID.fasta will already be present in the directory
	// make a reader for the file
	reader := GenerateFASTAReader("outputs/fastas/" + ensemblID + ".fasta")
	// make a scanner for the reader
	scanner := bufio.NewScanner(reader.file)
	// scan the first line of the file
	scanner.Scan()
	// split the first line into a slice of strings, separated by spaces
	firstLine := strings.Split(scanner.Text(), " ")
	// the uniprot ID is the string, between "uniprot_id: " and " ensembl_id"
	// extract the uniprot ID from the first line, which looks like:
	// >ENST00000252519 uniprot_id: Q9BYF1 ensembl_id: ENST00000252519 gene_name: ACE2 organism: Homo sapiens sequence_length: 805
	uniprotID := firstLine[2]
	// return the uniprot ID
	fmt.Println("uniprot for ensembl: " + ensemblID + " is: " + uniprotID)
	return uniprotID
}

// MakeOutputNames is a function that takes as its input a slice of Protein objects
// and returns a slice of strings containing strings of the form Protein_1, Protein_2, etc.
func MakeOutputNames(proteins []Protein) []string {
	// make a slice of strings
	outputNames := make([]string, 0)
	// iterate over the proteins
	for i := range proteins {
		// make a string of the form Protein_1, Protein_2, etc.
		outputNames = append(outputNames, "Protein_"+strconv.Itoa(i+1))
	}
	return outputNames
}

// DownloadPDBFiles is a function that takes as its input a slice of strings containing Uniprot IDs
// and downloads the PDB files for the proteins in the 3d_visualization_resources folder with the same name
func DownloadPDBFiles(uniprotIDs []string) {
	// iterate over the uniprot IDs
	for _, uniprotID := range uniprotIDs {
		// download the PDB file for the uniprot ID
		DownloadPDBFile(uniprotID)
	}
}

// DownloadPDBFile is a function that takes as its input a string containing a Uniprot ID
// and downloads the PDB file for the protein in the 3d_visualization_resources folder with the same name
func DownloadPDBFile(uniprotID string) {
	// make a string containing the URL for the PDB file
	// The format of URL: https://alphafold.ebi.ac.uk/files/AF- + uniprot_ID + -F1-model_v4.pdb

	url := "https://alphafold.ebi.ac.uk/files/AF-" + uniprotID + "-F1-model_v4.pdb"
	// download the PDB file from the URL
	DownloadFile("3d_visualization_resources/"+uniprotID+".pdb", url)
}

// DownloadFile is a function that takes as its input a string containing the name of the file to be downloaded
// and a string containing the URL of the file to be downloaded
func DownloadFile(filepath string, url string) {
	// make a HTTP request to the URL
	resp, err := http.Get(url)
	if err != nil {
		panic(err)
	}
	defer resp.Body.Close()

	// create the file
	out, err := os.Create(filepath)
	if err != nil {
		panic(err)
	}
	defer out.Close()

	// write the body of the HTTP request to the file
	_, err = io.Copy(out, resp.Body)
	if err != nil {
		panic(err)
	}
}

// MakeDNAOutputNames is a function that takes as its input a slice of Protein objects
// and returns a slice of strings containing strings of the form DNA_1, DNA_2, etc.
func MakeDNAOutputNames(proteins []Protein) []string {
	// make a slice of strings
	outputNames := make([]string, 0)
	// iterate over the proteins
	for i := range proteins {
		// make a string of the form DNA_1, DNA_2, etc.
		outputNames = append(outputNames, "DNA_"+strconv.Itoa(i+1))
	}
	return outputNames
}

// ALL VISUALIZATION FUNCTIONS END HERE

//Functions for identifying turns and for GOR start here.

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

	//Now that we have created a list of turns we merge the turns into the secondary structure
	//and return the new secndary structure.
	return MergeTurns(turns, secStruc)
}

// Functions to identify turns given a slice of turns and an existing ABHelixSheet data structure.
// It returns a new ABHelixSheet with turns included in the secondary structure.
// Input: A slice of turn data structures and an existing ABHelixSheet that has already identified
// Helices and Sheets.
// Output: A new modified ABHelixSheet datastructure with incorporated turns into the existing
// secondary structure so that we have identified secondary structures for coils, turns, helices
// and beta sheets.
func MergeTurns(identifiedTurns []Turn, secStruc []ABHelixSheet) []ABHelixSheet {

	//Range through the slice of identified turns.
	for i := range identifiedTurns {
		//Flag to break thorugh the range of the secondary structure.
		flagBreak := false

		//Set the index for the turn.
		indTurn := identifiedTurns[i].Index
		//Break label to break out of for loop once found the correct sheet.
	out:
		//Range through secondary structures.
		for j, structure := range secStruc {
			//Check if the index of the turn is after the end index of structure.
			if indTurn > structure.EndIndex {
				//Continue if the index is after the end index.
				continue
			}

			//If the structure is not a sheet continue to the next structure.
			if structure.typeAB != "sheet" {
				continue
			}

			//Create a new ABHelixSheet structure that is a turn.
			newTurn := new(ABHelixSheet)
			newTurn.StartIndex = indTurn
			newTurn.EndIndex = indTurn
			newTurn.typeAB = "Turn"

			//Check to see if the turn is in any of the secondary structure element.
			if inSecondStructure(indTurn, structure) {
				//Three cases
				//Index of the turn is at start index.
				if indTurn == structure.StartIndex {
					//Check to see if we are on the first structure.
					if j == 0 {
						//Increment the adjacent secondary structure start index by one.
						secStruc[j].StartIndex += 1

						//Append the new turn to the secondary structure.
						secStruc = append([]ABHelixSheet{*newTurn}, secStruc...)
					} else {
						//If we are not on the first structure we do the following.

						//Split the secondary structure slices into two arrays.
						upToTurn := make([]ABHelixSheet, len(secStruc[:j]))
						afterTurn := make([]ABHelixSheet, len(secStruc[j:]))
						upToTurn = CopySecondaryStructure(secStruc[:j])
						//Decrease the previous secondary structure by 1.
						upToTurn[len(upToTurn)-1].EndIndex -= 1
						afterTurn = CopySecondaryStructure(secStruc[j:])

						//Increase the start index after we insert the turn by 1.
						afterTurn[0].StartIndex += 1

						//Append/Insert the new turn into the secondary structure.
						secStruc = append(upToTurn, *newTurn)
						secStruc = append(secStruc, afterTurn...)
					}
				} else if indTurn == structure.EndIndex {
					//Index of the turn is at the end of the start index.

					//Create two slices,one up to where we will insert the new turn and after
					//the new turn.
					upToTurn := make([]ABHelixSheet, 0)
					afterTurn := make([]ABHelixSheet, 0)
					//Check to see if we are on the first structure.
					if j != 0 {
						//If we are not on the first structure we copy up until the current sheet
						//to the slice.
						upToTurn = CopySecondaryStructure(secStruc[:j+1])
						//Decrease the current sheet end index by one.
						upToTurn[j].EndIndex -= 1
					} else {
						//If we are on the first structure we only need the first sheet/first structure.
						upToTurn = []ABHelixSheet{secStruc[j]}
						//Decrease the end index of the current sheet by one.
						upToTurn[j].EndIndex -= 1
					}

					//Copy the secondary structure up to current sheet.
					afterTurn = CopySecondaryStructure(secStruc[j+1:])
					//Append the new turn to the copy.
					secStruc = append(upToTurn, *newTurn)
					//Append the remaining secondary structure to the copy and the turn.
					secStruc = append(secStruc, afterTurn...)
				} else {
					//split the secondary structure into two different sheets.

					//Create two ABHelixSheet data structures.
					firstStruc := *new(ABHelixSheet)
					secondStruc := *new(ABHelixSheet)

					//Create copies of the sheet we want to split in half.
					firstStruc = CopyHelixSheet(secStruc[j])
					secondStruc = CopyHelixSheet(firstStruc)

					//Adjust their start and end indeces respectively.
					firstStruc.EndIndex = indTurn - 1
					secondStruc.StartIndex = indTurn + 1

					//Add a panic to check that the indices are the right length.
					if secondStruc.StartIndex > secondStruc.EndIndex {
						panic("This shouldn't happen! Start index is bigger!")
					}

					//Create two slices of up to the turn and after the turn.
					upToTurn := make([]ABHelixSheet, len(secStruc[:j]))
					afterTurn := make([]ABHelixSheet, len(secStruc[j:]))

					//Copy the slice of the secondary structure up to the turn and after the turn.
					upToTurn = CopySecondaryStructure(secStruc[:j])
					afterTurn = CopySecondaryStructure(secStruc[j+1:])

					//Append the uptoturn with the first split of the secondary structure.
					secStruc = append(upToTurn, firstStruc)
					//Append the turn after the first split of the secondary structure.
					secStruc = append(secStruc, *newTurn)
					//Append the second split after the turn of the secondary structure.
					secStruc = append(secStruc, secondStruc)
					//Append the remaining elements to the uptoturn + firstsplit + turn + secondsplit.
					secStruc = append(secStruc, afterTurn...)

				}

				//Set the flag to true as it means that we have inserted our turn, because the turn
				//was found in the secondary structure.
				flagBreak = true
			}

			//Checks to see if the flag is true and if it is we don't have to keep searching for
			//the position to insert our turn into.
			if flagBreak {
				//Jump to the label out.
				break out
			}
		}

	}
	//Return the new secondary structure list, with turns inserted in the proper areas.
	return secStruc
}

// Function to copy secondary structure ABHelixSheet slice.
// Input: A slice of ABHelixSheet.
// Output: A value copy of the slice.
func CopySecondaryStructure(secStruc []ABHelixSheet) []ABHelixSheet {
	//Create a slice of the same length.
	newSecStruc := make([]ABHelixSheet, len(secStruc))
	//Range through the slice.
	for i := range secStruc {
		//Copy each individual data structure, ABHelixSheet using the CopyHelixSheet function.
		newSecStruc[i] = CopyHelixSheet(secStruc[i])
	}

	//Return the new secondary sturcture.
	return newSecStruc
}

// A function to copy secondary structure ABHelixSheet.
// Input: A specific ABHelixSheet data structure.
// Output: A value copy of that data structure.
func CopyHelixSheet(helixOrSheet ABHelixSheet) ABHelixSheet {
	//Creates a new variable/ABHelixSheet.
	newHelixSheet := new(ABHelixSheet)
	//Copies the parameters directly (there are no slice parameters.)
	newHelixSheet.typeAB = helixOrSheet.typeAB
	newHelixSheet.Score = helixOrSheet.Score
	newHelixSheet.StartIndex = helixOrSheet.StartIndex
	newHelixSheet.EndIndex = helixOrSheet.EndIndex

	//Returns the copy of the ABHelixSheet.
	return *newHelixSheet
}

// A function that will serve to help if a turn can be inserted into a sheet or another type
// of structure.
// Input: An integer index and a given ABHelixSheet data structure.
// Output: A boolean value if the index is within the ABHelixSheet range.
func inSecondStructure(index int, betaSheet ABHelixSheet) bool {
	//Define the start and end indeces of the given ABHelixSheet (presumably a sheet.)
	startInd, endInd := betaSheet.StartIndex, betaSheet.EndIndex
	//Checks to see if the index is in the correct range.
	if startInd <= index && index <= endInd {
		//Returns true if it meets these conditions.
		return true
	}
	//Otherwise returns false.
	return false
}

// Function to assess accuracy of predicted vs real string of secondary structures.
// Input: A prediction secondary structure string and a real secondary structure string.
// Output: Four float values, an accuracy float, a helix accuracy float, a sheet accuracy float,
// and a coil accuracy float.
func assessAccuracy(predSS, realSS string) (float64, float64, float64, float64) {
	//Check to see if the lengths of the two strings are the same.
	if len(predSS) != len(realSS) {
		print(len(predSS), len(realSS))
		panic("Not the same length secondary structure!")
	}

	//Create various floats for counting number of correctly identified secondary structures
	//and total secondary structures. One for each ss, helices, sheets, and coils.
	var count float64
	var countCorrectHelices, totalHelices float64
	var countCorrectSheets, totalSheets float64
	var countCorrectCoils, totalCoils float64

	//Range through the predicted secondary structure.
	for i := range predSS {
		//Check to see if the real and the predicted are the same.
		if predSS[i] == realSS[i] {
			//If they are increment the regular count.
			count++
			//Check to see what kind of secondary structure that was correctly identified is.
			if string(predSS[i]) == "H" {
				//Increment countCorrect helices.
				countCorrectHelices++
			} else if string(predSS[i]) == "S" {
				//Increment countCorrectSheets.
				countCorrectSheets++
			} else {
				//Otherwise you increment coils.
				countCorrectCoils++
			}
		}

		//If the predicted is not the same as the real, then wee just want to increment total
		//number /counts of those specific secondary structures.
		if string(realSS[i]) == "H" {
			//Increment total helices count.
			totalHelices++
		} else if string(realSS[i]) == "S" {
			//Increment total sheets count.
			totalSheets++
		} else {
			//Increment total coils count.s
			totalCoils++
		}

	}

	//Find total accuracy.
	acc := count / float64(len(predSS))

	//Find the individual accuracies of the various secondary structures and return as four
	//float values.
	accH := countCorrectHelices / totalHelices
	accS := countCorrectSheets / totalSheets
	accC := countCorrectCoils / totalCoils

	//Return the accuracies.
	return acc, accH, accS, accC
}

// Function to convert an integer slice to a string. This is a temporary datatype to convert
// ABHelixSheets to a slice of integers to an evnetual string.
// Input: A slice of integers.
// Output: A secondary structure string.
func convertAASecStrucToString(aaSecStruc []int) string {
	//Create a variable for containing the whole string.
	totalString := ""

	//Range through the slice of integers.
	for _, val := range aaSecStruc {
		//If the value is two, we add an S to the string.
		if val == 2 {
			totalString += "S"
		} else if val == 1 {
			//Otherwise we add an H for a helix.
			totalString += "H"
		} else {
			//If neither of those occurs then we add a loop or L.
			totalString += "L"
		}
	}

	//We then return the total string.
	return totalString
}

// Function that writes out the accuracies of the validation test sets given the directory of
// the validation tests, the gor parameters, and an index map for gor and the various parameters.
// Input: A directory string for the validation tests, a 2d float64 slice for the parameters, and
// an index map for the gor parameters.
// Output: Writes out a csv file that contains the accuracy for each different secondary structure.
func TestGorAndFasman(dir string, parameters [][]float64, aaIndexMap map[rune]int) {
	//Range through all of the results validation tests.
	files := ReadDirectory(dir)
	//Create a prediction array for both GOR and Chou-Fasman.
	ProteinPredArray := make([][]CFScore, len(files))
	ProteinPredGorArray := make([][]CFScore, len(files))

	//Read in the parameters for GOR from the validation test set (SecNet 2018)
	alphaParam := GORMethodInput("GorParams/alpha_helix.txt")
	sheetParam := GORMethodInput("GorParams/beta_strand.txt")
	turnParam := GORMethodInput("GorParams/beta_turn.txt")
	coilParam := GORMethodInput("GorParams/coil.txt")

	//Open the csv file and if there is an error panic.
	file, err := os.Create("results/Output/results.csv")
	//If error, panic.
	if err != nil {
		panic(err)
	}
	defer file.Close()
	//Print out the protein name, the accuracies for chou fasman, and then the accuracy for
	//Gor out into the csv file for header information.
	fmt.Fprintf(file, "Protein Name, accCF, accCFHelices, accCFSheets, accCFCoils, accGor, accGorHelices, accGorSheets, accGorCoils")
	//Range thorugh the different validation sets.
	for i, inputFile := range files {
		//Do Chou-Fasman on this input file.
		//Create a protein from this input file.
		var newProt Protein

		//Create a variable for the actual structure string.
		var actualStructure string
		newProt.Sequence, actualStructure = GORReadTests(dir + "/" + inputFile.Name())
		//Set the identifier using the input file name's first five letters.
		newProt.Identifier = inputFile.Name()[0:5]

		//Take the structure and map the structure to a comparable format.
		actualStructure = mapSSToLocalStructure(actualStructure)

		//Create the CFScore slice using the Chou-Fasman function method.
		ProteinPredArray[i] = ChouFasman(newProt, parameters, aaIndexMap)
		//Predict helices.
		helices := IdentifyHelicies(ProteinPredArray[i])

		//Predict beta sheets.
		betaSheets := IdentifySheets2(ProteinPredArray[i])

		//Reassign the helices and sheets as appropriate.
		reassignedABHelixSheet := AHelicalBSheetAssignment(append(helices, betaSheets...))
		reassignedABHelixSheet = FillGapsInSequence(len(ProteinPredArray), reassignedABHelixSheet)

		//Reassign structures for turns.
		reassignedABHelixSheet = IdentifyTurns(newProt, parameters, aaIndexMap, reassignedABHelixSheet)

		//Create a slice of integers.
		aaSecStructCFSingle := make([]int, len(newProt.Sequence))

		// Make int slice, which would store the information about the secondary structure of each position
		// (1 = helix, 2 = sheet, 3 = loop).

		// It should be initialized with 3s (as if not helix or sheet, then loop) and we assign helix and sheet below
		for i := 0; i < len(aaSecStructCFSingle); i++ {
			aaSecStructCFSingle[i] = 3
		}

		//Convert the reassignedABHelixSheet slice to aaSecStruct slice.
		aaSecStructSingle := ConvertABHelixSheetToAASecStruct(reassignedABHelixSheet, aaSecStructCFSingle)

		//Convert this aaSecStruct into a string that can be used to compare to the actual structure
		//string.
		cfstr := convertAASecStrucToString(aaSecStructSingle)
		//Generate the accuracies for chou-fasman!
		accCF, accCFHelices, accCFSheets, accCFCoils := assessAccuracy(cfstr, actualStructure)

		//Do GOR on this same file!
		//Create a 3d array of parameters for GOR.
		params := make([][][]float64, 4)

		//Set the various parameters for each different secondary structure.
		params[0], params[1], params[2], params[3] = alphaParam, sheetParam, turnParam, coilParam
		//Create an index map using the aa_index for gor parameters.
		aaIndexMap := ReadAAIndexMap("aa_index_map_gor.txt")

		//Create a protein prediction array using GORPrediciton method. This hijacks the CFScore data
		//structure to streamline the use of data-structures.
		ProteinPredGorArray[i] = GORPrediction(newProt, params, aaIndexMap)
		//Converts this temporary CFScore into an int array and then into a slice of ABHelixSheets.
		//This will then allow the previously outlined method for GOR to convert back into a slice of
		//integers.

		//Convert to an ABHelix Slice.
		ssStruc := GorPredictionConv(ConvertPredToArr(ProteinPredGorArray[i]))
		// Make int slice, which would store the information about the secondary structure of each position
		// (1 = helix, 2 = sheet, 3 = loop).
		//Create a slice of integers the same length as the protein sequence.
		aaSecStructGORSingle := make([]int, len(newProt.Sequence))

		// It should be initialized with 3s (as if not helix or sheet, then loop) and we assign helix and sheet below
		for i := 0; i < len(aaSecStructGORSingle); i++ {
			aaSecStructGORSingle[i] = 3
		}

		// Convert the reassignedABHelixSheet slice to aaSecStruct slice.
		aaSecStructSingle = ConvertABHelixSheetToAASecStruct(ssStruc, aaSecStructGORSingle)

		//Convert this into a string with S, H, or L.
		gorstr := convertAASecStrucToString(aaSecStructSingle)
		//Find the accuracies for GOR.
		accGor, accGorHelices, accGorSheets, accGorCoils := assessAccuracy(gorstr, actualStructure)

		//Write out all the accuracies into the csv file.
		//Convert all the accuracies into strings.
		accCFStr := strconv.FormatFloat(accCF, 'f', -1, 64)
		accCFHelicesStr := strconv.FormatFloat(accCFHelices, 'f', -1, 64)
		accCFSheetsStr := strconv.FormatFloat(accCFSheets, 'f', -1, 64)
		accCFCoilsStr := strconv.FormatFloat(accCFCoils, 'f', -1, 64)
		accGorStr := strconv.FormatFloat(accGor, 'f', -1, 64)
		accGorHelicesStr := strconv.FormatFloat(accGorHelices, 'f', -1, 64)
		accGorSheetsStr := strconv.FormatFloat(accGorSheets, 'f', -1, 64)
		accGorCoilsStr := strconv.FormatFloat(accGorCoils, 'f', -1, 64)

		//Print a new line into the csv file.
		fmt.Fprintf(file, "\n")
		//Print out all the accuracies to the csv file.
		fmt.Fprintf(file, newProt.Identifier+","+accCFStr+","+accCFHelicesStr+","+accCFSheetsStr+","+accCFCoilsStr+","+accGorStr+","+accGorHelicesStr+","+accGorSheetsStr+","+accGorCoilsStr)
	}

}

// This functtion takes in the inputs from the GorParams folder and turns them into a prediction
// protein string and the real secondary structure string.
// Input: A file string.
// Output: Two strings, one being the protein sequence, and one being the actual secondary structure
// of the protein.
func GORReadTests(file string) (string, string) {
	//Open the file.
	f, err := os.Open(file)
	if err != nil {
		//Panic if it can't open.
		log.Panicf("Cannot open file: %s", err)
	}
	//Create two new string variables.
	var ss string
	var seq string

	//Scan the file.
	scanner := bufio.NewScanner(f)

	//Create a counter variable.
	var counter int
	//Scan line by line.
	for scanner.Scan() {
		//Read in the line into a variable,
		line := scanner.Text()
		//Check to see if on first line.
		if counter == 0 {
			//If on first line we can set the protein sequence variable.
			seq = line
		}
		//If we are on second line check.
		if counter == 1 {
			//If we are on the second line we can set the real secondary structure to a variable.
			ss = line
		}
		//Increment the line counter.
		counter++
	}

	//Close the file after scanning.
	defer f.Close()

	//Return the two strings.
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

// Function that maps the sequence from the SECNET 2018 validation file into the same local
// program format.
// Input: A string that contains the Secnet 2018 secondary structure.
// Output: A locally formatted string that is in the same format as the predicted string.
func mapSSToLocalStructure(sequenceFromSS string) string {
	//Create a new string for the secondary structure.
	var newSequence string
	//Range through the real sequence from Secnet.
	for i := range sequenceFromSS {
		//IF the letter is a C for coil map to loop (treated the same as coil)
		if string(sequenceFromSS[i]) == "C" {
			newSequence += "L"
		} else if string(sequenceFromSS[i]) == "E" {
			//IF letter is E treat it as sheet.
			newSequence += "S"
		} else {
			//Otherwise treat it as a helix.
			newSequence += "H"
		}
	}
	//Return the local format secondary structure string.
	return newSequence
}

//Identify turn functions and GOR functions end here!
