package main

import (
	"fmt"
	"math"
	"os"
	"strconv"

	"github.com/fogleman/gg"
)

// Hey! Jon here- I commented this out because it was conflicting with functions.go Feel free to uncomment while you're working.
// Sure, No issues!

/*
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
func ChouFasman(seq string, windowSize int, parameters [][]float64, aaIndexMap map[rune]int) []int {

	// Define the prediction array
	predictionArray := make([]int, len(seq)) // 0: not determined, 1: helix, 2: sheet, 3: loop

	if len(seq) < windowSize+1 {
		fmt.Println("Sequence too short")
		return predictionArray
	}

	leftPointer := 0
	rightPointer := windowSize - 1

	for rightPointer < len(seq) {
		predictionArray[leftPointer+windowSize/2] = ChouFasmanWindow(seq[leftPointer:rightPointer+1], parameters, aaIndexMap)
		leftPointer += 1
		rightPointer += 1
	}

	fmt.Println(predictionArray)
	return predictionArray
}
*/

// VISUALIZATION FUNCTIONS BELOW

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
