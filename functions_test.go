package main

import (
	"fmt"
	"io/ioutil"
	"strconv"
	"strings"
	"testing"
)

// test ReadProteinsFASTA function
func TestReadProteinsFASTA(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/ReadProteinsFASTA/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/ReadProteinsFASTA/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)

		reader := GenerateFASTAReader(filePath)
		functionToString := strings.TrimSpace(fmt.Sprint(ReadProteinsFASTA(reader)))
		if outputFile != functionToString {
			t.Errorf("ReadProteinsFASTA() = %v, want %v", strings.TrimSpace(outputFile), functionToString)
		}
	}
}

// test TranslateDNA function
func TestTranslateDNA(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/TranslateDNA/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/TranslateDNA/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)
		reader := GenerateDNAReader(filePath)

		functionToString := fmt.Sprint(TranslateDNA(reader))
		if outputFile != functionToString {
			t.Errorf("TranslateDNA() = %v, want %v", outputFile, functionToString)
		}
	}
}

// test ReadParameters function
func TestReadParameters(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/ReadParameters/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/ReadParameters/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)
		reader := ReadParameters(filePath)
		if outputFile != fmt.Sprint(reader) {
			t.Errorf("ReadParameters() = %v, want %v", outputFile, fmt.Sprint(reader))
		}
	}
}

// test ReadAAIndexMap function
func TestReadAAIndexMap(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/ReadAAIndexMap/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/ReadAAIndexMap/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)
		reader := ReadAAIndexMap(filePath)
		if outputFile != fmt.Sprint(reader) {
			t.Errorf("ReadAAIndexMap() = %v, want %v", outputFile, fmt.Sprint(reader))
		}
	}
}

// test NamesToChar function
func TestNamesToChar(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/NamesToChar/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/NamesToChar/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)
		reader := NameToChar(filePath)
		if outputFile != fmt.Sprint(reader) {
			t.Errorf("NamesToChar() = %v, want %v", outputFile, fmt.Sprint(reader))
		}
	}
}

// test TranslateCodon function
func TestTranslateCodon(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/TranslateCodon/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/TranslateCodon/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)
		inputString, _ := readFileToString(filePath)
		reader := TranslateCodon(inputString)
		if reader.Identifier != outputFile {
			t.Errorf("TranslateCodon() = %v, want %v", outputFile, reader)
		}
	}
}

// test TranscribeDNA function
func TestTranscribeDNA(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/TranscribeDNA/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/TranscribeDNA/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)
		inputString, _ := readFileToString(filePath)
		reader := TranscribeDNA(inputString)
		//fmt.Println("\n", reader)
		if reader != outputFile {
			t.Errorf("TranscribeDNA() = %v, want %v", outputFile, reader)
		}

	}
}

// test ReadCIFToFasta function
func TestReadCIFToFasta(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/ReadCIFToFasta/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/ReadCIFToFasta/output_%d.txt", i)
		// Compare the expected result with the actual result from checkIfInBounds function
		outputFile, _ := readFileToString(outputFilePath)

		reader := GenerateCIFReader(filePath)
		protein, _ := ReadCIFToFasta(reader)

		if outputFile != fmt.Sprint(protein) {
			t.Errorf("ReadProteinsFASTA() = %v, want %v", outputFile, fmt.Sprint(protein))
		}
	}
}

// TestIdentifyHelicies tests the IdentifyHelicies function
func TestIdentifyHelicies(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/IdentifyHelicies/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/IdentifyHelicies/output_%d.txt", i)

		// Read the expected output file
		outputFile, _ := readFileToString(outputFilePath)

		// Read the input file and parse it into CFScoreArray
		inputString, _ := readFileToString(filePath)
		CFScoreArray, _ := parseProteinPredArray(inputString)

		// Run the IdentifyHelicies function and get the result
		reader := IdentifyHelicies(CFScoreArray)

		// Compare the result with the expected output
		if fmt.Sprint(reader) != outputFile {
			t.Errorf("IdentifyHelicies() = %v, want %v", outputFile, fmt.Sprint(reader))
		}
	}
}

// TestIdentifySheets2 tests the IdentifySheets2 function
func TestIdentifySheets2(t *testing.T) {
	// Loop through test cases, using indices from 0 to 3
	for i := 0; i <= 3; i++ {
		// Generate a file path for the input of each test case
		filePath := fmt.Sprintf("Tests/input/IdentifySheets2/input_%d.txt", i)

		// Generate a file path for the expected output of each test case
		outputFilePath := fmt.Sprintf("Tests/output/IdentifySheets2/output_%d.txt", i)

		// Read the expected output file
		outputFile, _ := readFileToString(outputFilePath)

		// Read the input file and parse it into CFScoreArray
		inputString, _ := readFileToString(filePath)
		CFScoreArray, _ := parseProteinPredArray(inputString)

		// Run the IdentifySheets2 function and get the result
		reader := IdentifySheets2(CFScoreArray)

		// Compare the result with the expected output
		if fmt.Sprint(reader) != outputFile {
			t.Errorf("IdentifyHelicies() = %v, want %v", outputFile, fmt.Sprint(reader))
		}
	}
}

func readFileToString(filePath string) (string, error) {
	// Read the file
	data, err := ioutil.ReadFile(filePath)
	if err != nil {
		// If there's an error, return the empty string and the error
		return "", err
	}

	// Convert the byte slice to a string and return it
	return string(data), nil
}

// parseProteinPredArray parses a string into a slice of CFScore structs.
func parseProteinPredArray(input string) ([]CFScore, error) {
	// Trim the leading and trailing brackets and split the string by '} {'
	trimmedInput := strings.Trim(input, "[]")
	splitInput := strings.Split(trimmedInput, "} {")

	var scores []CFScore

	for _, s := range splitInput {
		// Remove any leading or trailing braces or spaces
		trimmedScore := strings.Trim(s, "{} ")

		// Split the trimmed string by space to get individual scores
		scoreParts := strings.Fields(trimmedScore)

		if len(scoreParts) != 3 {
			return nil, fmt.Errorf("invalid score format: %s", s)
		}

		// Convert the string scores to float64 and create a CFScore struct
		helix, err := strconv.ParseFloat(scoreParts[0], 64)
		if err != nil {
			return nil, err
		}
		sheet, err := strconv.ParseFloat(scoreParts[1], 64)
		if err != nil {
			return nil, err
		}
		loop, err := strconv.ParseFloat(scoreParts[2], 64)
		if err != nil {
			return nil, err
		}

		scores = append(scores, CFScore{Helix: helix, Sheet: sheet, Loop: loop})
	}

	return scores, nil
}

func findDiffCharByChar(str1, str2 string) []int {
	var diffs []int
	minLen := len(str1)
	if len(str2) < minLen {
		minLen = len(str2)
	}

	for i := 0; i < minLen; i++ {
		if str1[i] != str2[i] {
			diffs = append(diffs, i)
		}
	}

	return diffs
}
