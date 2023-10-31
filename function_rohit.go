// Rohit Nandakumar

package main

import (
	"bufio"
	"log"
	"os"
	"strings"
)

// Protein structure holds the protein data
type Protein struct {
	Sequence   string
	Identifier string
}

type BetaSheet struct {
	StartIndex int
	EndIndex   int
	Score      float64
}

type FASTAReader struct {
	file *os.File
}

func GenerateFASTAReader(fileName string) *FASTAReader {
	file, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}

	// (Shashank) below line needed to be commented out to make the file non empty. Else, the ReadProteins
	// function will return an empty slice
	// defer file.Close()
	return &FASTAReader{file: file}
}

// ReadProteins reads proteins from given FASTA file
func ReadProteins(fastaReader *FASTAReader) []Protein {
	var proteins []Protein
	var currProtein Protein
	scanner := bufio.NewScanner(fastaReader.file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if currProtein.Identifier != "" {
				proteins = append(proteins, currProtein)
				currProtein = Protein{}
			}
			currProtein.Identifier = strings.TrimPrefix(line, ">")
		} else if line != "" {
			currProtein.Sequence += line
		}
	}
	return proteins
}
