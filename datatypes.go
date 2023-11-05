package main

import "os"

type PredArray []CFScore

type CFScore struct {
	Helix float64
	Sheet float64
	Loop  float64
}

type ABHelixSheet struct {
	StartIndex int
	EndIndex   int
	Score      float64
	typeAB     string
}

type Turn struct {
	Index int
}

// (Shashank) adding coordinates datatype for visualization

type Coordinate struct {
	X, Y, Z float64
}

// Protein structure holds the protein data
type Protein struct {
	Sequence   string
	Identifier string
}

type FASTAReader struct {
	file *os.File
}
