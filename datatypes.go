package main

type PredArray []CFScore

type CFScore struct {
	Helix float64
	Sheet float64
	Loop  float64
}

type Helix struct {
	StartIndex int
	EndIndex   int
	Score      float64
}
