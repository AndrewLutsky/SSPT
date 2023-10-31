package main

type PredArray []CFScore

type CFScore struct {
	Helix float64
	Sheet float64
	Loop  float64
	F_i   float64
	F_i1  float64
	F_i2  float64
	F_i3  float64
}

type Helix struct {
	StartIndex int
	EndIndex   int
	Score      float64
}

type Turn struct {
	StartIndex int
	EndIndex   int
	Score      int
}

// (Shashank) adding coordinates datatype for visualization

type Coordinate struct {
	X, Y, Z float64
}
