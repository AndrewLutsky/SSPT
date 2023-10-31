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
func IdentifyTurns(predArray PredArray) []Turn {
	turns := make([]Turn, 0)

	return turns
}
