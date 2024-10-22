# square-hexagon


This repository contains the numericals tools used in the following physics project : https://fr.overleaf.com/read/mjgrwfrzvmhy#facfeb


The goal was to explore the phase transition of the square-hexagon lattice using two numerical methods : Linear Spin Wave Theory (LSWT) and Plaquette Wave Theory (PWT).
It contains :

- "square_hexagon.m" : a Matlab script producing a visual interface that can produce the dispertion relations for both approaches for a given J_1'/J_1 value. 
		For some reason the J'/J slider doesn't always register mouse clicks for some slider values and drag/drop for others on some computers, if one doesn't work, try the other.


- "find_states.m" : a Matlab script giving information about the composition of the Bogoliubov boson for each curve in LSWT. Choose a J'/J value and a BZ point, the code then gives raw data of the computed energy levels at this point as well as the change of basis matrix to the bogoliubov matrix giving information on how the old Holstein-Primakoff bosons compose the new bogoliubov bosons (note this is for the 5x5 (A+B)(A-B) matrix in the report and hence is does not give an actual composition of the bogoliubov bosons).


- "total_energy.m" : given a set of grid sizes N (subdivision of the 1stBZ) and a J_1' value, computes and plots the total energy for LSWT.

- "counting.m" : an attempt at counting the number of each Holstein Primakoff Boson. Given a range of J_1' values will plot the number of bosons for each value, this might take time to process depending on the values. This was not included in the report.
