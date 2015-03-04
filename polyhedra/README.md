# polyhedra

Identifies ions that define the vertices of "coordination" polyhedra in a crystal lattice

Reads in a set of coordinates for lattice ions in a 3D periodic structure, and uses the
nearest neighbours for each ion (all other ions within a set cutoff radius) to construct 
lists of tetrahedra and octahedra.

This code was initially developed for analysis of close-packed lattices, and may give unexpected
results for lattices with irregular coordination polyhedra.

This is expected to work best for a geometry optimised "zero temperature" configuration, and
may well fail for thermally disordered structure.

## inputs
`polyhedra.inpt`    input paramaters (see below)
<coordinate file>   list of (x,y,z) Cartesian coordinates for the lattice ions

## example input

`<coordinates file>`  name of file containg the lattice ion (x,y,z) Cartesian coordinates.
`<natoms>`            number of lattice ions (number of entries in <coordinates file>).
`<cell_x, cell_y, cell_z>` cell lengths.
`<h (3,3)>`           (3x3) cell matrix of unit vectors for cell axes (as matrix rows).
`<x, y, z>`           vector along "close-packed" direction. Used to define the relative 
                    orientation of "up" and "down" oriented tetrahedra. 
`<r_cut>`             cutoff radius for constructing ion neighbour lists.

## output

### to stdout

TODO

### output files

TODO
