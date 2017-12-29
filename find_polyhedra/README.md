# `find_polyhedra`

Identifies ions that define the vertices of &ldquo;coordination&rdquo; polyhedra in a crystal lattice.

Reads in a set of coordinates for lattice ions in a 3D periodic structure, and uses the
nearest neighbours for each ion (all other ions within a set cutoff radius) to construct 
lists of tetrahedra and octahedra.

This code was initially developed for analysis of close-packed lattices, and may give unexpected results for lattices with irregular coordination polyhedra.

This is expected to work best for a geometry optimised &ldquo;zero temperature&rdquo; configuration, and may well fail for thermally disordered structure.

## inputs
`find_polyhedra.inpt`    input parameters (see below)
<coordinate file>   list of (x,y,z) Cartesian coordinates for the lattice ions

## example input

`<coordinates file>`  name of file containg the lattice ion (x,y,z) Cartesian coordinates.
`<natoms>`            number of lattice ions (number of entries in <coordinates file>).
`<ntet>`              number of expected tetrahedra.
`<noct>`              number of expected octahedra.
`<cell_x, cell_y, cell_z>` cell lengths.
`<h (3,3)>`           (3x3) cell matrix of unit vectors for cell axes (as matrix rows).
`<r_cut>`             cutoff radius for constructing ion neighbour lists.
`<cp_tet_assign>`     `.true.` or `.false.`
                      Whether to divide tetrahedra into two sets, according to their 
                      orientation &ldquo;up&rdquo or &ldquo;down&rdquo; with respect to
                      close-packed planes. If `.true.` then the close-packed vector is 
                      read from the next line.
`<x, y, z>`           vector along &ldquo;close-packed&rdquo; direction. Used to define the relative 
                      orientation of &ldquo;up&rdquo; and &ldquo;down&rqduo; oriented tetrahedra. 

e.g.

```
initial_positions.xyz
768                  natoms
49.06723742 49.06723742 49.06723742       
960 ntet
512 noct 
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
7.5                  rcut 
.false. assign tetrahedra according to close-packed direction (if .true. add close-packed vector)
0.0 0.0 1.0          close-packed vector. Only read if the previous line is .true.
```

## output

### to stdout

lists progress, and the numbers of tetrahedra and octahedra found, e.g.

```
 Creating neighbour lists
 searching for tetrahedra
 960 tetrahedra found
 searching for octahedra
 512 octahedra found
```

### output files

For each polyhedron type found, two output files are written:  
`<name>.list`     Ion numbers for the ions defining the polyhedron vertices, as a `n_vertex` &times; `n_polyhedra` array of integers.
`<name>_c.out`    Cartesian coordinates of the centres of each polyhedron, as a `3` &times; `n_polyhedra` array of floats.

e.g.

`oct.list`        Ion numbers defining the vertices of each octahedron
`tet_up.list`     Ion numbers defining the vertices of each tetrahedron in set 1 (pointing "up")
`tet_down.list`   Ion numbers defining the vertices of each tetrahedraon in set 2 (pointing "down")
`oct_c.out`       Cartesian coordinates of the centres of each octahedron
`tet_up_c.out`    Cartesian coordinates of the centres of each tetrahedron in set 1 (pointing "up")
`tet_down_c.out`  Cartesian coordinates of the centres of each tetrahedron in set 2 (pointing "down")
