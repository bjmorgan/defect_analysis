# siteid

D. Marrocchelli  
B. J. Morgan

Analyse a molecular dynamics trajectory to identify the sites occupied by mobile ions.

Reads in a simulation trajectory, and assigns ions to lattice sites for each frame.
Lattice sites can be defined as coordination octahedra or tetrahedra (calculated on the
fly using the sub-lattice geometry for each frame) or as spherical sites with fixed
centres (TODO: spherical sites are not implemented yet).

TODO: allowing for flexible site types might be a way to have e.g. two sets of tetrahdral sites (cf. up and down)

## inputs

TODO: document what inputs are needed.

- `defect_new.inpt` TODO: change name to `siteid.inpt`.

TODO: Other inputs? At present:
- file containing list of Cartesian coordinates (e.g. `poscart.out`)
- file containing list of cell geometries (cell lengths and matrices) (e.g. `cell_matrix.out`)
- input files that define polyhedra (e.g. `tet.list` and `oct.list`)

### example input

TODO: this will need to be restructured to allow for more flexible site types
```
poscart-t-LLZO.out
cell_matrix.out
tet.list
oct.list
9200      nsteps
4         nspec
768       natoms
192       natoms
128       natoms
448       natoms
1         species identity that makes up tetrahedral lattice
4         species identity of mobile ions
960       ntet
512       noct
.false.   variable cell
  49.81007063       49.79930828       48.57409755
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
```

## outputs

### to stdout

### output files

`sites_atoms.dat`: atom numbers occupying each site (NOTE: what happens in cases of multiple occupancy?) Empty sites are denoted `0`.
`atoms_sites.dat`: site numbers for each ion (NOTE: what happens if an atom is not assigned to a site?)
`npolyocc.dat`: Numbers of polyhedral types occupied. TODO: this should be more general to include all site types, and also write number of atoms not assigned to sites.

## History

This code was developed out of a series of codes written to identify defects in close-packed
structures. The first version was written by D. Marrocchelli during his DPhil, to identify vacant tetrahedra in fluorite-structured oxide-ion conductors \[1\]. TODO: add other Dario references.
This was adapted to identify octahedral defects in close-packed structured by D. Marrocchelli and B. J. Morgan \[2,3,4\], using the property that in a close-packed structure the set of tetrahedra and octahedra are space-filling, and any ion not in a tetrahedral site must occupy an octahedral site. 
The code was subsequently revised by B. J. Morgan to assign polyhedral occupations in non-close-packed lattices \[5\].
TODO: Addition of spherical site types (as fixed coordinates / centroids of sets of atoms?).

## References

1. D. Marrocchelli, P. A. Madden, S. T. Norberg, and S. Hull. &ldquo;Cation composition effects on oxide conductivity in the Zr<sub>2</sub>Y<sub>2</sub>O<sub>7</sub>-Y<sub>3</sub>NbO<sub>7</sub> system&rdquo; J. Phys.-Condens. Matt. **21** 405403 (2009).
2. B. J. Morgan and P. A. Madden &ldquo;Effects of Lattice Polarity on Interfacial Space Charges and Defect Disorder in Ionically Conducting AgI Heterostructures&rdquo; Phys. Rev. Lett. **107** 206102 (2011).
3. B. J. Morgan and P. A. Madden &ldquo;Absence of a space-charge-derived enhancement of ionic conductivity in &beta;|&gamma;-heterostructured 7H- and 9R-AgI&rdquo; J. Phys.-Conden. Matt. **24** 275303 (2012).
4. B. J. Morgan and P. A. Madden &ldquo;Relationships Between Atomic Diffusion Mechanisms and Ensemble Transport Coefficients in Crystalline Polymorphs&rqduo; Phys. Rev. Lett. **112** 145901 (2014).
5. M. Burbano, D. Carlier, F. Boucher, B. J. Morgan, and M. Salanne, &ldquo;Sparse Cyclic Excitations Explain the Low Ionic Conductivity of Stoichiometric Li<sub>7</sub>La<sub>3</sub>Zr<sub>2</sub>O<sub>12</sub>&rqduo; Phys. Rev. Lett. **116** 135901 (2016).

## inputs
`polyhedra.inpt`    input parameters (see below)
<coordinate file>   list of (x,y,z) Cartesian coordinates for the lattice ions

## example input

`<coordinates file>`  name of file containg the lattice ion (x,y,z) Cartesian coordinates.
`<natoms>`            number of lattice ions (number of entries in <coordinates file>).
`<cell_x, cell_y, cell_z>` cell lengths.
`<h (3,3)>`           (3x3) cell matrix of unit vectors for cell axes (as matrix rows).
`<x, y, z>`           vector along "close-packed" direction. Used to define the relative 
                    orientation of "up" and "down" oriented tetrahedra. 
`<r_cut>`             cutoff radius for constructing ion neighbour lists.

e.g.

```
initial_positions.xyz
768                  natoms
49.06723742 49.06723742 49.06723742        
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
0.0 0.0 1.0          close-packed vector 
7.5                  rcut 
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
 512
 960 tet1 0 tet2 512 oct
```

### output files

`oct.list`    Ion numbers defining the vertices of each octahedron
`tet1.list`   Ion numbers defining the vertices of each tetrahedron in set 1 (pointing "up")
`tet2.list`   Ion numbers defining the vertices of each tetrahedraon in set 2 (pointing "down")
`oct_c.out`   Cartesian coordinates of the centres of each octahedron
`tet1_c.out`  Cartesian coordinates of the centres of each tetrahedron in set 1 (pointing "up")
`tet2_c.out`  Cartesian coordinates of the centres of each tetrahedron in set 2 (pointing "down")
