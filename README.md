# expansion-cones



## Building

`expansion-cones` uses `cmake` for compilation.

### Requirements
This project relies on the following libraries:

* CGAL, for LP-solving and rational numbers representation
* OpenVolumeMesh, for everything mesh-related
* Eigen3

Make sure to install those in a way that `cmake` can find them with `find_package`


## Usage

The project generates an executable called `ShrinkAndExpand`. You can simply run it without any argument to see its usage.
Note that the input meshes for this project are in the `.ovm` format (OpenVolumeMesh). 
If you need to convert your mesh from another data format, please refer to `meshio` and Martin's fork supporting .ovm meshes: https://github.com/mheistermann/meshio





## TODOs

* Description
* link to publication / bibtex ref
* figure
* license
