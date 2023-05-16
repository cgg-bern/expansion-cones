# Expansion Cones: A Progressive Volumetric Mapping Framework

[Valentin Z. NIGOLIAN](https://cgg.unibe.ch/), [Marcel CAMPEN](http://graphics.cs.uos.de/), [David BOMMES](https://cgg.unibe.ch/)
ACM Transaction on Graphics (Proceedings of SIGGRAPH 2023)

![](/readme/cover-picture.png)
Our method generates a bijective map of a ball-topology tetrahedral mesh (a) to a star-shaped domain, e.g. a ball (b). Starting with all interior vertices
clustered inside the domain’s kernel (b), we iteratively split clusters by picking a subset of vertices whose 1-ring neighborhood union has a non-empty kernel.
By moving this subcluster into this kernel, some initially degenerate tetrahedra are expanded (c), without degenerating or inverting others. Repeating this
until no cluster remains, all tetrahedron images obtain positive volume, yielding a bijective map (d). Mesh refinement is applied adaptively in the process to
obtain the necessary degrees of freedom. A key invariant is that the intermediate maps never invert any tetrahedron

[`Project Page`](https://www.algohex.eu/publications/expansion-cones/)
[`Dataset`](todo)

## Abstract
Volumetric mapping is a ubiquitous and difficult problem in Geometry Processing and has been the subject of research in numerous and various directions. 
While several methods show encouraging results, the field still lacks a general approach with guarantees regarding map bijectivity. 
Through this work, we aim at opening the door to a new family of methods by providing a novel framework based on the concept of _progressive expansion_.
Starting from an initial map of a tetrahedral mesh whose image may contain degeneracies but no inversions, we incrementally adjust vertex images to expand degenerate elements. 
By restricting movement to so-called _expansion cones_, it is done in such a way that the number of degenerate elements decreases in a strictly monotonic manner, without ever introducing any inversion. 
Adaptive local refinement of the mesh is performed  to facilitate this process.
We describe a prototype algorithm in the realm of this framework for the computation of maps from ball-topology tetrahedral meshes to convex or star-shaped domains.
This algorithm is evaluated and compared to state-of-the-art methods, demonstrating its benefits in terms of bijectivity.
We also discuss the associated cost in terms of sometimes significant mesh refinement to obtain the necessary degrees of freedom required for establishing a valid mapping.
Our conclusions include that while this algorithm is only of limited immediate practical utility due to efficiency concerns, the general framework has the potential to inspire a range of novel methods improving on the efficiency aspect.


## Building

`expansion-cones` uses `cmake` for compilation.

### Requirements
This project relies on the following libraries:

* [CGAL](https://www.cgal.org/), for LP-solving and rational numbers representation
* [OpenVolumeMesh](https://www.graphics.rwth-aachen.de/software/openvolumemesh/), for everything mesh-related
* [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)

Make sure to install those in a way that `cmake` can find them with `find_package`


## Usage

The project generates an executable called `ShrinkAndExpand`. You can simply run it without any argument to see its usage.
Note that the input meshes for this project are in the `.ovm` format (OpenVolumeMesh). 
If you need to convert your mesh from another data format, please refer to `meshio` and Martin Heistermann's fork supporting .ovm meshes: https://github.com/mheistermann/meshio




