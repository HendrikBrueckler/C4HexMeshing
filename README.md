![teaser](https://github.com/HendrikBrueckler/C4HexMeshing/assets/38473042/f70ef0d1-e1f4-4e35-a024-ad9440baf83c)

# Collapsing Cubical Cell Complexes for Hex Meshing

`C4HexMeshing` is an implementation of [Collapsing Embedded Cell Collapses for Safer Hexahedral Meshing \[Brückler et al. 2023\]](http://graphics.cs.uos.de/papers/T-Collapsing_Bru%CC%88ckler_SA2023_Preprint.pdf) (SIGGRAPH Asia 2023) distributed under GPLv3.

If you make use of `C4HexMeshing` in your scientific work, please cite our paper. For your convenience,
you can use the following bibtex snippet:

```bibtex
    @article{C4HexMeshing,
        author     = {Hendrik Br{\"{u}}ckler and
                     Marcel Campen},
        title      = {Collapsing Embedded Cell Collapses for Safer Hexahedral Meshing},
        journal    = {ACM Trans. Graph.},
        volume     = {42},
        number     = {6},
        year       = {2023},
    }
```

***

## How does it work?

`C4HexMeshing` makes use of the [3D Motorcycle Complex](https://github.com/HendrikBrueckler/MC3D) to partition a tetrahedral mesh, equipped with a suitable seamless map, into blocks. A quantization of the cubical cell complex is then computed using [QGP3D](https://github.com/HendrikBrueckler/QGP3D), assigning non-negative integer lengths to the complex' edges.

To extract from this complex a valid integer-grid map and thus a hex mesh, all edges, patches and blocks of the complex which were quantized to zero extent are collapsed and their content redistributed among the remaining elements. On the remaining cell complex, whose cells are now blocks of strictly positive integer extents, an integer-grid map can be computed blockwise via (relatively simple) cube maps.
From this a hex mesh can be extracted.

Note, that (currently) no geometric optimization other than tentative untangling of the integer-grid map is performed, but some might be added in the future.

![Collapses](https://github.com/HendrikBrueckler/C4HexMeshing/assets/38473042/4fbb3c16-1baf-47e7-9d6b-58b7fc8555b3)

***

### Dependencies
- GMP (NOT included, must be installed on your system)
- GUROBI (NOT included, must be installed on your system)
- NLOPT (optional, for map optimization, NOT included, must be installed on your system)
- [MC3D](https://github.com/HendrikBrueckler/MC3D) (Included as submodule, together with all subdependencies)
- [QGP3D](https://github.com/HendrikBrueckler/QGP3D) (Included as submodule, together with all subdependencies)

### Building
In root directory

    mkdir build
    cd build
    cmake -DGUROBI_BASE=<path/to/gurobi/> ..
    make

### Usage
An example command-line application is included that reads a tetrahedral mesh including a seamless parametrization from a file in .hexex-format, as used and documented in [libHexEx](https://www.graphics.rwth-aachen.de/software/libHexEx/).
It can generate the output of several stages of the algorithm, including the original MC, collapsed MC, integer-grid map (unoptimized) and hex mesh (unoptimized).

After building, the CLI app can be found in ```build/Build/bin/cli``` .
For full information on its usage, execute

    c4hex_cli --help

Example input can be found in folder ```extern/QGP3D/extern/MC3D/tests/resources```.

### API
For details on the API of the library, check the headers in ```include```, they are thoroughly documented. Apart from that, ```cli/main.cpp``` demonstrates usage of the entire pipeline for both simple and advanced usage.
