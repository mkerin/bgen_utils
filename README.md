# bgen_utils

Collection of utility functions for use with bgen files. Not yet ready for public use.

Dependencies:
- installed version of the BGEN library
- boost (specifically, boost_iostreams to allow read/write to gzipped files)

In the CMakeLists.txt file I currently set paths to these dependencies in lines 9,17,25 (BGEN) and lines 13,21,29 (boost_iostreams).

To compile use:
```
mkdir bin
cd bin
cmake ..
cd ..
cmake --build bin -- -j 4
```

