LIBQIF
======

For LIBQIF _users_ you can download the library from the bin/ directory

For LIBQIF _contributers_ :
- There is an script for the instalation of the needeed libraries called installing_needeed_libs.sh
- There is an script for compile the LIBQIF library called compile.sh
- There is an script for generating the LIBQIF documentation with doxygen called create_documentation.sh
 There is an script for running the test cases called run_tests.sh
- There is an script for generation a new version of LIBQIF creating a .tar.gz on the bin/ directory.

The repository structure: 
- bin: It contains the .tar.gz files with the versions of LIBQIF
- inc: It contains the headers (.h files or includes)
- src: It contains the LIBQIF implementation files (.cpp files)
- html: It contains the documentation of the library generated automatically with Doxygen. The documentation file is index.html
- doxygen: It contains the configuration files for generating the library documentation
- tests: It contains the unit test cases that g-test uses for the library testing.
- samples: It contains some LIBQIF use examples.
- papers: It contains the papers that are referenced from the library documentation
- lib: It contains all the libraries needeed for developing LIBQIF 
- sci_files: It contains the files for generating the plots with .sci extension.
