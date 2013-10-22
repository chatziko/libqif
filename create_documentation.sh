#!/bin/sh

clear
# Documentation --------------------------------------------------------------------
echo "Creating Library Documentation"

# preparing doxygen files
cd doxygen/
mv Doxyfile footer.html init stylesheet.css mainpage.dox header.html ../
cd ..
#requires doxygen
doxygen Doxyfile

# restoring the doxygen files
mv Doxyfile footer.html init stylesheet.css mainpage.dox header.html doxygen/
