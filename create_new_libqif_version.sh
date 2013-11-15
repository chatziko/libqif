#!/bin/sh
# int-or-string.sh

clear
echo "Creating a new version of the LIBQIF Library"
echo "Insert the version number please" 
read version

mkdir LIBQIF-v$version

cp -R html/ LIBQIF-v$version/docs/
cp -R inc/ LIBQIF-v$version/inc/
cp -R sci_files/ LIBQIF-v$version/sci_files/
cp lib/liblibqif.a LIBQIF-v$version/liblibqif.a
cp licence.txt LIBQIF-v$version/licence.txt
cp -R samples/ LIBQIF-v$version/samples/
mkdir LIBQIF-v$version/thirdparty/
cp -R lib/armadillo-3.920.1/ LIBQIF-v$version/thirdparty/armadillo-3.920.1/
cp -R lib/glpk-4.52/ LIBQIF-v$version/thirdparty/glpk-4.52/

touch LIBQIF-v$version/INSTALL.sh
echo "#!/bin/sh" >> LIBQIF-v$version/INSTALL.sh
echo "# creating the libraries needeed by libqif" >> LIBQIF-v$version/INSTALL.sh
echo "echo Installing armadillo" >> LIBQIF-v$version/INSTALL.sh
echo "cd thirdparty/armadillo-3.920.1/" >> LIBQIF-v$version/INSTALL.sh
echo "./configure" >> LIBQIF-v$version/INSTALL.sh
echo "cmake ." >> LIBQIF-v$version/INSTALL.sh
echo "make" >> LIBQIF-v$version/INSTALL.sh
echo "sudo make install" >> LIBQIF-v$version/INSTALL.sh
echo "cd .." >> LIBQIF-v$version/INSTALL.sh
echo "cd .." >> LIBQIF-v$version/INSTALL.sh

echo "echo Installing glpk" >> LIBQIF-v$version/INSTALL.sh
echo "cd thirdparty/glpk-4.52/" >> LIBQIF-v$version/INSTALL.sh
echo "./configure" >> LIBQIF-v$version/INSTALL.sh
echo "make" >> LIBQIF-v$version/INSTALL.sh
echo "sudo make install" >> LIBQIF-v$version/INSTALL.sh
echo "cd .." >> LIBQIF-v$version/INSTALL.sh
echo "cd .." >> LIBQIF-v$version/INSTALL.sh

touch LIBQIF-v$version/Readme.txt
echo "This version was generated on" >> LIBQIF-v$version/Readme.txt
date >> LIBQIF-v$version/Readme.txt
echo "" >> LIBQIF-v$version/Readme.txt
cat lib_header.txt >> LIBQIF-v$version/Readme.txt

echo "DATE: " >> versions.txt
date >> versions.txt
echo "Version number: "$version >> versions.txt
echo "Features: " >> versions.txt
echo "Reported bugs: " >> versions.txt
echo "" >> versions.txt

tar czf LIBQIF-v$version.tar.gz LIBQIF-v$version
sudo rm -r LIBQIF-v$version/
mv LIBQIF-v$version.tar.gz bin/LIBQIF-v$version-linux-32bits.tar.gz

echo The new version is in the file bin/LIBQIF-v$version-linux-32bits.tar.gz