#!/bin/sh
# int-or-string.sh

clear
echo "Creating a new version of the LIBQIF Library"
echo "Insert the version number please" 
read version

mkdir LIBQIF-v$version

cp -R html/ LIBQIF-v$version/docs/
cp -R inc/ LIBQIF-v$version/inc/
cp lib/libqif.a LIBQIF-v$version/libqif.a
cp licence.txt LIBQIF-v$version/licence.txt
cp -R samples/ LIBQIF-v$version/samples/

touch LIBQIF-v$version/Readme.txt
echo "This version was generated on" >> LIBQIF-v$version/Readme.txt
date >> LIBQIF-v$version/Readme.txt
echo "" >> LIBQIF-v$version/Readme.txt
cat lib_header.txt >> LIBQIF-v$version/Readme.txt

echo "DATE: " >> versions.txt
date >> versions.txt
echo "Version number: " >> versions.txt
echo $version >> versions.txt
echo "Features: " >> versions.txt
echo "Reported bugs: " >> versions.txt
echo "" >> versions.txt

tar czf LIBQIF-v$version.tar.gz LIBQIF-v$version
sudo rm -r LIBQIF-v$version/
mv LIBQIF-v$version.tar.gz bin/LIBQIF-v$version.tar.gz

echo The new version is in the file bin/LIBQIF-v$version.tar.gz