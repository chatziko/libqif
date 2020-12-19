#!/usr/bin/perl -w
use strict;

use File::Spec;

system qq{
    pip install --upgrade delocate
};

# args are paths to .whl files
for my $path (@ARGV) {
	my ($volume, $dir, $file) = File::Spec->splitpath($path);

	$file =~ /(.*?)-cp/	or die "invalid file";
	my $base = $1;

	system qq{
		cd $dir

		# delocate-wheel does not work with .so in the top dir, so
		# we move .so files under "a" and create __init__.py
		wheel unpack $file
		mkdir $base/a
		touch $base/a/__init__.py
		mv $base/*.so $base/a

		wheel pack $base
		rm -rf $base

        # workaround delocate issue https://github.com/matthew-brett/delocate/issues/62
        ln -s /usr/local/lib \\\@loader_path
        ln -s /usr/local/lib \\\$ORIGIN

		# run delocate, put libs in ".qif.dylibs" instead of ".dylibs" to avoid conflicts with other packages
		delocate-wheel -v -L .qif.dylibs $file

        rm \\\@loader_path \\\$ORIGIN

		# move stuff under "a" to the root level
		wheel unpack $file
		mv $base/a/* $base/a/.qif.dylibs $base
		rmdir $base/a

		wheel pack $base
		rm -rf $base
	};
}