#!/bin/sh

astyle \
	--style=java \
	--indent=tab --indent-classes \
	--pad-oper --unpad-paren \
	--align-pointer=type --align-reference=type \
	-n \
	src/*.cpp inc/*.h tests/*.cpp samples/*.cpp

