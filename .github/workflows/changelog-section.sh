#!/bin/bash

# Prints the section of CHANGELOG.md for the given version (eg 1.5.0). Fails if it is not the latest version in the
# changelog, or if the section is empty. Used to check the changelog on release tags, and for the GitHub release notes.

set -e

version=$1
changelog=${CHANGELOG:-$(dirname "$0")/../../CHANGELOG.md}

latest=$(grep -m1 -oP '^## \[\K[0-9][^]]*' "$changelog")
if [ "$latest" != "$version" ]; then
	echo "CHANGELOG.md: the latest version is '$latest', expected '$version'" >&2
	exit 1
fi

# lines after the "## [version]" heading, up to the next "## " heading (or the link references at the end)
section=$(awk -v h="## [$version]" '
	index($0, h) == 1 { found = 1; next }
	found && (/^## / || /^\[[^]]+\]: /) { exit }
	found { print }
' "$changelog" | sed -e '/./,$!d' -e :a -e '/^\n*$/{$d;N;ba' -e '}')		# trim leading/trailing empty lines

if [ -z "$section" ]; then
	echo "CHANGELOG.md: the section of version $version is empty" >&2
	exit 1
fi

echo "$section"
