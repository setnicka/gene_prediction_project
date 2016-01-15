#!/bin/sh

for file in *.gb; do
	<$file grep "  gene  " | tr -s ' ' | cut -d" " -f3 >$file.genes
done
