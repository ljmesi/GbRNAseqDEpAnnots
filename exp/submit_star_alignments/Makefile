.ONESHELL:
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: make_test_dirs, \
help

## make_test_dirs: Make dirs that are useful for testing parsing sample names 
make_test_dirs:
	cat sample_names.txt \
	| xargs -n1 -p -I% sh -c "mkdir -p %; cd %; touch ReadsPerGene.out.tab"

## help: Show this message
help:
	@grep '^##' ./Makefile
