#!/usr/bin/env bash
{ [ -n "$1" ] && echo "$1" || cat; } | tr 'ACTGactg' 'TGACtgac' | rev
