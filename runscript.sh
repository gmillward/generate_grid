#!/bin/bash
# apex2000.f is self-contained (includes all library routines)
gfortran -std=legacy -ffixed-line-length-132 -w -O2 -o apex2000_prog apex2000.f

make

./apex2000_prog  < input_date > outfile
./apex_prog < input_date >> outfile
