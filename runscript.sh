ifort -o apex2000_prog apex2000.f
make
./apex2000_prog  < input_date > outfile
./apex_prog < input_date >> outfile
