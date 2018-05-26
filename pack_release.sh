#!/bin/bash

# This inserts a line after every line with the pattern:
# sed '/Sysadmin/a \ Linux Scripting' finename.txt

timestamp=`date +%Y-%b-%d`
targetdir=QMClecture.$timestamp
mkdir $targetdir

base="
H2_input.f90
paraHe_input.f90
reblocking.f90
orthoHe_input.f90
qmc.F90
rnd.f90
qmc_run.F90
types_const.f90
LICENSE
README
makefile
makefile.sys
"

doc="
man.pdf
man.tex
man.bbl
"

examples="
qmc_run.cfg
"

cfparser="
cfparser.f90
list_template.F90.inc
list_words.F90
makefile
"

for file in $base; do
  cp -p $file $targetdir/.
done

# compress fonts in the PDF manual
cd doc
mv man.pdf __XXmanXX.pdf
gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=man.pdf __XXmanXX.pdf
rm __XXmanXX.pdf
cd ..

mkdir $targetdir/doc
for file in $doc; do
  cp -p doc/$file $targetdir/doc/.
done

mkdir $targetdir/cfg_example
for file in $examples; do
  cp -p cfg_example/$file $targetdir/cfg_example/.
done

mkdir $targetdir/cfparser
for file in $cfparser; do
  cp -p cfparser/$file $targetdir/cfparser/.
done

cp -p makefile.release $targetdir/makefile
cp -pr make $targetdir/.

cd $targetdir
if [ `which hg 2>/dev/null | wc -w` -eq 1 ]; then
  hg tip | grep changeset > VERSION
fi
cd ..
