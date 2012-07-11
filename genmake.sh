#!/bin/bash
# ==========================================================================
# A script to generate makefile with module and #include dependencies for
# a Fortran 90/95 project. Perhaps one day there will be some 'dependmaker'
# inside the compilers themselves like 'gcc -MM'.
#
# There are many alternatives to this script in the web, for some reason
# I didn't find them soon enough. In no particular order:
#   http://www.geos.ed.ac.uk/homes/hcp/fmkmf
#   http://www.fortran.com/makemake.perl
#   http://personal.inet.fi/private/erikedelmann/makedepf90
#   http://www.gfdl.gov/~vb/mkmf.html
# --------------------------------------------------------------------------
# Copyright (c) J. Kolorenc <kolorenc@fzu.cz>, 2009
# ==========================================================================

# list of all source files to be considered
suffixes="f90 F90 f F"

# directories where sources reside (relative path, no trailing slash)
srcdirs="."

# directory to put executables (relative path, no trailing slash)
bindir="."

# when searching for modules in source files, I run these files through
# preprocessor first (because of my "templates"); define FPP=cat if the
# preprocessing is not desired
#FPP=fpp
FPP="cpp -traditional-cpp -Wno-endif-labels"

# commands to link an executable (need to start with TAB)
progline="\t@echo Linking \$@\n"
progline=$progline"\t@\$(FC) \$(LDFLAGS) -o \$@ \$< \$(OBJ) \$(LIBS)"

# strip suffix but keep the directory part of the filename
function strip_suffix {
  base=$1
  for suffix in $suffixes; do
    if [[ "$base" == "$1" ]]; then
      base=`dirname $1`"/"`basename $1 .$suffix`
    fi
  done
  echo $base
} # function strip_suffix

bin=""
obj=""
allobj=""
modfiles=""
deps=""
progs=""
module_file_pairs=""

# allow empty $bindir and $srcdir to represent the current directory
if [[ "$bindir" == "" ]]; then
  bindir="."
fi
if [[ "$srcdirs" == "" ]]; then
  srcdirs="."
fi

# combine $srcdirs and $suffixes into a set of masks
masks=""
for suffix in $suffixes; do
  for srcdir in $srcdirs; do
    masks=$masks" "$srcdir"/*."$suffix
  done
done

for f in $masks; do

  if [ -f $f ]; then
    
    base=`strip_suffix $f`
 
    # lines that start with 'use' or that have any number of spaces before it 
    modules=`$FPP $f| grep -i '^ *use ' | \
    # there can be several 'use' statements per line separated by semicolons
      awk -F';' '{ for (i=1; i<=NF; i++) print $i }' | \
    # after each module there can be a coma and 'only:' statement etc.
      awk -F',' '{ print $1 }' | \
    # finally get the module names and leave only unique ones
      awk '{ print $2 }' | sort -u`

    dep=""
    if [[ "$modules" != "" ]]; then
      for mod in $modules; do
	
	# relaxing the above requirement by searching for the module in
        # all source files; $module_file_pairs is in quotes in order the
	# leadind space to be printed (otherwise the first record is not
	# found)
	module_obj=`echo "$module_file_pairs" | sed "s/ $mod /:/g" | \
           awk -F':' '{ print $2 }' | awk '{ print $1 }'`
	if [[ "$module_obj" == "" ]]; then
	  for modf in $masks; do
	    if [ -f $modf ]; then
              modf_base=`strip_suffix $modf`
              modf_test=`$FPP $modf|grep -i -c -w "^ *module *$mod"`
              if [[ "$modf_test" -gt "0" ]]; then
		echo module \"$mod\" found in file \"$modf\" >&2
                module_file_pairs=$module_file_pairs" "$mod" "$modf_base.o
		module_obj=$modf_base.o
		# if not instructed otherwise, compilers put module files in
		# the directory where they are called from. Therefore no
		# directory is added here
		modfiles=$modfiles" "$mod.mod
		break
              fi
            fi # file $modf exists
          done
        fi
	if [[ "$module_obj" != "" ]]; then
	  # don't add the dependency if we already have it
          repetition_test=`echo $dep|grep -c $module_obj`
          if [[ "$repetition_test" -eq 0 ]]; then
	    # don't add dependency on itself (would happen when multiple
	    # modules defined in a single file)
	    if [[ "$module_obj" != "$base.o" ]]; then
	      dep=$dep" "$module_obj
            fi
          fi
        else
          unres=$unres" "$mod
        fi

      done # for mod in $modules
    fi

    # dependencies through '#include' statements
    includes_nodir=`grep "#include" $f | awk '{ print $2 }'|sed "s/\"//g"`

    # prepend includes with the directory where $f resides
    dir=`dirname $f`
    includes=""
    for inc in $includes_nodir; do
      includes=$includes" "$dir/$inc
    done
    
    if [[ "$includes" != "" ]]; then
      for inc in $includes; do
	  repetition_test=`echo $dep|grep -c $inc`
	  if [[ "$repetition_test" -eq 0 ]]; then
            dep=$dep" "$inc
          fi
      done
    fi

    if [[ "$dep" != "" ]]; then
      deps=$deps$base".o:"$dep"\n"
    fi 

    prog=`grep program $f | awk '{ print $1, $2 }' | awk '/^program/' | \
      awk '{ print $1 }'`
    if [[ "$prog" == "program" ]]; then
      # executables will be put into $bindir
	binbase=$bindir/`basename $base`  
        bin=$bin" "$binbase
	progs=$progs$binbase": "$base".o \$(OBJ)\n$progline\n\n"
    else
      obj=$obj" "$base.o
    fi
    allobj=$allobj" "$base.o

  fi  # file $f exists

done

if [[ "$unres" != "" ]]; then
  echo >&2
  echo unresolved modules: >&2
  echo $unres >&2
  echo >&2
fi

tab=`printf "%b" "\t"`

cat << EOF
FC =ifort
OPT=-O2 -pc80 -ip -openmp -openmp-report=2 #-reentrancy threaded
# -stand f03 , -stand f95 , -stand f90 ...to check compliance with standards
# note that -warn by default also makes those stupid files *__genmod.f90
FFLAGS  =\$(OPT) -warn -stand f03
F77FLAGS=\$(OPT) -warn nousage  # no warnings about obsolete Fortran constructs
LDFLAGS =-openmp #-reentrancy threaded
LIBS=

BIN      :=$bin
OBJ      :=$obj
ALLOBJ   :=$allobj
MODFILES :=$modfiles

all: bin

bin: \$(BIN)

`printf "%b" "$progs"`

`printf "%b" "$deps"`

%.o: %.F90
$tab\$(FC) \$(FFLAGS) -c \$< -o \$@

%.o: %.f90
$tab\$(FC) \$(FFLAGS) -c \$< -o \$@

%.o: %.F
$tab\$(FC) \$(F77FLAGS) -c \$< -o \$@

%.o: %.f
$tab\$(FC) \$(F77FLAGS) -c \$< -o \$@

clean:
${tab}rm -f \$(ALLOBJ)

mrproper: clean
${tab}rm -f \$(BIN) \$(MODFILES)

EOF
