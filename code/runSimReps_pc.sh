#! /usr/bin/env sh

NPARAMS=2
if [ "$#" -ne $NPARAMS ]; then
	echo
	echo "Must supply:"
	echo "(1) .properties file name"
	echo
	exit
fi

pwd

# Revised: get the output directory from the properties file
#OUTDIR=`sed -n -e '/^destinationDirectory/p' $1 | cut -d "=" -f 2 | xargs`
OUTDIR=$2
echo
echo "output directory = "$OUTDIR

# Copy the properties file
cp ${1} $OUTDIR
# Copy the site file
SITEFILE=`sed -n -e '/^sitefile/p' $1 | cut -d "=" -f 2 | xargs`
echo "sitefile = "$SITEFILE

cd $OUTDIR
pwd

c:/Users/sa04ts/.jdks/openjdk-17.0.1/bin/javaw -Xmx8192m -Xms4096m -jar ../../../../jar/particle_track.jar ${1} | tee stdout.txt
mkdir connectivity
mv connectivity_* connectivity/

