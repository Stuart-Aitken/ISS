#!/bin/sh
DIRLIBS="$PWD/lib/*.jar"
CLASSPATH="."
for i in ${DIRLIBS}
do
    if [ "$i" != "${DIRLIBS}" ] ; then
        CLASSPATH="$i":$CLASSPATH
    fi
done

echo $CLASSPATH
for i in {2..10}
do
    java -Xmx5000m -classpath $CLASSPATH NBO $1 $i $2  # run Java class NBO; usage random_seed [int] no_classes [int] data_file [string]
done
