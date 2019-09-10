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
javac -classpath $CLASSPATH *.java
