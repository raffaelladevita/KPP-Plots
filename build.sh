#!/bin/csh -f
#=================================================================
# BUILDING SCRIPT for COATJAVA PROJECT (first maven build)
# then the documentatoin is build from the sources and commited
# to the documents page
#=================================================================
# Maven Build

if(`filetest -e lib` == '0') then
    mkdir lib
endif

# ftCalCalib
echo "Building KPP-Plot..."
    rm lib/*
    mvn install
    mvn package
    cp target/KPP-Plots-3.2-jar-with-dependencies.jar lib/


# Finishing touches
echo ""
echo "--> Done building....."
echo ""
echo "    Usage : build.sh"
echo ""
