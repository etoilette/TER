# Update:

We are currently cleaning up the code.
It will be available from 20/07 on at the latest.

# pn2HDA

pn2HDA is a tool to convert different types of Petri Nets into (partial) Higher Dimensional Automata.

It can be used as a headed only library or as the stand-alone tool pn2HDA.

## Build

The project depends on the symmetri library https://github.com/thorstink/Symmetri.
Otherwise, to build the project in a folder called "build"
cmake -DCMAKE_BUILD_TYPE="Release" -S . -B build
cmake --build build
This is also what the helper scripts for plotting expect.

## Recreating the examples
The examples from the paper can be translated to ST-automata calling
python3 helpers/show_pHDA.py
This will create a folder called "examples" at the root containing the dot
and pdf files corresponding to each of the examples.

## Recreating the statistics plot
To recreate a plot similar to the one given in the paper showing the 
number of cells as a function of their dimension.
The plot created by default here contains more examples and additional information.
For instance, we mark now also the number of maximal cells, that is cells 
of maximal concurrency that do not appear as the upper or lower face of another cell.
Warning, to run all examples, we recommend at least 16Gb of free ram
and it will take a couple of minutes to generate.
python3 --runstats