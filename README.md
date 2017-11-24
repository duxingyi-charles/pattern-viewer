# Pattern Viewer
This is a viewer for the quadrangulation patterns in the database proposed in the publication

> **Data-Driven Interactive Quadrangulation**<br/>
> Giorgio Marcias, Kenshi Takayama, Nico Pietroni, Daniele Panozzo, Olga Sorkine-Hornung, Enrico Puppo, Paolo Cignoni<br/>
> In *ACM Transactions on Graphics (Proceedings of SIGGRAPH 2015)*<br/>
> [Project page](http://igl.ethz.ch/projects/ddq/)

Patterns can be browsed by the number of sides and their index in the database. A simplified pattern is shown if further polychord collapses are possible for the original pattern. Polychord expansions are also demonstrated by changing the width of independent polychords.

This program is mainly based on the Sketch Retopo software. You need to agree to the [software License](http://igl.ethz.ch/projects/sketch-retopo/sketch-retopo-license.html) before using this program.

# Screenshot
![screenshot](https://github.com/dohoney/pattern-viewer/raw/master/resources/screenshot.png)


# Usage
First download the pattern database from the [official site](http://igl.ethz.ch/projects/ddq/).

Install the prerequisites ``qt4`` and ``sqlite3`` and then compile using CMake.

On MacOS, compile with the following commands

	git clone --recursive https://github.com/dohoney/pattern-viewer.git
	cd pattern-viewer
	mkdir build
	cd build
	cmake ..
	make -j 4

Launch the binary with a parameter specifying the path of the pattern database, for example

	./PatternViewer ~/patches.db


