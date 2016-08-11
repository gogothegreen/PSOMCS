PSOMCS
======

A tool for calculating cMCS optimizing multiple reactions in a metabolic network using mixed integer linear programming (MILP) with particle swarm optimization.

Table of Contents
=================
1. [Installation](#Installation)
2. [Examples](#Examples)
3. [The program](#the program)
   * psomcs.pl
   * fba.pl
4. Author
5. Copyright and License

# <a name="Installation"></a>Installation
Scripts to run the program are in the folder **_scripts_**. Running the program requires Christian Jungreuthmayer's [Math::CPLEX](https://homepage.boku.ac.at/jungreuthc/) perl module to be installed. It's also necessary to have the [IBM ILOG CPLEX Optimization Studio](http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud) installed on the system running this program. The scripts can be executed without compilation. Use the -h option for the help page.

# <a name="Examples"></a>Examples
A relatively small example using the *E. coli* core network (Trinh 2008) modified to grow anaerobically on glucose. The relevant files are in the folder **_examples_**. The output files will also be written to this folder.

# <a name="the program"></a>The program

**psomcs.pl**

```
This script runs the PSO algorithm and solves the MILP using CPLEX. The number of threads to be used can be specified using the 'threads' parameter, default is the number of cores - 1.
```

**fba.pl**

```
A modified FBA calculator to get the mimimum and maximum fluxes for the objective reactions while the uptake reactions are limited by the 'normalization_constraints' parameter.
```

# <a name="author"></a>Author
Govind Nair [govind.nair@boku.ac.at](mailto:govind.nair@boku.ac.at)

# <a name="copyright and license"></a>Copyright and License

MIT License

Copyright (c) 2016 Govind Nair and Christian Jungreuthmayer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
