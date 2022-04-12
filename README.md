# PyPySLALib

### Number converted: 68 /  209

## Introduction
An attempt to convert the very good [SLALIB](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html) into 
pure python without any wrappers using [pyslalib](https://github.com/scottransom/pyslalib) as a base.

### Workflow
The workflow I created is using the file `convert.py` to do common string conversions from fortran to python then 
use [pycharm community edition](https://www.jetbrains.com/pycharm/download/) refactoring tools to clean it up and fix
the formatting. While this requires some hand's on work, the `convert.py` reduces the workload substantially.

### Transpiling?
Writing a full transpiler, would probably take 
more work than using string substitution and manual re-writing, as fortran can do some crazy stuff; and while 
it's easy for a human to see what the program is doing, having a program figure it out would be a huge endeavor. 
It runs into the problem of [automating re-writing is slower than just re-writing it with string substitutions 
for some common manipulations.](https://xkcd.com/1319/)

### Doc-Strings?
While [pyslalib](https://github.com/scottransom/pyslalib) does have a docstring converter, I am currently
more focused on getting the functionality up and running. When converting a new function, I am leaving the docstring
(except for the license, which is now in a dedicated `SLALIB_LICENSE` file so it is removed) as is. This is so that 
if docstring conversion does come around, I could borrow some work done by [pyslalib](https://github.com/scottransom/pyslalib).

## Legal
The original `SLALIB` was released with [GNUv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html); while
[pyslalib](https://github.com/scottransom/pyslalib) doesn't have a LICENSE file. Based on the conditions 
in the GNUv2 I am assuming `pyslalib` was released under GNUv2 as well, and I am releasing this code under 
GNUv3. The various `LICENSE` files are for each set of code being used and manipulated. I also took the liberty of 
updating the `SLALIB` license to the modern GNUv2 as opposed to the very outdated version (It had the wrong address).
