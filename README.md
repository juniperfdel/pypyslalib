# PyPySLALib

### Number converted: 12 /  209

## Introduction
An attempt to convert the very good [SLALIB](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html) into 
pure python without any wrappers using [pyslalib](https://github.com/scottransom/pyslalib) as a base. The reasons for this
conversion is:
1. gfortran is not very portable and requires some workarounds to get working on Windows
2. Moving the code into a more readable format allows people to understand what specifically it is doing
3. Python has some really nice features in it that make the library simpler
4. As we move into the future, fortran 77 is probably not going to follow us very well.
5. If Numpy is effectively utilized, then the slow-down from moving from fortran is minimized while gaining the benefits.  

### Workflow
The workflow I created is using the file `convert.py` to do common string conversions from fortran to python then 
use [pycharm community edition](https://www.jetbrains.com/pycharm/download/) refactoring tools to clean it up and fix
the formatting. While this requires some hand's on work, the `convert.py` reduces the workload substantially.

### Transpiling?
Writing a full transpiler, while possible (convert Fortan 77 to an AST then back to python code); would probably take 
more work than using simple string substitution and manual re-writing; as some stuff like `GOTO` (and recognizing which
arguments are being used as output) requires a bit of cleverness to properly convert and/or work around. 
It runs into the problem of [automating re-writing is slower than just re-writing it with some string substitutions 
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

## Help
Feel free to help out, it doesn't take much time especially if you are already have a passing familiarity with python,
just ensure the formatting follows everything else and make a pull request.
