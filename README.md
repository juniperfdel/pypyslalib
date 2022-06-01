# PyPySLALib

### Number converted: 79 /  209

## Description
An attempt to convert the very good [SLALIB](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html) into 
pure python without any wrappers using [pyslalib](https://github.com/scottransom/pyslalib) as a base.

## Workflow
`convert.py` identifies the components of a given Fortran line of code with `pyparsing` then attempts to 
convert it into a line of python code and then finally uses `black` to attempt to format the entire transformed 
function. While this process is not perfect, it is close enough that a human should be able to make it into python with
little effort.


## Legal
The original `SLALIB` was released with [GNUv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html); while
[pyslalib](https://github.com/scottransom/pyslalib) doesn't have a LICENSE file. Based on the conditions 
in the GNUv2 I am assuming `pyslalib` was released under GNUv2 as well, and I am releasing this code under 
GNUv3. The various `LICENSE` files are for each set of code being used and manipulated. I also took the liberty of 
updating the `SLALIB` license to the modern GNUv2 as opposed to the very outdated version (It had the wrong address).
