# PyPySLALib

### Number converted: 82 /  209

## Description
An attempt to convert the very good [SLALIB](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html) into
pure python without any wrappers using [pyslalib](https://github.com/scottransom/pyslalib) as a base.

## Workflow
`convert.py` identifies the components of a given Fortran line of code with `pyparsing` then attempts to
convert it into a line of python code and then finally uses `black` to attempt to format the entire transformed
function. While this process is not perfect, it is close enough that a human should be able to make it into decent
python with little effort and fix any errors.


## Legal
The original `SLALIB` was released with [GPLv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html); and the README in
[pyslalib](https://github.com/scottransom/pyslalib) indicates that it too is released under GPLv2 as well. I am releasing this code under
GPLv3 per the conditions in GPLv2. The various `LICENSE` files are for each set of code being used and manipulated. I also took the liberty of
updating the `SLALIB` license to the modern GNUv2 as opposed to the outdated version with the wrong address.
