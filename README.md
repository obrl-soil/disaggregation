## modified functions from R DSMART

These scripts are modified versions of the core functionality in the dsmart package for R. All props to Nathan Odgers who [originated the concept](http://dx.doi.org/10.1016/j.geoderma.2013.09.024) and coded it in python, and Brendan Malone who ported it to R . If you want the original python code, hit @nathanodgers up on twitter. If you want the official R package, contact Brendan at brendan.malone@sydney.edu.au. That code is under a GPL-2 license, which this material inherits.

These modifications started with a desire for area-proportional sampling within polygons rather than a flat rate, and then I couldn't leave well enough alone. Along the way I've managed to improve the main function's stability, if not its speed, and I've made the outputs a little easier to work with. As I write this, its happily chugging through 6000 polygons and a 14-raster covariate stack with about 1.4M cells per layer.

### main changes

* area-proportional sampling, with a minimum rate that still allows minor soil classes a look-in
* improved memory management
* all inputs start in sp or raster objects and stay there
* outputs handle NA values correctly e.g. masked out coastal areas
* option of outputs just as rds objects
