# Usage

The software has lately changed from bad S3 to probably worse S4.
On the positive site you can now install this software like that:

library(devtools)

install_github('stela2502/StefansExpressionSet')

or

install_github('stela2502/StefansExpressionSet', build_vignettes=FALSE)

# Install on minimal server

I had extreme problems installing this lib on a headless shiny server as it depends on the rgl package.

I always missed was the libpng-devel package!!