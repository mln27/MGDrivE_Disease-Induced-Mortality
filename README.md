# MGDrivE - K.M.

### Repository Organization

* [MGDrivE2](./MGDrivE2) - Adapted version of MGDrivE2
* [MGDrivE2_1.1.0.tar.gz](MGDrivE2_1.1.0.tar.gz) - CRAN version from which this adaptation is based
* [Projects](./Projects) - Directory containing scripts for current projects

### Package Summary

This project is a expansion of the basic [MGDrivE](https://cran.r-project.org/package=MGDrivE2) 
package. Upgrades expand significantly the biological functionality and provide 
computational improvements.

15 new inheritance patterns have been added to this package. All of these gene drives 
are versions of the self-excision mechanism ([S.E.M.](https://doi.org/10.1098/rstb.2019.0804)) - 
one potential direction for autonomous removal of a gene-drive system. These inheritance 
patterns represent the state-of-the-art for controllable gene drives. There are 5 
mechanistically-unique designs with 3 computationally-different versions of each. 
Additionally, the possibility for disease-impacts on mortality has been included.

Beyond biological extensions, several software updates have been performed. 
Code duplication has been reduced, function APIs and documentation has been cleaned 
and improved, and novel plotting and analysis routines have been added. Finally, 
algorithm improvements provide a 5-10x improvement in simulation run time. 

### Building this Package

There are several methods for building R packages, I will cover quickly cover 
the three main methods I use. Additionally, here are references from people better 
than me.  

There may be packages missing, specifically `devtools`, and the build may continually 
fail and tell you what packages are required for completion. Install any required 
packages using `install.packages("pkg-name")`.

NB: To use this repository cannot be used in conjunction with the current CRAN 
version of the MGDrivE2 package (make sure to first uninstall the package if 
already installed or not install it even if prompted).
* RStudio Method
  * This is the easiest method, and can be done all in the GUI.
  1. Download this repository.
  2. Navigate into the package directory - this is the `~/MGDrivE2/` directory.
  3. Double click the `MGDrivE2.Rproj` file - this opens the packge in RStudio.
  4. Go the to `Build` tab at the top of the window.
  5. click 'load all'
  6. Click `clean and rebuild` - this will build and install the package, and the 
    basic documentation, but not the vignettes.
* CMDline Method (in R)
  * This is very similar to the RStudio method, and makes use of the `devtools` package.
  1. Download this repository.
  2. In the CMDline, navigate to the package directory  - this is the `~/MGDrivE2/` directory.
  3. Begin R, and type `devtools::document()` - this builds the basic documentation.
  4. When that is finished, type `devtools::build()` - this actually builds the package.
  5. Finally, type `devtools::install()` - this finishes the installation process, 
    everything is ready to go.
* CMDline Method (bash)
  * This is a more difficult method, but does not involve extra packages such as `devtools`
  1. Download this repository.
  2. Navigate into the package directory - this is the `~/MGDrivE2/` directory.
  3. Type `R CMD build ./` - this builds the package and all documentation, including 
    vignettes by default.
  4. Type `cd ../` (alternatively, you can stay 1 level higher in step 2, and replace 
    `./` with `MGDrivE2`)
  5. Type `R CMD INSTALL MGDrivE2.tar.gz` - this finishes the installation.
  
### Further Reading

Further references and good places to learn basics or more advanced concepts.

* [MIT](https://web.mit.edu/insong/www/pdf/rpackage_instructions.pdf)
* [Coding Club](https://ourcodingclub.github.io/tutorials/writing-r-package/)
* [Prestevez](https://www.prestevez.com/post/r-package-tutorial/)
* [kbroman](https://kbroman.org/pkg_primer/pages/build.html)
  * I like everything this guy has blogged.
* [r-pkgs](https://r-pkgs.org/index.html)
  * This one I can vouch for, it's also an entire book. 
* [adv-r](https://adv-r.hadley.nz/)
  * This will take you zero-to-hero, but it's a deep dive.
* [r4ds](https://r4ds.had.co.nz/index.html)
  * R for data science, also includes projects and packages.
* [More Advanced](https://support.rstudio.com/hc/en-us/articles/200486488-Developing-Packages-with-the-RStudio-IDE)
  * Several links to other resources, not the most clear, but very complete.
  * Includes some info for windows vs linux.
* [Video](https://www.youtube.com/watch?v=79s3z0gIuFU)
  * 15-min walkthrough of setting up a new package through RStudio.
* [Slide Deck](https://johnmuschelli.com/smi_2019/index.html#1)
