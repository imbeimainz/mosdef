# mosdef 1.6.0

* Slight edits to have better behaved `de_volcano()` and `go_volcano()` functions, with filters/looks depending on the adjusted p-values, and a better management of the color values to be used when plotting the different subsets

# mosdef 1.4.0

* `gene_plot()` defaults now to NULL in the `intgroup` parameter, which translates into using the first `colData` item

# mosdef 1.2.0

* New functionality added for general scope tasks, such as feature annotation and summaries on the expression tables
* The `gene_plot()` function now is (again) able to handle the specification of more than one grouping variable

# mosdef 1.0.0

* mosdef is on Bioconductor!

# mosdef 0.99.4

* Added persistent location to the script that generates the data objects

# mosdef 0.99.3

* Working on the size and time constraints for the vignette/package to build and check on the BBS

# mosdef 0.99.2

* Better specification of data objects format and location in the installed folder

# mosdef 0.99.1

* This version contains the newly implemented changes as a response to the Bioconductor review.
  In brief, this includes:
  - renaming the functions (and parameters) to a more consistent style
  - reduction of the dependencies and runtime of checks/tests
  - more details on the exported data objects, and on the outputs (detailed in the vignette)
  - better structure for the vignette with cross-references among sections
  - modularization of some functions to avoid repetitive code
  - implementation of an API which is already framework-agnostic, to later accommodate e.g. edgeR/limma.
    Mainly, this implies the renaming of the `dds` to a more generic `de_container`, whereas the `res_de`
    parameter stays constant.
  
  For a full list of all changes implemented, please refer to the PR on the `mosdef`
  repository https://github.com/imbeimainz/mosdef/pull/11/
  

# mosdef 0.99.0

* Ready for submission to Bioconductor!

# mosdef 0.4.0

* Using `longtests` to reduce checking time, as per Bioc guidelines

# mosdef 0.3.0

* Reordered the functions thematically into a better file naming structure
* Renaming of functions to better mirror their aim - e.g. `topGOtable()` -> 
  `run_topGO()`, while also homogenizing the case

# mosdef 0.2.0

* Vignette contains worked out examples for the whole package
* Added the `de_table_painter()` function

# mosdef 0.1.0

* Changed the examples to macrophage, fixed minor bugs and added some parameters to the buttonifier function

# mosdef 0.0.2 

* Added new features (go_volcano, input_checks, de_volcano). Styled up the package and maybe more

# mosdef 0.0.1

* Added the main functionality

# mosdef 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
