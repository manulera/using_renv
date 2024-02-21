# Using renv, a dependency manager for R

## What is a dependency manager and why it's good

A dependency manager is a program that:
* Documents the exact libraries/dependencies your project uses as you install them.
* Allows someone else (or yourself in another machine) to re-create the same environment (same version of all packages). That makes code execution more deterministic and reproducible.

## How to use renv

Follow the instructions in this website (summarised below): https://rstudio.github.io/renv/index.html

In the R console:

```R
# install rev
install.packages("renv")

# run this command on the R console, you will have to type y to accept
renv::init()

```

This will create the following new files in current directory, described [here](https://rstudio.github.io/renv/articles/renv.html#getting-started).

```
├── .Rprofile
├── .Rproj.user
├── renv
└── renv.lock
```

You can then, from the R console, install a package normally, or using the `renv:install` command.

```
renv::install("ggplot2")
```

Now, when you share your project with someone else, they can install the exact version of the packages you had by calling in the project directory.

```
renv::restore()
```

## More details

You can install packages from Bioconductor, github and others, see [full documentation](https://rstudio.github.io/renv/reference/install.html#bioconductor). Below an example to install from bioconductor

```
renv::install("bioc::S4Vectors")
```

