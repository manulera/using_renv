# Using renv, a dependency manager for R

## What is a dependency manager and why it's good

A dependency manager is a program that:
* Documents the exact libraries/dependencies your project uses as you install them.
* Allows someone else (or yourself in another machine) to re-create the same environment (same version of all packages). That makes code execution more deterministic and reproducible.

## How to use renv

Full instructions in [this website](https://rstudio.github.io/renv/index.html) and summarised below:

### Setting up the project

* Go to the folder you want to create your project in.
* In the R console:
    ```R
    # install rev (if you haven't)
    install.packages("renv")

    # run this command on the R console, you will have to type y to accept
    renv::init()

    ```
* This will create the following new files in current directory, described [here](https://rstudio.github.io/renv/articles/renv.html#getting-started).
    ```
    ├── .Rprofile
    ├── .Rproj.user
    ├── renv
    └── renv.lock
    ```

### Add new dependencies

🚨**Very important**🚨. Before you install anything:
* By default `renv` will not update `renv.lock` (the file that contains the info of installed dependencies). Normally, you have to call `renv::snapshot()` to update it. However, you can make it update the file by default (desirable behaviour) by setting the following option:
    ```R
    options(renv.config.auto.snapshot = TRUE)
    ```
* Notice that by defaulr `renv` will only list packages that are imported by the code. If you install a package before you include `library(package-name)` in the code, it will not include it in the lock file (see [docs](https://rstudio.github.io/renv/articles/faq.html#why-isnt-my-package-being-snapshotted-into-the-lockfile)). You can make `renv`` include all the dependencies you install by running:
    ```R
    renv::settings$snapshot.type("all")
    ```

To add a dependency, you can install a package normally, or using the `renv:install` command.

```R
# Here I specify that I want the package to be installed from binary if possible.
renv::install("ggplot2", type = "binary")
```

### Re-create the same environment somewhere else

Now, when you share your project with someone else, they can install the exact version of the packages you had by using the command below:

```R
renv::restore()
```

You can try this in your computer by cloning the repository

```bash
git clone https://github.com/manulera/using_renv
cd using_renv
```

Then in the R terminal:

```R
install.packages("renv")
renv::restore()
```

You can verify that it works by running the `dummy.R` script.

### More details

You can install packages from Bioconductor, github and others, see [full documentation](https://rstudio.github.io/renv/reference/install.html#bioconductor). Below an example to install from bioconductor

```R
renv::install("bioc::S4Vectors")
```

