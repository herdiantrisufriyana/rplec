# rplec: An R package of placental epigenetic clock to estimate aging by DNA-methylation-based gestational age

`rplec` is an R package designed to estimate placental aging based on gestational age using DNA methylation levels, so called placental epigenetic clock (PlEC). We developed a PlEC for the 2024 Placental Clock DREAM Challenge (https://www.synapse.org/Synapse:syn59520082/wiki/628063). Our PlEC achieved the top performance based on an independent test set. PlEC can be used to identify accelerated/decelerated aging of placenta for understanding placental dysfunction-related conditions, e.g., great obstetrical syndromes including preeclampsia, fetal growth restriction, preterm labor, preterm premature rupture of the membranes, late spontaneous abortion, and placental abruption.


## Features

- **Normalize DNA methylation values**: We provided normalization feature based on beta mixture quantile (BMIQ) method before using a PlEC.

- **Estimate DNA-methylation-based gestational age**: A PlEC is available for estimating gestational age.

- **Perform quality control**: A user can evaluate the PlEC accuracy based on calibration plot, root mean squared error (RMSE), mean absolute error (MAE), and correlation coefficient (Pearson's r) before interpreting the DNA-methylation-based gestational age to identify placental aging.

- **Identify placental aging**: Compare placental aging based on DNA-methylation-based gestational age between condition of interest and control.


## Installation

You can install `rplec` from CRAN with:

```r
install.packages("rplec")
```

You can install the development version of `rplec` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("herdiantrisufriyana/rplec")
```


## Quick Start

Load necessary packages.

```r
library(rplec)
```

Load our example data.

```r
beta_values_case <- load_beta_values_case()
```

Normalize DNA methylation values.

```r
norm_beta_values_case <- bmiq_norm_450k(beta_values_case)
```

Estimate DNA-methylation-based gestational age.

```r
dnam_ga_case <- plec(norm_beta_values_case)
```


## Vignettes

Explore detailed examples and methodologies in the following vignettes:

- [**Placental Aging Analysis**](https://herdiantrisufriyana.github.io/rplec/doc/placental_aging_analysis.html): A data analysis pipeline to identify placental aging using `rplec`.

- [**Reference Manual**](https://github.com/herdiantrisufriyana/rplec/blob/master/extras/rplec_0.1.2.pdf): Comprehensive documentation of all functions and features available in `rplec`. Ideal for detailed reference and advanced use cases.


## License

`rplec` is licensed under the MIT license. See the LICENSE file for more details.


# Citation

If you use `rplec` in your research, please consider citing it:

```bibtex
@misc{rplec2025,
  author = {Herdiantri Sufriyana and Emily Chia-Yu Su},
  title = {rplec: An R package of placental epigenetic clock to estimate aging by DNA-methylation-based gestational age},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\\url{https://github.com/herdiantrisufriyana/rplec}}
}
```


## Contact

For questions or support, please contact herdi[at]nycu.edu.tw.


# Programming Environment

## System requirements

Install Docker desktop once in your machine. Start the service every time you build this project image or run the container.

## Installation guide

Build the project image once for a new machine (currently support AMD64 and ARM64).

```{bash}
docker build -t rplec --load .
```

Run the container every time you start working on the project. Change left-side port numbers for either Rstudio or Jupyter lab if any of them is already used by other applications.

In terminal:

```{bash}
docker run -d -p 8787:8787 -p 8888:8888 -v "$(pwd)":/home/rstudio/project --name rplec_container rplec
```

In command prompt:

```{bash}
docker run -d -p 8787:8787 -p 8888:8888 -v "%cd%":/home/rstudio/project --name rplec_container rplec
```

## Instructions for use

### Rstudio

Change port number in the link, accordingly, if it is already used by other applications.

Visit http://localhost:8787.
Username: rstudio
Password: 1234

Your working directory is ~/project.

### Jupyter lab

Use terminal/command prompt to run the container terminal.

```{bash}
docker exec -it rplec_container bash
```

In the container terminal, run jupyter lab using this line of codes.

```{bash}
jupyter-lab --ip=0.0.0.0 --no-browser --allow-root
```

Click a link in the results to open jupyter lab in a browser. Change port number in the link, accordingly, if it is already used by other applications.






