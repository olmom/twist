<!-- badges: start -->
<!--<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205155"><img src="https://img.shields.io/badge/Data-GSE205155-green.svg?style=plastic" alt="" /></a>-->
<!--[![](https://img.shields.io/badge/Data-10.1101/2022.06.03.494693-blue.svg)](https://doi.org/10.1101/2022.06.03.494693)-->
[![](https://img.shields.io/badge/Preprint-10.1101/2023.02.12.528191-yellow.svg)](https://www.biorxiv.org/content/10.1101/2023.05.17.541139v1)
<!--[![](https://img.shields.io/badge/Data-10.1101/2022.06.03.494693-blue.svg)](https://doi.org/10.1101/2022.06.03.494693)-->
 <!-- badges: end -->

# twist

This repository contains the reproducible code for the generation of data, analysis and figures in the manuscript "Are circadian amplitudes and periods correlated? A new twist in the story" ([preprint](https://www.biorxiv.org/content/10.1101/2023.05.17.541139v1)). The Supplementary Material of the manuscript is provided in this repository. 

To execute this code:

1. Clone this repository at a suitable location. This will place the code within a directory named **twist**
<!--2. Download all the simulated data from [here](https://www.zenodo.org/) (under the `results` folder) (or alternatively generate all the simulated data using the *main.py* script)-->
2. Make sure you have Python 3.7+ and the following libraries installed: `numpy`, `scipy`, `pandas`, `astropy`, `matplotlib`
3. To reproduce the figures, the Python files can now be executed in order *1_harmonic.py*, *2_goodwin.py,* ... within the project 

The Python script used to generate results and figures rely on classes and objects that are stored under `utils/`.

<!--To reproduce the figures, the Python files can now be executed in order *fig1.py, fig2.py,* ... within the project.-->
