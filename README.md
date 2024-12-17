# boostDM manuscript analyses

Source code to reproduce the figures of the paper:

> **In silico saturation mutagenesis of cancer genes**<br> 
  Ferran Muiños, Francisco Martinez-Jimenez, Oriol Pich, Abel Gonzalez-Perez, Nuria Lopez-Bigas<br>
  DOI: https://doi.org/10.1038/s41586-021-03771-1

## Content

This repo contains the source code to reproduce the main and extended figures of the paper.<br>
Each figure has its own jupyter notebook to render the figure's panels.<br>

#### Main figures
Figure 1: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Figure1/display_panels_Figure1.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Figure1.pdf)]<br>
Figure 2: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Figure2/display_panels_Figure2.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Figure2.pdf)]<br>
Figure 3: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Figure3/display_panels_Figure3.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Figure3.pdf)]<br>
Figure 4: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Figure4/display_panels_Figure4.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Figure4.pdf)]<br>

#### Extended Figures

Extended Figure 1: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_1/display_panels_Extended_Figure_1.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure1.pdf)]<br>
Extended Figure 2: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_2/display_panels_Extended_Figure_2.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure2.pdf)]<br>
Extended Figure 3: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_3/display_panels_Extended_Figure_3.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure3.pdf)]<br>
Extended Figure 4: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_4/display_panels_Extended_Figure_4.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure4.pdf)]<br>
Extended Figure 5: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_5/display_panels_Extended_Figure_5.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure5.pdf)]<br>
Extended Figure 6: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_6/display_panels_Extended_Figure_6.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure6.pdf)]<br>
Extended Figure 7: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_7/display_panels_Extended_Figure_7.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure7.pdf)]<br>
Extended Figure 8: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_8/display_panels_Extended_Figure_8.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure8.pdf)]<br>
Extended Figure 9: [[ipynb](https://github.com/bbglab/boostdm-analyses/blob/master/Extended_Figure_9/display_panels_Extended_Figure_9.ipynb)] [[pdf](https://github.com/bbglab/boostdm-analyses/blob/master/figures_paper/Extended_Figure9.pdf)]<br>

## Complementary content

You can access to boostDM source code and documentation in the [boostDM pipeline repository](https://github.com/bbglab/boostdm-pipeline).<br>
You can explore and download the main outputs of boostDM in the [boostDM website](https://www.intogen.org/boostdm).<br>

## Requirements

#### Download source data

All the code features in this repo feeds on source data. 

Make sure that you download a stable copy of the source data from zenodo and keep it in the root of the repo
from [zenodo](https://zenodo.org/) as follows:

```
$ pip install zenodo_get
$ bash get.sh
$ tar -xvf source-data/source-data-zenodo.tar.gz
$ cp -r source-data/boostdm-analyses .
```

#### Run notebooks with singularity

The notebooks must be run on a jupyter-notebook or jupyter-lab session launched from 
[Singularity](https://sylabs.io/) image that already satisfies all the dependencies for the notebooks to run.<br>

Follow these steps:

* [Install](https://sylabs.io/guides/3.0/user-guide/installation.html#) the latest Singularity release<br>

* Create a singularity image using the [Singularity](https://github.com/bbglab/boostdm-analyses/blob/master/Singularity) recipe:

```
$ sudo singularity build boostdm-analyses.sif Singularity
```

* Now you can run the notebooks from singularity:

```
$ singularity exec boostdm-analyses.sif jupyter-lab
```
