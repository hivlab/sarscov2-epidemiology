# Nextstrain build for Sars-Cov-2

This repository provides the data and scripts for Sars-CoV-2 Nextstrain build and **visualisation**.

## Install

To view Sars-CoV-2 data visualisations in your computer you need to install **auspice**.

To get **auspice** up and running follow these steps:

- If you don't have miniconda already in your computer, install miniconda3, following instructions on miniconda website at <https://docs.conda.io/en/latest/miniconda.html>, otherwise proceed to next step.

- In terminal, create conda environment.

```bash
conda create -n next2 python=3
```

- Enter to conda environment.

```bash
conda activate next2
```

- Install `nodejs` into conda environment.

```bash
conda install -c conda-forge nodejs=10
```

- Install **auspice** server.

```bash
npm install --global auspice
```

> Now you should have software environment ready.

- Next, you need data, to get this dataset, (fork and) clone this repository to your computer.

For example, let's clone this repo to `Downloads` folder in your home directory, and then enter this cloned folder.

For the next step (Running) you need to be inside sars-cov2-est folder.

```bash
cd Downloads
git clone https://github.com/avilab/sars-cov2-est.git
cd sars-cov2-est
```

## Running

To run visualisations, either type and enter:

```bash
npm start
```

OR

Run auspice on your dataset directly:

```bash
auspice view --datasetDir auspice
```

Now, if everything goes as expected, you should see instructions in your terminal that direct you to web browser where you find Sars-CoV-2 dataset and **start exploring**.

## Build

To rebuild this dataset, in case you append new sequences,
you need also install Augur and dependencies as described in <https://nextstrain.org/docs/getting-started/local-installation#install-augur--auspice-with-conda-recommended>.

