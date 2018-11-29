2018-tms-netsci
==============================

Analysis of various tms induced networkss


Data Description
========================

```
| L_Fp | L_IFG_Anat | L_IFG_BLA | L_VmFp | L_aMFG | L_pMFG | R_FEF | R_Fp | R_IFG_Anat | R_IFG_BLA | R_IFJ | R_IPL | R_M1 | R_VmFp | R_aMFG | R_pMFG | R_preSMA |
|:----:|:----------:|:---------:|:------:|:------:|:------:|:-----:|:----:|:----------:|:---------:|:-----:|:-----:|:----:|:------:|:------:|:------:|:--------:|
|  64  |     29     |     29    |   29   |   75   |   78   |   76  |  67  |     27     |     39    |   74  |   49  |  77  |   31   |   78   |   77   |    62    |
```

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── results            <- Any trained models or results of analysis
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── results         <- Scripts used to write experiments, analyse data and generate results
    │   │   │      
    │   │   ├── analysis.py
    │   │   
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.testrun.org


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
