# MolWall (aka Wall of Molecules)

Streamlit app that displays a "wall of molecules" interface 
with possibility to rate molecules.

## Setup
Create environment using conda or [mamba](https://github.com/mamba-org/mamba).
```shell
conda install mamba
mamba env create --file environment.yml
```

There is also `requirements.txt` that we use for [rsconnect-python](https://docs.rstudio.com/rsconnect-python/).
It can be used for pip, virtualenv or Poetry too.

## Running

```shell
conda activate molwallst
streamlit run molwall.py
```

## License
Apache 2.0.