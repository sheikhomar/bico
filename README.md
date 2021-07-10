# BICO

Code from the [BICO website](https://ls2-www.cs.tu-dortmund.de/grav/en/bico#references).

## Getting Started

Use the pre-made `Makefile` in the `bico/build` directory to build the project:

```bash
make -C bico/build
```

## Downloading Datasets

- US Census Data (1990)

    ```bash
    mkdir data
    curl https://archive.ics.uci.edu/ml/machine-learning-databases/census1990-mld/USCensus1990.data.txt --output data/USCensus1990.data.txt
    ```
