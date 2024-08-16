workflow:

- verify desired region is correctly mapped by running `plots/BasinGeo.py`
- make sure `EOF.csv` is downloaded from [here]() and in the same folder as `plots/BasinAnalysis.ipynb`
- go through the cells of `BasinAnalysis.ipynb` to subsample data to desired region
- go through the cells `data/merge_eof_slope.ipynb` to generate a csv containing data for t/s/bathy/ssh for your region (slowest part)
    - Southern pacific data is available [here]() and does not have to be generated
- go through the cells of `train.ipynb` to train a model on generated data (connect to a GPU!)