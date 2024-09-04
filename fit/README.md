### Raw yield extraction

Careful: the current `extract_raw_yield.py` macro uses the version of `flarefly`'s dev branch (where `plot_mass_fit()` method returns a tuple and not a single `matplotlib.figure.Figure`).

One can clone the git repo of flarefly and then run the `pip install -e path_to_local_flarefly_git_repo` or directly install it from the online branch location using:
```
pip install git+git@github.com:flarefly/flarefly.git@dev
```