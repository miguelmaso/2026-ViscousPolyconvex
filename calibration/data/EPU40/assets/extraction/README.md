Reproduce (Python, `pip install pymupdf numpy matplotlib`), from `data/EPU40`:

    python assets/extraction/extract.py  assets/EPU40NonLin2020.pdf   # -> raw.json (strain, MPa, all vertices)
    python assets/extraction/resample.py raw.json .                    # -> quasi-static.csv, cyclic-loading.csv
    python assets/extraction/overlay.py  assets/EPU40NonLin2020.pdf raw.json .   # -> check_overlay.png

The drawing indices of the curves and axes boxes in `extract.py` are specific to this PDF.
Stretch increments are set at the top of `resample.py` (D_QS, D_CYC).
