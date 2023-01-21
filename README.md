# AmorProt: Amino Acid Molecular Fingerprints Repurposing-based Protein Fingerprint

* PyPI (pip):

```console
$ pip install amorprot
```

* Example (python):

```python
from amorprot import AmorProt

ap = AmorProt(maccs=True, ecfp4=True, ecfp6=True, rdkit=False)
fp = ap.fingerprint('MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHL')
```


```python
from amorprot import AmorProt
import pandas as pd
import numpy as np
import parmap

def make_fp(inputs):
    ap, sq = inputs
    return ap.fingerprint(sq).tolist()

df = pd.read_csv('./data/example.csv')
ap = AmorProt(maccs=True, ecfp4=True, ecfp6=True, rdkit=True)
fps_list = parmap.map(make_fp, [[ap, sq] for sq in df['sequence'].tolist()], 
                      pm_pbar=True, pm_processes=20)
fps = np.array(fps_list)
fps
```
