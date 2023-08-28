# AmorProt

<img src="https://github.com/mhlee216/AmorProt/blob/main/main.png">

#### AmorProt: Amino Acid Molecular Fingerprints Repurposing-based Protein Fingerprint

<a href="https://doi.org/10.1021/acs.biochem.3c00253">https://doi.org/10.1021/acs.biochem.3c00253</a>

Myeonghun Lee<sup>+,</sup>\*, Kyoungmin Min\*

Biochemistry

* PyPI:

```console
$ pip install amorprot
```

* Example:

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
```
