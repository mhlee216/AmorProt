from setuptools import setup, find_packages

setup(name='amorprot',
      version="1.0.0",
      url='https://github.com/mhlee216/AmorProt',
      packages=find_packages(),
      author='Myeonghun Lee',
      author_email="leemh216@gmail.com",
      description='AmorProt',
      long_description='AmorProt: Amino Acid Molecular Fingerprints Repurposing-based Protein Fingerprint',
      install_requires=["numpy >= 1.19.0", "rdkit >= 2021.09.2"],
      classifiers=['License :: OSI Approved :: MIT License'])
