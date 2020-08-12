from setuptools import setup, find_namespace_packages
### Metagenome Operational Taxanomic Unit (OTU) Data Transformer
### install locally with 'pip install .'
setup(name="motupy",
        version="0.1dev1",
        description='Tools to manipulate and transform microbiome/metagenomics dataset.',
        author='Chuan Fu Yap',
        author_email='yapchuanfu@gmail.com',
        license='MIT',
        packages=find_namespace_packages(include=['motupy.dataprocessing', "motupy.utils", "motupy.timeseries"]),
        install_requires=['pandas', 'ete3', 'numpy', 'scipy', 'scikit-bio'])