from setuptools import setup
###Metagenome Operational Taxanomic Unit (OTU) Data Transformer
### install locally with 'pip install .'
setup(name="motupy",
        version="0.1dev0",
        description='Tools to manipulate and transform microbiome/metagenomics dataset.',
        author='Chuan Fu Yap',
        author_email='yapchuanfu@gmail.com',
        license='MIT',
        packages=["motupy"],
        install_requires=['pandas', 'ete3', 'numpy'])