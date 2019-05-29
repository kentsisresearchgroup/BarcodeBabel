from setuptools import setup, find_packages

setup(
    name='BarcodeBabel',
    version='0.0',
    description='Generate peptide barcodes',
    author='Nicole McNeer, Josh Fass, Alex Kentsis, John Chodera',
    author_email='mcneern@mskcc.org',
    packages=find_packages(),
    url='https://github.com/kentsisresearchgroup/BarcodeBabel',
    license='AGPL-3.0',
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3'],
)