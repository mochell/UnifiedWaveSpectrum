from setuptools import setup, find_packages

setup(
    name='UnifiedWaveSpectrum',
    version='0.1',
    packages=find_packages(),
    description='A package to compute the unified wave spectra from Elifouhaily et al. (1997)',
    author='Momme Hell',
    author_email='mhell@ucar.edu',
    url='https://github.com/mochell/UnifiedWaveSpectra',  # if you want to link to a GitHub repo
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
    install_requires=[
        'numpy',
        # Add other packages that your module depends on
    ],
)