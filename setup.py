f***REMOVED*** setuptools import setup
f***REMOVED*** setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="TS_pipeline",
    version="1.0",
    author="Roman Peshkov",
    author_email="***REMOVED***an.peshkov@kemi.uu.se",
    description=" (Semi)automatic pipeline for TS search in molecular photoswitches",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CrespiLab/TS_pipeline",
    install_requires=[
        'numpy>=2.2.2',
        'scipy>=1.15.1',
        'lmfit>=1.3.2',
        'pandas>=2.2.3',
        'PyQt5',
        'matplotlib>=3.10.0'
    ],
    
    #extras_require={
    #},
    
    packages=find_packages(
		where='.',
		include=[
		'src',
		'templates',
		'tools']),
    package_dir={'':'.'},
  
    ## try out later: to make a command-line script (need def main() in main.py)
    #entry_points={
    #    'console_scripts': [
    #        'tsp=main:main',
    #    ],
    #},
    
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating Sy***REMOVED***m :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Computational Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    python_requires='>=3.1',
    keywords=["science", "transition state", "computational"],
)
