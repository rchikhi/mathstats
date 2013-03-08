from setuptools import setup, find_packages

setup(
    name='mathstats',
    version='0.1dev',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    classifiers=[
        "Development Status :: 1 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bioinformatics"
      ],
    #scripts=['runBESST'],
    description='Statistical functions, goodness-of-fit tests and special and special distributions not implemented in scipy/numpy .',
    author='Kristoffer Sahlin',
    author_email='kristoffer.sahlin@scilifelab.se',
    url='https://github.com/ksahlin/mathstats',
    license='GPLv3',
    long_description=open('README.md').read(),
    #install_requires=['pysam==0.6',
    #                  'networkx>=1.4'],
    #platforms=['Unix', 'Linux', 'Mac OS']
)
