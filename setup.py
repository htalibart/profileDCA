from setuptools import setup, find_packages
from glob import glob

print(find_packages())

setup(
    name='profileDCA',
    version='1.0',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'biopython>=1.78', 'pandas', 'seaborn'],  
    package_data={
        'profileDCA_align':['ppalign_solver.so'],
    },
    entry_points={
        'console_scripts':['profileDCA_viz = profileDCA_viz.__main__:main', 'profileDCA_build = profileDCA_build.__main__:main', 'profileDCA_align = profileDCA_align.__main__:main', 'fasta2csv = profileDCA_utils.fasta2indices:main', 'profileDCA_3Dviz = profileDCA_3Dviz.vizpymol:main']
        },
)
