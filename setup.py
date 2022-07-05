from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="xrlesionfinder",
    version="1.0",
    description='Helper scripts for use in various bioinformatics projects',
    long_description_content_type="text/markdown",
    url='https://github.com/bmorledge-hampton19/xrlesionfinder',
    author='Ben Morledge-Hampton',
    author_email='b.morledge-hampton@wsu.edu',
    license='MIT',
    python_requires='>=3.7',
    packages=find_packages(),
    package_data={"xrlesionfinder": ["alignment_and_formatting/*.bash", "XR-seq_primers.fa",
                                     "alignment_and_formatting/exact_unique_match_custom_args.txt"]},
    entry_points=dict(
        console_scripts=['xrlesionfinder=xrlesionfinder.Main:main']
    ),
    install_requires=["benbiohelpers"]
)