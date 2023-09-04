import os, sys
from pathlib import Path
import setuptools

HERE = Path(os.path.realpath(__file__)).parent
NAME = "fabfos".lower()

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open(HERE.joinpath(f"src/{NAME}/version.txt")) as f:
    version = f.read()

setuptools.setup(
    name=NAME.title(),
    version=version,
    author="Tony Liu, Connor Morgan-Lang, Avery Noonan, and Steven Hallam",
    author_email="shallam@mail.ubc.ca",
    description="Analysis pipeline for pooled fosmids",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=f"https://github.com/hallamlab/{NAME}",
    project_urls={
        "Bug Tracker": f"https://github.com/hallamlab/{NAME}/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    # packages=pks,
    package_data={
        "":[ # "" is all packages
            "version.txt",
            "usearch",
        ],
        # examples
        # "package-name": ["*.txt"],
        # "test_package": ["res/*.txt"],
    },
    entry_points={
        'console_scripts': [
            'ffs = fabfos.cli:main',
        ]
    },
    python_requires=">=3.11",
    install_requires=[
        "pyfastx >=0.8.4",
        "packaging >=21.0",
    ]
)