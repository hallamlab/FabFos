import os, sys
from pathlib import Path
import setuptools

HERE = Path(os.path.realpath(__file__)).parent
NAME = "fabfos".lower()
ENTRY_POINTS =  [
    'ffs = fabfos.cli:main',
]
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open(HERE.joinpath(f"src/{NAME}/version.txt")) as f:
    VERSION = f.read()

if __name__ == "__main__":
    setuptools.setup(
        name=NAME,
        version=VERSION,
        author="Tony Liu, Connor Morgan-Lang, Avery Noonan, Zach Armstrong, and Steven J. Hallam",
        author_email="shallam@mail.ubc.ca",
        description="Analysis pipeline for pooled fosmids",
        long_description=long_description,
        long_description_content_type="text/markdown",
        license_files = ('LICENSE',),
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
                "deinterleave_fastq.sh",
            ],
            # examples
            # "package-name": ["*.txt"],
            # "test_package": ["res/*.txt"],
        },
        entry_points={
            'console_scripts': ENTRY_POINTS,
        },
        python_requires=">=3.11",
        install_requires=[
            "packaging >=21.0",
        ]
    )