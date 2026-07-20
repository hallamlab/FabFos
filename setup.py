import os, sys
from pathlib import Path

HERE = Path(os.path.realpath(__file__)).parent
# src-layout: the package lives in ./src, this file sits at the repo root.
sys.path = [str(p) for p in {HERE / "src"} | set(sys.path)]
import setuptools
from fabfos import NAME, USER, ENTRY_POINTS, __version__, SHORT_SUMMARY

README = HERE / "README.md"
long_description = README.read_text(encoding="utf-8") if README.exists() else SHORT_SUMMARY

if __name__ == "__main__":
    setuptools.setup(
        name=NAME,
        version=__version__,
        author="Tony Liu, Connor Morgan-Lang, Avery Noonan, Zach Armstrong, and Steven J. Hallam",
        author_email="shallam@mail.ubc.ca",
        description=SHORT_SUMMARY,
        long_description=long_description,
        long_description_content_type="text/markdown",
        license_files=("LICENSE",),
        url=f"https://github.com/{USER}/{NAME}",
        classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: Unix",
        ],
        package_dir={"": "src"},
        packages=setuptools.find_packages(where="src"),
        package_data={
            "": [
                "version.txt",
                # the bundled metasmith library shipped with conda installs
                "_library/**/*",
            ],
        },
        include_package_data=True,
        entry_points={"console_scripts": ENTRY_POINTS},
        python_requires=">=3.12",
        # metasmith provides the planner/executor; it is a conda dependency
        # (see conda_recipe), not a pip one.
        install_requires=[],
    )
