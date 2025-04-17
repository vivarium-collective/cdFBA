import re
from setuptools import setup, find_packages


VERSION = '0.0.2'


with open("README.md", "r") as readme:
    description = readme.read()
    # Patch the relative links to absolute URLs that will work on PyPI.
    description2 = re.sub(
        r']\(([\w/.-]+\.png)\)',
        r'](https://github.com/vivarium-collective/cdFBA/raw/main/\1)',
        description)
    long_description = re.sub(
        r']\(([\w/.-]+)\)',
        r'](https://github.com/vivarium-collective/cdFBA/blob/main/\1)',
        description2)

setup(
    name="cdFBA",
    version=VERSION,
    author="Tasnif Rahman",
    author_email="trahman@uchc.edu",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    license_files=["LICENSE"],
    url="https://github.com/vivarium-collective/cdFBA",
    # packages=find_packages(),
    packages=[
        'cdFBA',
        'cdFBA.processes',
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.9",
    install_requires=[
        "vivarium-interface",
        "process_bigraph==0.0.33",
        "bigraph-schema==0.0.54",
        "cobra",
        "matplotlib",
        "ipdb",
        "pytest"
    ]
)
