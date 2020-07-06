# coding: utf-8
from setuptools import setup

setup(
    name="GenSeek",
    version="0.0.0.1",
    description="A python library and tool to generate and explore conformational spaces of materials",
    keywords="Structure Search",
    url="https://github.com/maksimovdmitrii/genseek/",
    author="Dmitrii Maksimov",
    author_email="maksimov@fhi-berlin.mpg.de",
    classifiers=[
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    py_modules=[
        "genseek",
    ],
    entry_points={
        "console_scripts": [
            "genseek = genseek:main",
        ],
    },
    install_requires=[
        "numpy>=1.12.0",
        "scipy>=0.19.1",
        "ase"

    ]
)
