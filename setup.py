
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="frag-lib-filter",
    version="0.0.1",
    author="Helena Martin",
    author_email="helena.martin@nostrumbiodiscovery.com",
    description="Fragment library package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hemahecodes/fragment_library",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
