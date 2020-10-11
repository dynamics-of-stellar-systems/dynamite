import setuptools

# this loads the version number from the dynamite/version.py module
version = open("dynamite/_version.py")
version = version.readlines()[-1].split()[-1].strip("\"'")

# load the readme as long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# load the package requirements from requirements.txt
with open("requirements.txt", "r") as fp:
    required = fp.read().splitlines()

setuptools.setup(
    name="dynamite",
    version=version,
    author="Prashin Jethwa, Sabine Thater, Thomas Maindl",
    author_email="prashin.jethwa@univie.ac.at",
    description="dynamics, age and metallicity indicators tracing evolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.univie.ac.at/dynamics/dynamite_docs/index.html",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    project_urls={
        "Source": "https://github.com/dynamics-of-stellar-systems/dynamite/",
        "Documentation": "https://www.univie.ac.at/dynamics/dynamite_docs/index.html",
    },
    python_requires=">=3.6",
    # use the already parsed requirements from requirements.txt
    install_requires=required,
    # extra requirements for testing
    extras_require={
        "testing": [
            "pytest",
            "coverage",
        ]
    },
    package_data={
        "dynamite": [
            "../legacy_fortran/modelgen",
            "../legacy_fortran/orbitstart",
            "../legacy_fortran/orblib",
            "../legacy_fortran/partgen",
            "../legacy_fortran/triaxmass",
            "../legacy_fortran/triaxmassbin",
            "../legacy_fortran/triaxnnls_CRcut",
            "../plotbin/*",
        ]
    },
)
