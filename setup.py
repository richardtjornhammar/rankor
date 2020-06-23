import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "rankor",
    version = "0.2.0",
    author = "Richard Tjörnhammar",
    author_email = "richard.tjornhammar@gmail.com",
    description = "Rankor Analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/richardtjornhammar/rankor",
    packages = setuptools.find_packages('src'),
    package_dir = {'rankor':'src/rankor','quantification':'src/rankor','reducer':'/src/reducer','contrasts':'/src/contrasts','imputation':'/src/imputation'},
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
