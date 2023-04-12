from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='altmap',
    packages=['altmap'],
    package_dir={'altmap': './altmap'},
    package_data={'altmap': ['data/*.ndjson','data/*.csv']},

    version='0.0.1',

    license='MIT',

    install_requires=['numpy','tqdm','requests','gpxpy','ndjson','scipy'],
    author='obscraft23',
    author_email='obscraft23@gmail.com',

    url='https://github.com/obscraft23/altmap',

    description='A package for alpine climbers to handle GSI raw data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='GIS GSI map',

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
)
