from setuptools import setup

setup(
    name='exoTras',
    version='0.1',
    author='Ruiqiao He',
    author_email='ruiqiaohe@gmail.com',
    packages=['exoTras'],
    url='http://pypi.python.org/pypi/exoTras/',
    description='exosome-containing droplets identification and source tracking in scRNA-seq data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        "scanpy",
        "numpy",
        "pandas",
        "scipy",
        "statsmodels.api",
        "copy"
        "sys",
        "os",
        "pickle",
        "multiprocessing",
        "gseapy",
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'exoTras=exoTras.main:exoTras_command',
        ]
    }
)