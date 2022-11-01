from setuptools import setup

setup(
    name='exoTras',
    version='0.2',
    author='Ruiqiao He',
    author_email='ruiqiaohe@gmail.com',
    packages=['exoTras'],
    license="GPL",
    url='http://pypi.python.org/pypi/exoTras/',
    description='exosome-containing droplets identification and source tracking in scRNA-seq data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        "scanpy==1.6.0",
        "numpy==1.21.5",
        "pandas==1.1.2",
        "scipy==1.5.4",
        "statsmodels==0.12.1",
        "copy"
        "sys",
        "os",
        "multiprocessing",
        "pickle==4.0",
        "gseapy>=0.10.4",
    ],
    python_requires='==3.7',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
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
