from setuptools import setup, find_packages
from setuptools.command.install import install

class CustomInstallCommand(install):
    def run(self):
        # Call the default install method
        install.run(self)
setup(
    name='Azurify',
    version='0.9.9',
    packages=find_packages(),
    author='Ashkan Bigdeli',
    author_email='ashkan.bigdeli@pennmedicine.upenn.edu',
    description='Azurify is a ML based classifier that predicats small variant pathogencity',
    long_description='This is an example Python package with a setup.py file.',
    url='https://github.com/yourusername/example_package',
    license='GNU',
    classifiers=[
        'License :: OSI Approved :: GNU License',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    install_requires=[
        'pandas', # Required for data manipulation
        'catboost', # Required for machine learning
        'tqdm', # Required for progress bars
        'liftover', # Required for genome conversion

    ],
    cmdclass={
        'install': CustomInstallCommand,
    },
)
