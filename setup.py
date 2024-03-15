from setuptools import setup, find_packages
from setuptools.command.install import install
import boto3

class CustomInstallCommand(install):
    def run(self):
        # Call the default install method
        install.run(self)

        # Download files from S3 after installation
        self.download_files_from_s3()

    def download_files_from_s3(self):
        # Initialize AWS S3 client
        s3 = boto3.client('s3')

        # Specify the S3 bucket name and file key
        bucket_name = 'your-s3-bucket'
        file_key = 'example_file.txt'

        # Specify the local directory to save the downloaded file
        local_directory = 'example_directory'

        # Download file from S3
        s3.download_file(bucket_name, file_key, f'{local_directory}/{file_key}')

setup(
    name='example_package',
    version='1.0.0',
    packages=find_packages(),
    author='Your Name',
    author_email='your.email@example.com',
    description='An example Python package',
    long_description='This is an example Python package with a setup.py file.',
    url='https://github.com/yourusername/example_package',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    install_requires=[
        'boto3',  # Required for interacting with AWS services
        'pandas', # Required for data manipulation
    ],
    cmdclass={
        'install': CustomInstallCommand,
    },
)