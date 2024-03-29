from setuptools import setup

setup(name='opticalmaterialspy',
      version='0.3.4',
      description='Python library with optical material properties.',
      url='https://github.com/jtambasco/opticalmaterialspy',
      author='Jean-Luc Tambasco',
      author_email='an.obscurity@gmail.com',
      license='MIT',
      install_requires=[
          'numpy',
          'scipy',
          'future'
      ],
      packages=['opticalmaterialspy'],
      zip_safe=False)
