from setuptools import setup, find_packages


def get_requires():
    requires = []
    with open("requirements.txt", "rt") as req_file:
        for line in req_file:
            requires.append(line.rstrip())
    return requires


setup(name='opentargets_pharmgkb',
      version='0.2.0',
      packages=find_packages(),
      package_data={
          'opentargets_pharmgkb': ['OT_SCHEMA_VERSION']
      },
      install_requires=get_requires(),
      tests_require=get_requires(),
      setup_requires=get_requires(),
      test_suite='tests'
      )
