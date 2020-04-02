Contributing
============

Contribute to this project according to the following guide:

Fork this repository
--------------------

`Fork this repository before contributing`_. It is a better practice, possibly even enforced, that only Pull Request from forks are accepted, this creates a cleaner representation of the whole `contributions to the project`_.

Install for developers
----------------------

First, clone the repository as described in the :ref:`section above<Fork this repository>`.

Create a dedicated Python environment where to develop the project.

If you are using :code:`pip` follow the official instructions on `Installing packages using pip and virtual environments`_, most likely what you want is:

::

    python3 -m venv idpconfgen
    source idpconfgen/bin/activate

If you are using `Anaconda`_ go for:

::

    conda create --name idpconfgen python=3.7
    conda activate idpconfgen

Either under *pip* or *conda*, install the package in :code:`develop` mode, and also :ref:`tox<Uniformed Tests>`.

::

    python setup.py develop
    # for pip
    pip install tox bumpversion
    # for conda
    conda install tox bumpversion -c conda-forge

Under this configuration the source you edit in the repository git folder is automatically reflected in the development installation.

Continue your implementation following the development guidelines described bellow.

Clone your fork
---------------

Indeed the first thing to do is to clone your fork, and keep it `up to date with the upstream`_:

::

    git clone https://github.com/YOUR-USERNAME/IDPConformerGenerator.git
    cd into/cloned/fork-repo
    git remote add upstream https://github.com/julie-forman-kay-lab/IDPConformerGenerator.git
    git fetch upstream
    git pull upstream develop

New feature
-----------

To work on a new feature, branch out from the ``latest`` branch:

::
    
    git checkout develop
    git checkout -b feature_branch

Develop the feature and keep regular pushes to your fork with comprehensible commit messages.

Push to latest
--------------

To see your development accepted in the main project, you should create a `Pull Request`_ to the ``latest`` branch following the `PULLREQUEST.rst`_ guidelines.

**Before submitting a Pull Request, verify your development branch passes all tests as** :ref:`described bellow<Uniformed Tests>` **. If you are developing new code you should also implement new test cases.**


To official contributors
------------------------

Release Branches
~~~~~~~~~~~~~~~~

Create a short lived branch to prepare for the release candidate, in this example ``release/0.1.0``.

::
    
    git checkout latest
    git checkout -b release/0.1.0

Fix the final bugs, docs and minor corrections, and finally bump the version.

::
    
    # first commit and push your changes
    # then bump
    bumpversion patch|minor|major
    git push origin release/0.1.0

Finally, merge to ``master`` AND from ``master`` to ``latest``.

::
    
    git checkout master
    git merge --no-ff release/0.1.0
    git push origin master --tags
    git checkout latest
    git merge --no-ff master

Hotfixes from master
~~~~~~~~~~~~~~~~~~~~

The hotfix strategy is applied when a bug is identified in the production version that can be easily fixed.

::
    
    git checkout master
    git checkout -b hotfix_branch

Work on the fix...

::

    # push yours commits to GitHub beforehand
    git push origin hotfix_branch  
    bumpversion patch
    git push origin hotfix_branch --tags
    git checkout master
    git merge --no-ff hotfix_branch
    git push origin master
    git checkout latest
    git merge --no-ff master 
    git push origin latest


Uniformed Tests
---------------

Thanks to `Tox`_ we can have an uniform testing platform where all developers are forced to follow the same rules and, above all, all tests occur in a controlled Python environment.

::

    pip install tox tox-conda
    # or
    conda install tox tox-conda -c conda-forge


Before creating a Pull Request from your branch, certify that all the tests pass correctly by running:

::
    
    tox

Also, you can run individual environments if you wish to test only specific functionalities, for example:

::
    
    tox -e check  # code style and file compatibility
    tox -e spell  # spell checks documentation
    tox -e docs  # only builds the documentation
    # to run only pytest
    tox -e py37-nocov




.. _read further about src layout: https://python-project-skeleton.readthedocs.io/en/latest/configuration.html#project-layout
.. _Tox: https://tox.readthedocs.io/en/latest/install.html
.. _Fork this repository before contributing: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/network/members
.. _up to date with the upstream: https://gist.github.com/CristinaSolana/1885435
.. _contributions to the project: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/network
.. _Gitflow Workflow: https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
.. _Pull Request: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/pulls
.. _PULLREQUEST.rst: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/blob/develop/PULLREQUEST.rst
.. _1: https://git-scm.com/docs/git-merge#Documentation/git-merge.txt---no-ff
.. _2: https://stackoverflow.com/questions/9069061/what-is-the-difference-between-git-merge-and-git-merge-no-ff
.. _bumpversion: https://pypi.org/project/bumpversion/
.. _versioneer: https://github.com/warner/python-versioneer
.. _Installing packages using pip and virtual environments: https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment
.. _Anaconda: https://www.anaconda.com/
