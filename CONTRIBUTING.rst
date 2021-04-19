Contributing
============

How to contribute to this project.

Fork this repository
--------------------

`Fork this repository before contributing`_.

Clone your fork
~~~~~~~~~~~~~~~

Next, clone your fork to your local machine, keep it `up to date with the
upstream`_, and update the online fork with those updates.

::

    git clone https://github.com/YOUR-USERNAME/IDPConformerGenerator.git
    cd IDPConformerGenerator
    git remote add upstream git://github.com/julie-forman-kay-lab/IDPConformerGenerator.git
    git fetch upstream
    git merge upstream/master
    git pull origin master

Install for developers
----------------------

Create a dedicated Python environment where to develop the project.
If you are using `Anaconda`_ go for::

    conda env create -f requirements.yml

Install the package in the repository::

    python setup.py develop --no-deps

This configuration, together with the use of the ``src`` folder layer, guarantee
that you will always run the code after installation. Also, thanks to the
``develop`` flag, any changes in the code will be automatically reflected in the
installed version.

Make a new branch
-----------------

From the ``master`` branch create a new branch where to develop the new code::

    git checkout master
    git checkout -b new_branch


Develop the feature and keep regular pushes to your fork with comprehensible
commit messages.

::

    git status
    git add (the files you want)
    git commit (add a nice commit message)
    git push origin new_branch

While you are developing, you can execute ``tox`` as needed to run your
unittests or inspect lint, etc. ``tox`` usage is described in the last
section of this page.

Update CHANGELOG
~~~~~~~~~~~~~~~~

Update the changelog file under :code:`docs/CHANGELOG.rst` with an explanatory
bullet list of your contribution. Add that list right after the main title and
before the last version subtitle::

    Changelog
    =========

    * here goes my new additions explain them shortly and well

    vX.X.X (1900-01-01)
    -------------------

Also add your name to the authors list at :code:`docs/AUTHORS.rst`.

Pull Request
~~~~~~~~~~~~

Once you are finished, you can Pull Request you additions to the main
repository, and engage with the community. Please read the
`docs/PULLREQUEST.rst`_ guidelines first, you will see them when you open a PR.

**Before submitting a Pull Request, verify your development branch passes all
tests as** :ref:`described bellow<Uniformed Tests with tox>` **. If you are
developing new code you should also implement new test cases.**


Uniformed Tests with tox
------------------------

Thanks to `Tox`_ we can have a unified testing platform where all developers are
forced to follow the same rules and, above all, all tests occur in a controlled
Python environment.

Before creating a Pull Request from your branch, certify that all the tests pass
correctly by running:

::

    tox

These are exactly the same tests that will be performed online in the Github
Actions.

Also, you can run individual environments if you wish to test only specific
functionalities, for example:

::

    tox -e lint  # code style
    tox -e build  # packaging
    tox -e docs  # only builds the documentation
    tox -e prreqs  # special requirements before Pull Request
    tox -e py37  # runs code tests in Python 3.7 environment


.. _Fork this repository before contributing: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/network/members
.. _up to date with the upstream: https://gist.github.com/CristinaSolana/1885435
.. _Anaconda: https://www.anaconda.com/
.. _Tox: https://tox.readthedocs.io/en/latest/
.. _docs/PULLREQUEST.rst: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/blob/master/docs/PULLREQUEST.rst
