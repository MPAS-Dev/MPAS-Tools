.. _dev_releasing:

***********************
Releasing a New Version
***********************

This document describes the steps for maintainers to tag and release a new
version of ``mpas_tools``, and to update the conda-forge feedstock.

Version Bump and Dependency Updates
===================================

1. **Update the Version Number**

   - Open a pull request (PR) to update the version number in the following
     two files:
     - ``conda_package/mpas_tools/__init__.py``
     - ``conda_package/recipe/meta.yaml``

   - Make sure the version follows semantic versioning (see
     https://semver.org/).
     For release candidates, use versions like ``1.3.0rc1`` (no ``v`` prefix).

2. **Check and Update Dependencies**

   - Ensure that dependencies and their constraints are up-to-date and
     consistent in:

     - ``conda_package/recipe/meta.yaml`` (dependencies for the conda-forge
       release)
     - ``conda_package/pyproject.toml`` (dependencies for PyPI; used as a
       sanity check)
     - ``conda_package/dev-spec.txt`` (development dependencies; should be a
       superset of those for the conda-forge release)

   - The dependencies in ``meta.yaml`` are the ones that will be used for the
     released package on conda-forge. The dependencies in ``pyproject.toml``
     are for PyPI and should be kept in sync as much as possible but are only
     there as a sanity check when we run ```pip check``. The ``dev-spec.txt``
     file should include all dependencies needed for development and testing.

   - Review and update dependency versions and constraints as needed.

3. **Make a PR and merge it**

Tagging and Publishing a Release Candidate
==========================================

4. **Tagging a Release Candidate**

   - For release candidates, **do not create a GitHub release page**. Just
     create a tag from the command line:

     - Make sure your changes are merged into ``master`` (or your update
       branch) and your local repo is up to date.

     - Tag the release candidate (e.g., ``1.3.0rc1``):

       ::

           git checkout master
           git fetch --all -p
           git reset --hard origin/master
           git tag 1.3.0rc1
           git push origin 1.3.0rc1

       (Replace ``1.3.0rc1`` with your actual version, and ``master`` with
       your branch if needed.)

     **Note:** This will only create a tag. No release page will be created
     on GitHub.

5. **Updating the conda-forge Feedstock for a Release Candidate**

   - The conda-forge feedstock does **not** update automatically for release
     candidates.
   - You must always create a PR manually, and it must target the ``dev``
     branch of the feedstock.

   Steps:

   - Download the release tarball:

     ::

         wget https://github.com/MPAS-Dev/MPAS-Tools/archive/refs/tags/<version>.tar.gz

   - Compute the SHA256 checksum:

     ::

         shasum -a 256 <version>.tar.gz

   - In the ``meta.yaml`` of the feedstock recipe:
     - Set ``{% set version = "<version>" %}``
     - Set the new ``sha256`` value
     - Update dependencies if needed

   - Commit, push to a new branch, and open a PR **against the ``dev`` branch**
     of the feedstock:
     https://github.com/conda-forge/mpas_tools-feedstock

   - Follow any instructions in the PR template and merge once approved

Publishing a Stable Release
===========================

6. **Publishing a Stable Release**

   - For stable releases, create a GitHub release page as follows:

     - Go to https://github.com/MPAS-Dev/MPAS-Tools/releases

     - Click "Draft a new release"

     - Enter a tag (e.g., ``1.3.0``; **do not** include a ``v`` in front)

     - Set the release title to the version prefixed with ``v`` (e.g.,
       ``v1.3.0``)

     - Generate or manually write release notes

     - Click "Publish release"

7. **Updating the conda-forge Feedstock for a Stable Release**

   - Wait for the ``regro-cf-autotick-bot`` to open a PR at:
     https://github.com/conda-forge/mpas_tools-feedstock

   - This may take several hours to a day.

   - Review the PR:
     - Confirm the version bump and dependency changes
     - Merge once CI checks pass

   **Note:** If you are impatient, you can accelerate this process by creating
   a bot issue at: https://github.com/conda-forge/mpas_tools-feedstock/issues
   with the subject ``@conda-forge-admin, please update version``.  This
   will open a new PR with the version within a few minutes.

   - If the bot PR does not appear or is too slow, you may update manually (see
     the manual steps for release candidates above, but target the ``main``
     branch of the feedstock).

Post Release Actions
====================

8. **Verify and Announce**

   - Install the package in a clean environment to test:

     ::

         conda create -n test-mpas-tools -c conda-forge mpas_tools=<version>

   - Optionally announce the release on relevant communication channels

   - Update any documentation or release notes as needed
