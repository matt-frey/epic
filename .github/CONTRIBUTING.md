## How to contribute to EPIC
Hi there. The core developers of EPIC thank you for your contribution. We have some guidelines though which we like you to follow before opening pull request.

<!-- 29 November 2021 -->
<!-- https://stackoverflow.com/questions/11948245/markdown-to-create-pages-and-table-of-contents?page=1&tab=votes#tab-top -->
#### Table of Contents
1. [Unit Testing](#unit-testing)
2. [Coding Style](#coding-style)
2. [Opening a pull request](#pull-request)

### Unit Testing <a name="unit-testing"></a>
Please write unit tests for newly added features and routines. All unit tests are located in the directory `unit-tests`.

### Coding Style <a name="coding-style"></a>
Please follow our coding style when contributing to EPIC.
* use an indentation of 4 spaces
* use lower case letters only
* names of variables and subroutines / functions that consists of multiple words should be separated by underscores, e.g. `this_is_my_subroutine`
* add a whitespace character before and after mathematical operations, e.g. `i = i + 1`
* add a whitespace character after `if` and `while`, e.g. `if (i == 0)`


### Opening a pull request <a name="pull-request"></a>
Developers should fork the main repo and work with the branch `develop`. We only accept pull requests into `develop` of this repository. From time to time we pull all changes in `develop` into `main` and publish a release.