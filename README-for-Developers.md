# README for developers
This file contains some notes for the developers to be aware of.
This is by far not exhaustive or very structured.
See it as a dump for links and information relevant for Julia package development!

## Julia specific conventions
Developing packages in Julia there are certain conventions.
Although not aware of all of them, the following paragraphs outline some resources helpful in learning them.

__Code Style__: https://github.com/invenia/BlueStyle

__Testing__: _Unit tests_ cover single components of the code (developers perspective), while _Integration tests_ (a.k.a. _Functional tests_) cover entire use cases (user perspective). See https://eth-vaw-glaciology.github.io/course-101-0250-00/lecture6/#unit_testing_and_reference_tests_in_julia for some hands-on explanations.

__Debugging/Coding workflow__: some tutorial to get into the hang of using VSCode Studio for Julia package development:
https://techytok.com/lesson-workflow/

## DevOps

___Git Flow__: Main project is `develop`.
Features are to be developed in branches and merged.
For further information, see: https://discourse.julialang.org/t/good-practices-for-package-development-in-the-julia-ecosystem/8175/2

__Commit messages__: There are standards to commit messages.
Various info is summarized in the file `.git-commit-msg-template`.
It also serves as a template for commit messages: execute `git config commit.template .git-commit-msg-template` to make it the standard commit message for the repository, which will be shown each time you `git commit`.

__Version numbering__: https://semver.org

__Before pushing__: you **could** run test suite locally and test generation of documentation.
To do so use `julia --project=docs docs/make.jl` or `julia --project=. test/runtests.jl`
(as could be inferred from `.github/workflows/Documentation.yml`)

