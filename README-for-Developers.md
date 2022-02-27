# README for developers
This file contains some notes for the developers to be aware of.
This is by far not exhaustive or very structured.
See it as a dump for links and information relevant for Julia package development!

## Julia specific conventions
Developing packages in Julia there are certain conventions.
Although not aware of all of them, the following paragraphs outline some resources helpful in learning them.

_Code Style_: https://github.com/invenia/BlueStyle

_Testing_: _Unit tests_ cover single components of the code (developers perspective), while _Integration tests_ (a.k.a. _Functional tests_) cover entire use cases (user perspective).

_Debugging/Coding workflow_: some tutorial to get into the hang of using VSCode Studio for Julia package development:
https://techytok.com/lesson-workflow/

## DevOps

_Git Flow_: Main project is `develop`.
Features are to be developed in branches and merged.
For further information, see: https://discourse.julialang.org/t/good-practices-for-package-development-in-the-julia-ecosystem/8175/2

_Commit messages_: https://www.conventionalcommits.org/en/v1.0.0/

_Version numbering_: https://semver.org

