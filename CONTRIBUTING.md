# Contributing to Assemblerflow

Thank you for your interest in contributing to Assemblerflow. All kinds of 
contributions are welcome :tada:!

## Issues

Feel free to [submit issues](https://github.com/assemblerflow/assemblerflow/issues)
and enhancement requests.

## Git branch convention

Contributions with new code (not documentation), should follow this standard procedure:

1. Create a new branch for the new feature/bug fix.
2. One the new code is finished and **passes all automated tests**, it will be 
merged into the `dev` branch. This branch is where all the new code lives and 
serves as an incubator stage while field tests are performed to ensure that everything
is working correctly.
3. Merging the `dev` code into `master` is associated with a new release. Therefore, 
the `master` branch is basically the same of the latest official release in PyPI. 

## Contributing

In general, we follow the "fork-and-pull" Git workflow.

 1. **Fork** the repo on GitHub
 2. **Clone** the project to your own machine
 3. **Commit** changes to your own branch
 4. **Push** your work back up to your fork
 5. Submit a **Pull request** so that we can review your changes. Pull requests will be merged first into the `dev` branch to perform some field tests before being merged into `master` 

NOTE: Be sure to merge the latest from "upstream" before making a pull request!
  