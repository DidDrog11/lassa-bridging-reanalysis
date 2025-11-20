# r-reproducible-repo

The repository provides a template for reproducible R projects. It uses [`{targets}`](https://books.ropensci.org/targets/) to manage code execution (among other things) and [`{renv}`](https://rstudio.github.io/renv/index.html) to manage packages. Compute environments are managed via [docker](https://rocker-project.org/) and continuous integration happens in [github actions](https://docs.github.com/en/actions).

## How to use the template

1)  On [github.com](https://github.com/viralemergence/r-reproducible-repo) click "Use this template" or run the following code from the command line.

```         
gh repo create my-new-project --template viralemergence/r-reproducible-repo --private --clone
```

2)  [Clone the repo](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) from Github to your local machine

```         
# in location where the project should live
git clone https://github.com/user-name/repo-name.git
```

3)  Open `scripts/project_startup.R` and follow the script.

4)  Write code in {targets}, use `renv::snapshot()` to keep dependencies up to date, and share your project with github!

    -   `usethis::use_r("my-function")` will create a new file in your /R folder
    -   [`{fnmate}`](https://github.com/MilesMcBain/fnmate) is super handy for this as well.

## Includes

-   \_targets.R
-   packages.R
-   .github/workflow
-   .env
-   /R
-   /Scripts/project_startup.R
-   License.md
-   Readme.md

## Nice conventions

-   List all your packages in `packages.R` 
-   All targets are nouns and all functions are verbs
-   Each function lives in its own file.
-   Each function is well documented with roxygen tags (press `cmd + shift + r` with cursor inside the function parenthesis)
-   Use the scripts folder for code that works outside the targets framework and that you want to save (e.g. testing/debugging functions, data wrangling that only has to happen once, etc.)
-   Use psuedo-code/in line comments to make `_targets.R` easier to follow
-   Use an R formatter like [`{air}`](https://www.tidyverse.org/blog/2025/02/air/) to keep your code tidy or use built-in IDE functions e.g. `cmd + shift + i` in rstudio.
-   Commit often. Keep commit messages brief.
-   Use a branching strategy and issues when working in GIT to keep work manageable. Each branch should address a specific issue.
    -   enhancement - branch name used for improving the code base (re-factoring, improving documentation, etc)
    -   fix - branch name used to repair some broken piece of code previously commited to `main`
    -   feature - branch name used for adding new elements to the code base (e.g. new modeling approach, new visualizations, etc)
    -   hotfix - branch name used to repair something that needs immediate attention.

### renv

At some point `renv` will cause you a headache. Its a good idea to talk to collaborators about how you will manage the lock file.

Ideally one person runs `renv::init()` when the project is started and everyone can just use `renv::restore()` to keep their package versions in sync. In practice there is always one package that causes an issue. This becomes more of a problem as R versions diverge among collaborators.

Possible solutions: burn it to the ground run `renv::deactivate(clean = TRUE)` and then re-initialize `renv`. Less drastic measures include using `renv::hydrate()` to load whatever version of the package is already on your machine or `renv::snapshot()` to update the lockfile with new package versions because some package won't compile.

## Running Github Actions locally (e.g. not wasting time/server minutes)

The test-targets.yml file includes conditions for running github actions locally using [ACT](https://nektosact.com/introduction.html).

`act` is a command line tool so all instructions should be run in a terminal.

These instructions are for Mac/Linux machines - windows users should be able to follow along with [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

### Installation

0)  install [homebrew](https://brew.sh/) for package management
1)  install docker `brew install docker` or use [desktop application](https://www.docker.com/products/docker-desktop/)
2)  install act `brew install act`

### The fundamentals

### List workflows

```         
act -l

### returns
# Stage  Job ID        Job name      Workflow name  Workflow file     Events                
# 0      test-targets  test-targets  test-targets   test-targets.yml  push,workflow_dispatch
```

This repo has a single workflow file called `test-targets.yml` that lives in `.github/workflows`.

The `Job ID` comes from line 22 in the yml file and since we did not explicitly set a name, `Job name` is the same as `Job ID`. `Job Name` is what will be displayed on github when the job is run.

Workflow name comes from line 1 and is the name of the workflow that will be displayed on github when the action is run.

`Events` tells us under what conditions the workflow will be run. In this case, on push and via manual dispatch (clicking a button in github).

#### Filtering the list

The following command will list all workflows that are triggered on pull request.

```         
act -l pull_request
```

Since none of the actions in the repo are setup to run after a pull request is initiated, no items are returned.

By default, `act` looks for workflows with `push` triggers.

### run workflows

Note: on an intel mac I had to make some minor changes to two files

```         
### change credsStore to credStore
vim ~/.docker/config.json

### add the following text to avoid errors with mounting docker images
### --container-daemon-socket -
~ % vim ~/.actrc
```

#### Run all workflow files that have an "on push" trigger

```         
act
```

#### run a specific job

```         
act -j test-targets
```

#### run a specific job and re-use containers

Caches containers to make things more efficient.

```         
act -r -j test-targets
```

### Documentation

To see all the options you can set with `act` run

```         
act --help
```

For some worked examples see the[usage guide](https://nektosact.com/usage/index.html).
