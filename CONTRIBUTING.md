# Contributing

Here you will find all the information and steps necessary to contribute to UnCoVar.

## Snakemake

We use `snakemake` as the backbone of this project. Have a look at the [official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

If you already have `conda` installed, run:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

For faster installation, you can use [`mamba`](https://github.com/mamba-org/mamba), a reimplementation of the conda package manager in C++.

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Once installed, activate your `snakemake` environment by running `conda activate snakemake`.

## Pre-commit

Git hook scripts help identify simple issues before submission to code review. We run our hooks on every commit to automatically point out problems in code such as missing semicolons, trailing whitespace, and debug statements. Pointing these issues out before code review allows a code reviewer to focus on the architecture of a change while not wasting time with trivial style nitpicks.

Before you can run hooks, you need to have the [`pre-commit`](https://pre-commit.com/) package manager installed.

Using pip:

```bash
pip install pre-commit
```

Using homebrew:

```bash
brew install pre-commit
```

Using conda (via conda-forge):

```bash
conda install -c conda-forge pre-commit
```

Run pre-commit install to set up the git hook scripts

```bash
$ pre-commit install
pre-commit installed at .git/hooks/pre-commit
```

Now, `pre-commit` will run automatically on `git commit`.

## Get a Copy of the Workflow

Follow the following steps to acquire a local copy of UnCoVar:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
1. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system (to a different place than where you would run analysis).

## Before submitting

Before submitting your code, please do the following:

1. Implement your changes.
1. Add tests, if possible, for the new changes.
1. Edit documentation if you have changed something significant.
1. Make sure the pre-commit hooks had run (`precommit run --all`).
1. [Commit](https://git-scm.com/docs/git-commit) and [push](https://git-scm.com/docs/git-push) your changes to your fork.
1. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.tests`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).

## Other help

You can contribute by spreading the word about this project.
Writing a short article on using this project would also be a considerable contribution.
You can also share your best practices with us.
