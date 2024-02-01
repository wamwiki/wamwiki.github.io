---
title: "Contribution Guidelines"
linkTitle: "Contribution Guidelines"
weight: 100
description: >
  How to contribute to this website
---

{{% pageinfo %}}

We use **Hugo**, an open-source static site generator.

Hugo takes **Markdown** or **HTML** pages, applies a theme, and wraps everything up into a website.

We prefer to write documentation pages in **R-Markdown** so that we can easily combine text, equations, R-code and R-output.
We then knit the **R-Markdown** pages into **Markdown** so that Hugo can use them.

All submissions, including submissions by project members, require review. We
use **GitHub pull requests** for this purpose. Consult
[GitHub Help](https://help.github.com/articles/about-pull-requests/) for more
information on using pull requests.

**If you are not familiar with Git, Github, R-Markdown you can send your contributions by email to <wambamwiki@gmail.com>.**

{{% /pageinfo %}}

## Quick start

Here's a quick guide to updating the Reference code and tutorial pages. It assumes you're familiar with the
GitHub workflow and explain [how to configure the automated preview](#setup-hugo) of your doc updates.


1. Fork the [wamwiki.github.io  repo](https://github.com/wamwiki/wamwiki.github.io) on GitHub.
1. Create a folder for your tutorial located within `content/en/docs/` or a subdirectory within.
1. Add a `_index.Rmd` file, `data` and `images` files as needed within the same folder **not** in subfolders
1. Knit the `.Rmd` to get the `.md` `.html` and plot files
1. Check that everything looks good on a [local site](#check-local) or at least on the `.html` file produced.
1. Create a pull request to get your changes merged into the main branch and update the website.Your PR should include at least the `.Rmd`, `.md` files as well as the `data`, `plots` and `images` files if needed

## Workflow from .Rmd to a working webpage

A brief explanation of all the steps needed to get started or refresh your memory.

### Create a new folder for your tutorial

First question is where?

When building the website `Hugo` is looking for `.md` within the `/content/en/` directory and subdirectories within. All files needed for a page  (plots, images) should be located at the same level as the `.md` file referring to them.  The url to a page is `{site url}/path_within_content_en/name_of_file _without_md_extension`. So if you create a "Unicorn_magic" folder including a `h2.md` file, the page url will be `{site url}/unicorn_magic/h2`. The exception is files "_index.md" are seen as root and their url would be`{site url}/path_within_content_en/`. If your file is `_index.md` then its url will be `{site url}/unicorn_magic`.

In most cases, having a folder with `_index.md` is the easiest for navigation and folder structure.

We thus suggest to create a folder in the appropriate location within `content/en` containing a file `_index.Rmd`.

### Initialise a .Rmd file 

To make your life easier, you can download a [short `.Rmd` file](/docs/Starter.Rmd) with the suggested structure of the frontmatter and suggested knitr options [here](/docs/Starter.Rmd).

Your file frontmatter should include at least the following properties in the YAML

  * title (what appears at the top of the page, within quotation marks)
  * linkTitle (what url to the page look like, within quotation marks)
  * author: your name
  * weight (the higher the weight the lower the page link in side menus, just a number)
  * description (what textual links to the page look like, **no** quotation marks but start with ">" then starts a new line, see example below or in the source Markdown files)
  * math: true/false if you have any equations
  * output: should be html and specify the option to keep the md file. You must insure that a Markdown file is kept after kniting as Hugo works on Markdown, not R-Markdown. In RStudio, it can be done if you go to `Edit the R Markdown format options for the current file` (that is the little cogwheel right of the knit button) -> `Output Options...` -> `Advanced` and check `Keep markdown source file`.

Here is what the frontmatter of your `.Rmd` file should look like

```yaml
---
title: "Name of your QG tutorial"
linkTitle: "Heritability of Magical beast"
author: "Wam community extraordinaire"
weight: 4
math: true
description: >
  Fitting a simple animal model of unicorn color preference.
output: 
  html_document: 
    keep_md: yes
---
```

### Write a working `_index.Rmd` file

  * Reading data. Save your `data` files directly in your working folder (*i.e.* folder where you have your `.Rmd` file). All your call to read data will be of the type `read.csv("mydata.csv")` (assuming you are using `.csv`). Since data files are used only when knitting and not shown on pages, they can be saved in a subfolder.
  * Including images: Save them in your working folder and refer to them using either Markdown notation (`![](cute_image.jpg)`) or with a R code chunk as:

````md
```{r}
knitr::include_graphics("cute_image.jpg")
```
````

  * Including R plots: You need to save plots generated by R somewhere that is accessible by `hugo` . The easiest place is simply your working folder, so include the following R chunk at the beginning of your file just below the frontmatter. It will generate the R plots and save them in your working directory.

````md
```{r setup_fig, include=FALSE}
knitr::opts_chunk$set(fig.path = "")
```
````

### Knit your file

* Make sure that your file knit correctly without any error or warning. Since you will be potentially knitting your file multiple times, I recommend caching the output or your R code chunks so that knitting is faster. 
When option `cache = TRUE` for a code chunk then R code is run once and the output are saved (*cached*). The R code chunks are rerun only when you edit the code chunks or the cached fies are compromised.

You can setup the cache option globally for all code chunks with 

````md
```{r setup_cache, include=FALSE}
knitr::opts_chunk$set(cache = TRUE)
```
````

### Check the output on local site {#check-local}

You need to have a working installation of `Hugo`. If you don't, see [setting-up Hugo](#setup-hugo)

In a terminal, from your fork root directory, run:

```bash
hugo server
```

By default your site will be available at <http://localhost:1313/>. Now that you're serving your site locally, Hugo will watch for changes to the content and automatically refresh your site. Now that you're serving your site locally, Hugo will watch for changes to the content and automatically refresh your site each you save a `.md` file (*i.e.* each time you knit your `.Rmd`).

### Create a Pull Request to github

1. While editing your tutorial, continue with the usual GitHub workflow to edit files, and commit them. You can also push them regularly to your fork
1. Once you are happy with your tutorial/edits, make sure your knit the latest version of your `_index.Rmd` file, then include all latest changes in a commit.
1. Make sure you have pushed to your fork through (one or multiple) commit(s) your `.Rmd` and rendered `.md` file as well as any `data`, `images` and `plots` files as needed.
1. Create a pull request to merge the changes in your fork to the branch `master` on wamwiki


## Updating a single page

If you've just spotted something you'd like to change while using the docs, here is a shortcut for you:

* Click **Edit this page** in the top right hand corner of the page.

You can also edit the page on your fork (please make sure that your fork is up-to-date with the latest commit on the master branch) and make a pull request of course

## Creating an issue

If you've found a problem in the Wamwiki docs, but you're not sure how to fix it yourself, please create an issue in the [wamwiki.github.io  repo](https://github.com/wamwiki/wamwiki.github.io/issues). You can also create an issue about a specific page by clicking the **Create documentation issue** button in the top right hand corner of the page.

## Setting-up Hugo {#setup-hugo}

If you want to run your own local Hugo server to preview your changes as you work you need to have the **extended** version of `Hugo` installed on your computer including its dependencies `git`, `Go` and `Dart Sass`. For information on how to install follow the links below:

* [Hugo](https://gohugo.io/installation/). You need the extended version
* [Go](https://go.dev/doc/install)
* [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
* [Dart Sass](https://gohugo.io/hugo-pipes/transpile-sass-to-css/#dart-sass)

## Useful resources about the Docsy theme, Hugo and Github

* [Docsy user guide](https://www.docsy.dev/docs/): All about Docsy, including how it manages navigation, look and feel, and multi-language support.
* [Hugo documentation](https://gohugo.io/documentation/): Comprehensive reference for Hugo.
* [Github Hello World!](https://guides.github.com/activities/hello-world/): A basic introduction to GitHub concepts and workflow.
