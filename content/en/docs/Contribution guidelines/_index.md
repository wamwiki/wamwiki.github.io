---
title: "Contribution Guidelines"
linkTitle: "Contribution Guidelines"
weight: 100
description: >
  How to contribute to this website
---

{{% pageinfo %}}

{{% /pageinfo %}}

We use Hugo, an open-source static site generator that provides us with templates, 
content organisation in a standard directory structure, and a website generation 
engine. 
We write documentation pages in R-Markdown and knit them to Markdown so that Hugo can use them.
You can also write the pages directly in Markdown, or HTML if you want.
Hugo takes all the Markdown files, applies the theme and wraps everything up into a website.


All submissions, including submissions by project members, require review. We
use GitHub pull requests for this purpose. Consult
[GitHub Help](https://help.github.com/articles/about-pull-requests/) for more
information on using pull requests.

If you are not familiar with Git or Github you can send your contributions by email to one of us.


## Quick start

Here's a quick guide to updating the docs. It assumes you're familiar with the
GitHub workflow and you're happy to use the automated preview of your doc
updates:

1. Fork the [wamwiki.github.io  repo](https://github.com/wamwiki/wamwiki.github.io) on GitHub.
1. To be continued


## Workflow from .Rmd to a working webpage

A brief summary of all the steps needed to get started or refresh your memory.

1. Initialise a .Rmd file with at least the following properties in the YAML
	* title (what appears at the top of the page, within quotation marks)
	* linkTitle (what url to the page look like, within quotation marks)
	* weight (the higher the weidht the lower the page link in side menus, just a number)
	* description (what textual links to the page look like, **no** quotation marks but start with ">" then starts a new line, see example below or in the source Markdown files)

Go to Output 
1. Write a working .Rmd file, that maybe runs models and shows results.
1. Knit the file

## Updating a single page

If you've just spotted something you'd like to change while using the docs, here is a shortcut for you:

1. Click **Edit this page** in the top right hand corner of the page.


## Previewing your changes locally

If you want to run your own local Hugo server to preview your changes as you work:

1. Follow the instructions in [Getting started](/docs/getting-started) to install Hugo and any other tools you need. You'll need at least **Hugo version 0.45** (we recommend using the most recent available version), and it must be the **extended** version, which supports SCSS.
1. Fork the [wamwiki.github.io  repo](https://github.com/wamwiki/wamwiki.github.io) repo into your own project, then create a local copy using `git clone`. Don’t forget to use `--recurse-submodules` or you won’t pull down some of the code you need to generate a working site.

    ```
    git clone --recurse-submodules --depth 1 https://github.com/wamwiki/wamwiki.github.io.git
    ```

1. Run `hugo server` in the site root directory. By default your site will be available at http://localhost:1313/. Now that you're serving your site locally, Hugo will watch for changes to the content and automatically refresh your site.
1. Continue with the usual GitHub workflow to edit files, commit them, push the
  changes up to your fork, and create a pull request.

## Creating an issue

If you've found a problem in the Wamwiki docs, but you're not sure how to fix it yourself, please create an issue in the [wamwiki.github.io  repo](https://github.com/wamwiki/wamwiki.github.io/issues). You can also create an issue about a specific page by clicking the **Create documentation issue** button in the top right hand corner of the page.

## Useful resources about the Docsy theme, Hugo and Github

* [Docsy user guide](https://www.docsy.dev/docs/): All about Docsy, including how it manages navigation, look and feel, and multi-language support.
* [Hugo documentation](https://gohugo.io/documentation/): Comprehensive reference for Hugo.
* [Github Hello World!](https://guides.github.com/activities/hello-world/): A basic introduction to GitHub concepts and workflow.


