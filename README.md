##README

This folder contains an R package with data and code sufficient to reproduce the analyses in (the as-yet-unpublished Maehrlab scRNA thymus atlas manuscript).

###Software practices

####File structure

#####Code 

To reproduce our analyses, start in 

- `calls_draft_early_2017.Rmd` is the master script for the paper draft.
- Everything other `.Rmd` is called by the master script to produce a figure or part of a figure.

#####Not code

- `old_scripts`: is purely for convenience. Since each results folder already has a copy of the associated script, you can modify `old_scripts` in any way you like.
- `key_intermediates`: holds stable results that have been discussed with Rene and are not expected to change. 
- `results`: everything important is in here, grouped around themes such as thymus atlas paper draft #5 in `2017_jan_thy_organogenesis_5`. 
- `tables`: gene lists that are necessary for proper functioning of code in `thymus_functions.Rmd`: human-mouse orthologs, receptor-ligand pairs, handpicked genes important in thymus development, et cetera.


####Tools/Methods

#####Workflow

The workflow uses the [`freezr` package](https://github.com/ekernf01/freezr) to save code, results, and session info.

#####Version control

This folder contains hidden files that belong to a version-control system called Git. If you are new to version control, I like the intro [here](https://betterexplained.com/articles/a-visual-guide-to-version-control/). For more information on these files, type `git status`, `git log`, or `git branch`. To bake recent changes in, use `git add .` followed by `git commit`. (Upon committing, a Vim window will appear. Briefly describe your changes, save, and quit.) There is a [backup of this repository](https://bitbucket.org/maehrlab/scrna) on BitBucket; see it using `git remote show origin` and push new changes to it via `git push -u origin --all`.

#####R markdown

I chose to move the code over to R markdown (`.Rmd`) documents whenever possible. R markdown is a set of tools for presenting and annotating code as part of a "living document". It helps me stay sane by splitting my code into bite-size (byte-size?) chunks with formatted comments. The editor RStudio and the R package `knitr` together have machinery to:

- edit the code and insert formatted comments
- run the code
- "knit" the code and documentation into a nice PDF or HTML document
- "purl" the code into a plain R script

The first time I saw R markdown, I had no idea how to parse it. It was much easier to disentangle once I learned about (plain) [markdown](https://daringfireball.net/projects/markdown/). 