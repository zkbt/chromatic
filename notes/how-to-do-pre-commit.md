# how to set up pre-commit

[pre-commit](https://pre-commit.com) seems to be a neat little tool for making sure some tools get run every time someone tries to commit to the repository. Here are some quick notes from Zach trying to get it set up.

I followed the instructions at https://pre-commit.com, including running `pre-commit sample-config` to get a sample start to the config file. I copied and pasted this into `.pre-commit-config.yaml`. I made some modifications based on the example at [exoplanet](https://github.com/exoplanet-dev/exoplanet/blob/main/.pre-commit-config.yaml) to add `black`, apparently both for the code and for the documentation notebooks.

I ran `pre-commit autoupdate` to assign the current latest versions of the tools. It sounds like these won't naturally update, so if there are changes to them, we should rerun `pre-commit autoupdate` at some point to update the version numbers.

I ran `pre-commit install` to connect it to the repository. I iterated `autoupdate` and `install` a few times to get rid of some warnings about mutable versions.

I tried to run `pre-commit run --all-files` to run it on everything that had already been committed. We should redo this if we change the pre-commit rules at any point. The first time I ran it there were some complaints because by the style-checker thingamabobs, but then I put `black` before them, reran, and everything passed happily.

Seems to work!?
